#!/usr/bin/env python3
"""
te_autoscan.py - hands-off nucleophile detection and catalytic-His geometry
scanning for AlphaFold3 models of NRPS thioesterase (TE) domain complexes.

What it does, with no per-job input from the user:

  1. Walks an AF3 output tree and finds every predicted model (CIF), pairing
     each with its seed/sample identity and confidence sidecars.
  2. Classifies chains into enzyme (polymer) and substrate (ligand / short
     peptide) without being told which is which.
  3. Locates the Ser-His-Asp catalytic triad geometrically, tolerating a
     Ser->Ala active-site mutation (uses Ala CB as the Ser OG proxy) and
     histidine ring-flip ambiguity.
  4. Perceives every chemically plausible nucleophile on the substrate from
     connectivity, so arbitrary CCD atom naming (C1/N7/O12...) does not matter.
     Amides, guanidines, phenol vs alcohol, thioester sulfur, esters and
     carboxylate oxygens are all distinguished and either scored or excluded.
  5. Measures nucleophile -> His CE1 distance (plus NE2, ND1, flip-invariant
     minimum, elbow atom, and attack distance / Burgi-Dunitz angle to the
     electrophilic carbon) and assigns the near-attack / intermediate / too-far
     bands.
  6. Emits tidy long-format and summary CSVs, a JSON run report, and optional
     PyMOL scripts.

Design notes:
  * Bond perception prefers the authoritative _chem_comp_bond block if the CIF
    carries one; otherwise it falls back to covalent-radius distance criteria.
  * Nothing is keyed on residue names or atom names for the substrate. Enzyme
    side chains do use standard PDB naming, which AF3 always writes correctly.
  * Every heuristic that fires is recorded in the output, so a suspicious row
    can be traced back without rerunning.

Requires: gemmi (pip install gemmi). numpy optional but recommended.
Author: written for the Cryle Lab TE cyclisation pipeline.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import os
import re
import sys
import traceback
from collections import Counter, defaultdict
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Sequence, Set, Tuple

try:
    import gemmi
except ImportError:  # pragma: no cover
    sys.stderr.write(
        "ERROR: gemmi is required.  conda activate af3_ligand && pip install gemmi\n"
    )
    raise

try:
    import numpy as np
except ImportError:  # pragma: no cover
    np = None

__version__ = "4.6.0"

# ---------------------------------------------------------------------------
# Chemistry constants
# ---------------------------------------------------------------------------

# Cordero covalent radii (Angstrom), enough elements for peptides + ligands.
COVALENT_RADII: Dict[str, float] = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "SI": 1.11, "P": 1.07, "S": 1.05, "CL": 1.02, "BR": 1.20, "I": 1.39,
    "SE": 1.20, "AS": 1.19, "FE": 1.32, "ZN": 1.22, "MG": 1.41, "MN": 1.39,
    "CA": 1.76, "NA": 1.66, "K": 2.03, "CU": 1.32, "NI": 1.24, "CO": 1.26,
}
DEFAULT_RADIUS = 1.20
BOND_TOLERANCE = 1.25          # bond if d < (r1 + r2) * tolerance
MIN_BOND_DIST = 0.55

# Bond-length ceilings used to call a C=O when explicit bond orders are absent.
CARBONYL_CO_MAX = 1.32         # C=O (and carboxylate C-O) upper bound
HYDROXYL_CO_MIN = 1.32         # C-OH lower bound
AROMATIC_BOND_RANGE = (1.28, 1.46)

# Maximum heavy-atom valence. Perception that exceeds these caps is wrong, so the
# longest offending bonds are dropped. This is what stops a compressed pose from
# inventing a second C-O bond and destroying a carboxyl group.
VALENCE_MAX = {"C": 4, "N": 4, "O": 2, "S": 3, "P": 5, "B": 4, "SE": 3,
               "F": 1, "CL": 1, "BR": 1, "I": 1}
TETHER_MAX_DIST = 1.95        # ligand-to-enzyme contact treated as a covalent tether

METAL_OR_ION_ELEMENTS = {
    "NA", "K", "MG", "CA", "MN", "FE", "CO", "NI", "CU", "ZN", "CD", "HG",
    "CL", "BR", "I", "F",
}

# Solvent / cryo / buffer junk that should never be treated as the substrate.
NON_SUBSTRATE_RESNAMES = {
    "HOH", "DOD", "WAT", "SO4", "PO4", "GOL", "EDO", "PEG", "DMS", "ACT",
    "ACY", "MES", "TRS", "EPE", "IMD", "CIT", "FMT", "NO3", "IOD", "CL",
    "NA", "MG", "ZN", "CA", "K", "MN", "FE", "FE2", "CU", "NI", "UNX", "UNL",
}

STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE", "SEC", "PYL", "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "CSO",
    "CYX", "ASH", "GLH", "LYN", "TPO", "SEP", "PTR", "MLY", "M3L", "HYP",
}
HIS_ALIASES = {"HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "NEP", "HIC"}
ACID_RESIDUES = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2"),
                 "ASH": ("OD1", "OD2"), "GLH": ("OE1", "OE2")}
# Residue -> atom that occupies the nucleophile-elbow position.
# ALA is included on purpose: a Ser->Ala active-site mutant leaves CB pointing
# where OG was, and that is the standard tethering construct in this project.
ELBOW_ATOMS = {
    "SER": ("OG",), "CYS": ("SG",), "THR": ("OG1",), "ALA": ("CB",),
    "SEP": ("OG",), "CSO": ("SG",), "GLY": ("CA",),
}
ELBOW_PRIORITY = {"SER": 0, "CYS": 1, "THR": 2, "SEP": 0, "CSO": 1,
                  "ALA": 3, "GLY": 6}

# Distance bands for the nucleophile -> His CE1 metric (project convention).
BAND_NEAR_ATTACK = 3.5
BAND_INTERMEDIATE = 6.0


def band_for(distance: Optional[float]) -> str:
    if distance is None:
        return "no_distance"
    if distance < BAND_NEAR_ATTACK:
        return "near_attack"
    if distance <= BAND_INTERMEDIATE:
        return "intermediate"
    return "too_far"


def elem(atom) -> str:
    return atom.element.name.upper()


def radius_of(symbol: str) -> float:
    return COVALENT_RADII.get(symbol.upper(), DEFAULT_RADIUS)


# ---------------------------------------------------------------------------
# Lightweight molecular graph
# ---------------------------------------------------------------------------


@dataclass
class GraphAtom:
    key: Tuple[str, int, str, str]     # (chain, seqid, resname, atom name)
    chain: str
    label_chain: str
    seqid: int
    resname: str
    name: str
    element: str
    pos: "gemmi.Position"
    bfactor: float
    is_polymer: bool
    idx: int = -1

    @property
    def label(self) -> str:
        return f"{self.chain}/{self.resname}{self.seqid}/{self.name}"


class MolGraph:
    """Bond graph over a chosen set of atoms, with simple perception."""

    def __init__(self, atoms: Sequence[GraphAtom]):
        self.atoms: List[GraphAtom] = list(atoms)
        for i, a in enumerate(self.atoms):
            a.idx = i
        self.adj: Dict[int, Set[int]] = defaultdict(set)
        self.order: Dict[Tuple[int, int], str] = {}
        self.aromatic_bond: Dict[Tuple[int, int], bool] = {}
        self.rings: List[List[int]] = []
        self.ring_aromatic: List[bool] = []
        self.aromatic_atoms: Set[int] = set()
        self.bond_source = "none"

    # -- construction --------------------------------------------------
    def add_bond(self, i: int, j: int, order: str = "SING", aromatic: bool = False) -> None:
        if i == j:
            return
        self.adj[i].add(j)
        self.adj[j].add(i)
        k = (min(i, j), max(i, j))
        # Do not let a later SING overwrite a known DOUB.
        if k not in self.order or self.order[k] == "SING":
            self.order[k] = order
        self.aromatic_bond[k] = self.aromatic_bond.get(k, False) or aromatic

    def bond_order(self, i: int, j: int) -> str:
        return self.order.get((min(i, j), max(i, j)), "SING")

    def drop_bond(self, i: int, j: int) -> None:
        self.adj[i].discard(j)
        self.adj[j].discard(i)
        self.order.pop((min(i, j), max(i, j)), None)
        self.aromatic_bond.pop((min(i, j), max(i, j)), None)

    def cap_valence(self) -> List[Tuple[str, str, float]]:
        """Enforce per-element valence limits, dropping the longest bonds first.

        A conformation that squeezes two atoms together can otherwise produce a
        bond that is geometrically plausible but chemically impossible, and one
        spurious C-O bond is enough to make a carboxylic acid look like an ester.
        """
        dropped: List[Tuple[str, str, float]] = []
        for i in range(len(self.atoms)):
            cap = VALENCE_MAX.get(self.atoms[i].element)
            if cap is None:
                continue
            nbrs = [j for j in self.adj[i] if self.atoms[j].element != "H"]
            while len(nbrs) > cap:
                worst = max(nbrs, key=lambda j: self.atoms[i].pos.dist(self.atoms[j].pos))
                d = self.atoms[i].pos.dist(self.atoms[worst].pos)
                dropped.append((self.atoms[i].name, self.atoms[worst].name, round(d, 3)))
                self.drop_bond(i, worst)
                nbrs = [j for j in self.adj[i] if self.atoms[j].element != "H"]
        return dropped

    def build_by_distance(self, tolerance: float = BOND_TOLERANCE,
                          intra_only: bool = True) -> List[Tuple[str, str, float]]:
        """Perceive bonds from interatomic distances.

        With intra_only (the default) a bond is only formed between two atoms of
        the same entity, except for genuine short ligand-to-enzyme contacts which
        are kept so a covalent tether is still detected. Letting enzyme atoms into
        the ligand's own bond graph is what makes perception conformation
        dependent, since a close packing contact reads as a bond.
        """
        n = len(self.atoms)
        for i in range(n):
            ai = self.atoms[i]
            ri = radius_of(ai.element)
            for j in range(i + 1, n):
                aj = self.atoms[j]
                d = ai.pos.dist(aj.pos)
                if d < MIN_BOND_DIST:
                    continue
                if d >= (ri + radius_of(aj.element)) * tolerance:
                    continue
                if ai.element in METAL_OR_ION_ELEMENTS and aj.element in METAL_OR_ION_ELEMENTS:
                    continue
                if intra_only and ai.is_polymer != aj.is_polymer and d > TETHER_MAX_DIST:
                    continue
                self.add_bond(i, j, self._guess_order(ai, aj, d))
        self.bond_source = "distance"
        return self.cap_valence()

    def apply_consensus_bonds(self, pairs: Sequence[Tuple[str, str]],
                              scope: Optional[Set[int]] = None) -> int:
        """Replace perceived intra-entity bonds with a voted, pose-independent set.

        `pairs` are atom-name pairs agreed across the ensemble. Bonds involving
        atoms outside `scope` (the substrate) are left alone, so tether detection
        still works.
        """
        scope = scope if scope is not None else set(range(len(self.atoms)))
        byname: Dict[str, int] = {}
        for i in scope:
            byname[self.atoms[i].name] = i
        for (i, j) in [(i, j) for i in scope for j in list(self.adj[i]) if j in scope]:
            self.drop_bond(i, j)
        applied = 0
        for a, b in pairs:
            if a in byname and b in byname:
                i, j = byname[a], byname[b]
                d = self.atoms[i].pos.dist(self.atoms[j].pos)
                self.add_bond(i, j, self._guess_order(self.atoms[i], self.atoms[j], d))
                applied += 1
        self.bond_source = "ensemble_consensus"
        return applied

    def bond_name_pairs(self, scope: Set[int]) -> List[Tuple[str, str]]:
        out = []
        for i in scope:
            for j in self.adj[i]:
                if j in scope and i < j:
                    a, b = self.atoms[i].name, self.atoms[j].name
                    out.append((a, b) if a <= b else (b, a))
        return out

    @staticmethod
    def _guess_order(ai: GraphAtom, aj: GraphAtom, d: float) -> str:
        pair = {ai.element, aj.element}
        if pair == {"C", "O"} and d <= CARBONYL_CO_MAX:
            return "DOUB"
        if pair == {"C", "N"} and d <= 1.28:
            return "DOUB"
        if pair == {"C", "S"} and d <= 1.68:
            return "DOUB"
        return "SING"

    def apply_chem_comp_bonds(self, comp_bonds: Dict[str, List[Tuple[str, str, str, bool]]]) -> int:
        """Overlay authoritative intra-residue connectivity from _chem_comp_bond."""
        by_res: Dict[Tuple[str, int, str], Dict[str, int]] = defaultdict(dict)
        for a in self.atoms:
            by_res[(a.chain, a.seqid, a.resname)][a.name] = a.idx
        applied = 0
        for (chain, seqid, resname), lookup in by_res.items():
            for a1, a2, order, arom in comp_bonds.get(resname.upper(), []):
                if a1 in lookup and a2 in lookup:
                    self.add_bond(lookup[a1], lookup[a2], order, arom)
                    applied += 1
        if applied:
            self.bond_source = ("chem_comp+distance" if self.bond_source == "distance"
                                else "chem_comp")
        return applied

    # -- perception ----------------------------------------------------
    def neighbours(self, i: int) -> List[int]:
        return sorted(self.adj[i])

    def heavy_degree(self, i: int) -> int:
        return len([j for j in self.adj[i] if self.atoms[j].element != "H"])

    def find_rings(self, max_size: int = 7) -> None:
        """Enumerate small cycles by bounded BFS from each atom."""
        seen: Set[frozenset] = set()
        rings: List[List[int]] = []
        for start in range(len(self.atoms)):
            paths = [[start]]
            while paths:
                nxt = []
                for path in paths:
                    for j in self.adj[path[-1]]:
                        if self.atoms[j].element == "H":
                            continue
                        if j == start and len(path) >= 3:
                            fs = frozenset(path)
                            if fs not in seen:
                                seen.add(fs)
                                rings.append(list(path))
                            continue
                        if j in path or j < start:
                            continue
                        if len(path) < max_size:
                            nxt.append(path + [j])
                paths = nxt
        self.rings = rings
        self._flag_aromatics()

    def _plane_rms(self, idxs: Sequence[int]) -> float:
        pts = [(self.atoms[i].pos.x, self.atoms[i].pos.y, self.atoms[i].pos.z) for i in idxs]
        n = len(pts)
        cx = sum(p[0] for p in pts) / n
        cy = sum(p[1] for p in pts) / n
        cz = sum(p[2] for p in pts) / n
        shifted = [(p[0] - cx, p[1] - cy, p[2] - cz) for p in pts]
        if np is not None:
            m = np.array(shifted)
            try:
                _, _, vt = np.linalg.svd(m, full_matrices=False)
                normal = vt[-1]
                dev = m @ normal
                return float(math.sqrt(float((dev ** 2).mean())))
            except Exception:
                return 0.0
        # numpy-free fallback: average of triangle-normal projections.
        best = 0.0
        for k in range(2, n):
            ax, ay, az = shifted[0]
            bx, by, bz = shifted[1]
            nx = ay * bz - az * by
            ny = az * bx - ax * bz
            nz = ax * by - ay * bx
            norm = math.sqrt(nx * nx + ny * ny + nz * nz)
            if norm < 1e-6:
                continue
            nx, ny, nz = nx / norm, ny / norm, nz / norm
            dev = [abs(p[0] * nx + p[1] * ny + p[2] * nz) for p in shifted]
            best = max(best, sum(dev) / len(dev))
            break
        return best

    def ring_fingerprint(self, ring: Sequence[int]) -> str:
        """Topological identity of a ring: size, elements and heavy degrees.

        Canonical over rotation and reflection, and independent of coordinates and
        of atom naming, so the same ring in two different ligands gets the same
        key and can therefore be given the same aromaticity verdict.
        """
        toks = [(self.atoms[i].element, self.heavy_degree(i)) for i in ring]
        best: Optional[Tuple] = None
        for seq in (toks, toks[::-1]):
            for k in range(len(seq)):
                cand = tuple(seq[k:] + seq[:k])
                if best is None or cand < best:
                    best = cand
        return f"{len(ring)}|" + ",".join(f"{e}{d}" for e, d in (best or ()))

    def aromaticity_votes(self) -> List[Tuple[str, bool]]:
        return [(self.ring_fingerprint(r), a)
                for r, a in zip(self.rings, self.ring_aromatic)]

    def apply_consensus_aromaticity(self, verdicts: Dict[str, bool]) -> int:
        """Replace geometric aromaticity with the run-wide verdict for each ring."""
        changed = 0
        for k, ring in enumerate(self.rings):
            v = verdicts.get(self.ring_fingerprint(ring))
            if v is not None and v != self.ring_aromatic[k]:
                self.ring_aromatic[k] = bool(v)
                changed += 1
        self.aromatic_atoms = set()
        for ring, a in zip(self.rings, self.ring_aromatic):
            if a:
                self.aromatic_atoms.update(ring)
        return changed

    def _flag_aromatics(self) -> None:
        self.ring_aromatic = [False] * len(self.rings)
        for ring in self.rings:
            if not 4 < len(ring) < 7:
                continue
            if any(self.atoms[i].element not in {"C", "N", "O", "S"} for i in ring):
                continue
            flagged = [self.aromatic_bond.get((min(a, b), max(a, b)), False)
                       for a, b in zip(ring, ring[1:] + ring[:1])]
            if all(flagged):
                self.ring_aromatic[self.rings.index(ring)] = True
                self.aromatic_atoms.update(ring)
                continue
            lengths = [self.atoms[a].pos.dist(self.atoms[b].pos)
                       for a, b in zip(ring, ring[1:] + ring[:1])]
            planar = self._plane_rms(ring) < 0.12
            conjugated = all(AROMATIC_BOND_RANGE[0] <= L <= AROMATIC_BOND_RANGE[1]
                             for L in lengths)
            if planar and conjugated:
                self.ring_aromatic[self.rings.index(ring)] = True
                self.aromatic_atoms.update(ring)

    def is_aromatic(self, i: int) -> bool:
        return i in self.aromatic_atoms

    def is_carbonyl_carbon(self, i: int) -> Tuple[bool, Optional[int]]:
        """True if atom i is a C bearing a double-bonded O (or short C-O)."""
        if self.atoms[i].element != "C":
            return False, None
        for j in self.adj[i]:
            if self.atoms[j].element != "O":
                continue
            if (self.bond_order(i, j) == "DOUB"
                    or self.atoms[i].pos.dist(self.atoms[j].pos) <= CARBONYL_CO_MAX):
                return True, j
        return False, None

    def oxygen_count(self, i: int) -> int:
        return len([j for j in self.adj[i] if self.atoms[j].element == "O"])

    def nitrogen_count(self, i: int) -> int:
        return len([j for j in self.adj[i] if self.atoms[j].element == "N"])


# ---------------------------------------------------------------------------
# Nucleophile perception
# ---------------------------------------------------------------------------

# Lower priority number = more likely to be the operative nucleophile in a
# TE-mediated macrocyclisation.  Priority is reported and used only as a
# tie-break; ranking is by distance unless --rank-by priority is given.
CLASS_PRIORITY = {
    "primary_amine": 0,
    "thiol": 1,
    "aliphatic_hydroxyl": 2,
    "secondary_amine": 3,
    "phenol": 4,
    "aryl_amine": 5,
    "tertiary_amine": 6,
    "carboxyl_OH": 7,
    "heteroaromatic_N": 8,
}
# Coarse grouping used for series comparison. The fine distinction between a
# phenol and an aliphatic hydroxyl depends on aromaticity perception, which is the
# least certain judgement the classifier makes, so the coarse group is the safer
# key when comparing one substrate against another.
CLASS_GROUP_OF = {
    "primary_amine": "amine", "secondary_amine": "amine",
    "aryl_amine": "amine", "tertiary_amine": "amine",
    "phenol": "hydroxyl", "aliphatic_hydroxyl": "hydroxyl",
    "thiol": "thiol", "carboxyl_OH": "carboxyl",
    "heteroaromatic_N": "aromatic_n",
}

CLASS_GROUPS = {
    "amine": {"primary_amine", "secondary_amine", "aryl_amine", "tertiary_amine"},
    "hydroxyl": {"aliphatic_hydroxyl", "phenol"},
    "thiol": {"thiol"},
    "carboxyl": {"carboxyl_OH"},
    "aromatic_n": {"heteroaromatic_N"},
}
DEFAULT_CLASSES = {"primary_amine", "secondary_amine", "aryl_amine",
                   "aliphatic_hydroxyl", "phenol", "thiol"}


@dataclass
class NucCandidate:
    atom_idx: int
    chain: str
    label_chain: str
    seqid: int
    resname: str
    atom_name: str
    element: str
    nuc_class: str
    accepted: bool
    reason: str
    plddt: float
    priority: int = 99
    neighbours: str = ""
    nuc_residue: object = ""
    nuc_site: str = ""
    nuc_residue_c: object = ""
    class_group: str = ""
    class_source: str = "pose"
    class_agreement: Optional[float] = None
    pose_class: str = ""
    presence: Optional[float] = None
    # geometry, filled later
    d_ce1: Optional[float] = None
    d_ne2: Optional[float] = None
    d_nd1: Optional[float] = None
    d_imid_min: Optional[float] = None
    d_elbow: Optional[float] = None
    d_attack: Optional[float] = None
    bd_angle: Optional[float] = None
    band: str = "no_distance"


def classify_nucleophile(g: MolGraph, i: int) -> Tuple[str, bool, str]:
    """Return (class, accepted, reason) for a candidate heteroatom."""
    a = g.atoms[i]
    e = a.element
    nbrs = [j for j in g.neighbours(i) if g.atoms[j].element != "H"]
    deg = len(nbrs)

    # Any covalent link to the enzyme means this atom is the tether, not a
    # free nucleophile.
    if any(g.atoms[j].is_polymer != a.is_polymer for j in nbrs):
        return ("tethered", False, "covalently linked to the other entity (tether)")

    if e == "N":
        carbonyl_nbrs = [j for j in nbrs if g.is_carbonyl_carbon(j)[0]]
        # guanidine / amidine: N on a carbon carrying two or more nitrogens
        for j in nbrs:
            if g.atoms[j].element == "C" and g.nitrogen_count(j) >= 2:
                return ("guanidine_N", False, "guanidinium/amidine nitrogen (Arg-like)")
        if g.is_aromatic(i):
            return ("heteroaromatic_N", False, "aromatic ring nitrogen")
        if deg == 0:
            return ("isolated_N", False, "no detected bonds")
        if carbonyl_nbrs:
            if deg == 1:
                return ("amide_N", False, "primary amide NH2 (Asn/Gln-like or C-term amide)")
            if deg == 2:
                return ("amide_N", False, "secondary amide (peptide bond)")
            return ("amide_N", False, "tertiary amide (N-alkylated peptide bond)")
        if deg == 1:
            if g.is_aromatic(nbrs[0]):
                return ("aryl_amine", True, "aniline-type NH2 on an aromatic carbon")
            return ("primary_amine", True, "free primary amine NH2")
        if deg == 2:
            if any(g.is_aromatic(j) for j in nbrs):
                return ("aryl_amine", True, "secondary aryl amine")
            return ("secondary_amine", True, "secondary amine NH (N-alkyl / Pro-type)")
        return ("tertiary_amine", True, "tertiary amine, sterically hindered")

    if e == "O":
        if deg == 0:
            return ("isolated_O", False, "no detected bonds")
        if any(g.atoms[j].element == "P" for j in nbrs):
            return ("phosphate_O", False, "phosphate/phosphonate oxygen")
        if deg >= 2:
            heavy = [g.atoms[j].element for j in nbrs]
            if heavy.count("C") >= 2:
                return ("ester_or_ether_O", False, "bridging ester/ether oxygen")
            return ("bridging_O", False, "oxygen with two heavy neighbours")
        j = nbrs[0]
        if g.atoms[j].element != "C":
            return ("other_O", False, f"oxygen bonded to {g.atoms[j].element}")
        d = a.pos.dist(g.atoms[j].pos)
        is_double = g.bond_order(i, j) == "DOUB" or d <= CARBONYL_CO_MAX
        n_ox = g.oxygen_count(j)
        if is_double:
            if n_ox >= 2:
                return ("carboxyl_C_O", False, "carbonyl oxygen of a carboxylate/ester")
            return ("carbonyl_O", False, "ketone/amide carbonyl oxygen")
        # single-bonded O-H
        if n_ox >= 2:
            return ("carboxyl_OH", True, "carboxylic acid hydroxyl (anhydride chemistry only)")
        if g.is_aromatic(j):
            return ("phenol", True, "phenolic hydroxyl")
        return ("aliphatic_hydroxyl", True, "aliphatic/secondary hydroxyl")

    if e == "S":
        if deg == 0:
            return ("isolated_S", False, "no detected bonds")
        if deg == 1:
            j = nbrs[0]
            if g.is_carbonyl_carbon(j)[0]:
                return ("thioacid_S", False, "thiocarbonyl sulfur")
            return ("thiol", True, "free thiol SH")
        if any(g.atoms[j].element == "S" for j in nbrs):
            return ("disulfide_S", False, "disulfide sulfur")
        if any(g.is_carbonyl_carbon(j)[0] for j in nbrs):
            return ("thioester_S", False, "thioester sulfur (this is the electrophile side)")
        return ("thioether_S", False, "thioether sulfur")

    return (f"non_nucleophilic_{e}", False, "element is not N/O/S")


def enumerate_nucleophiles(g: MolGraph, substrate_idxs: Sequence[int],
                           chain: Optional["PeptideChain"] = None) -> List[NucCandidate]:
    out: List[NucCandidate] = []
    for i in substrate_idxs:
        a = g.atoms[i]
        if a.element not in {"N", "O", "S"}:
            continue
        nuc_class, accepted, reason = classify_nucleophile(g, i)
        nb = ",".join(
            f"{g.atoms[j].name}({g.atoms[j].element})" for j in g.neighbours(i)
            if g.atoms[j].element != "H"
        )
        pos, site = (chain.atom_site.get(i, ("", "unassigned"))
                     if chain is not None else ("", ""))
        pos_c = chain.atom_pos_c.get(i, "") if chain is not None else ""
        out.append(NucCandidate(
            atom_idx=i, chain=a.chain, label_chain=a.label_chain, seqid=a.seqid,
            resname=a.resname, atom_name=a.name, element=a.element,
            nuc_class=nuc_class, accepted=accepted, reason=reason,
            plddt=a.bfactor, priority=CLASS_PRIORITY.get(nuc_class, 99),
            neighbours=nb, nuc_residue=pos, nuc_site=site, nuc_residue_c=pos_c,
            class_group=CLASS_GROUP_OF.get(nuc_class, ""),
        ))
    return out



# ---------------------------------------------------------------------------
# Peptide backbone traversal (no residue names, no atom names, no templates)
# ---------------------------------------------------------------------------


@dataclass
class LigandResidue:
    idx: int
    atoms: Set[int] = field(default_factory=set)
    in_link: Optional[int] = None      # heteroatom bonded to the previous residue
    out_link: Optional[int] = None     # carbonyl carbon bonded to the next residue


@dataclass
class PeptideChain:
    n_units: int = 0
    n_fragments: int = 0
    coverage: float = 0.0
    cyclic: bool = False
    n_ester_links: int = 0
    n_branches: int = 0
    c_term_acyl: Optional[int] = None
    c_term_o: Optional[int] = None
    n_term_amine: Optional[int] = None
    atom_site: Dict[int, Tuple[object, str]] = field(default_factory=dict)
    atom_pos_c: Dict[int, object] = field(default_factory=dict)
    note: str = ""


def _linkage_bonds(g: MolGraph, S: Set[int]) -> List[Tuple[int, int, bool]]:
    """Amide and ester bonds that join one residue to the next.

    A linkage is a carbonyl carbon bonded to an N or O that itself continues to
    another carbon. A terminal NH2 or a free acid OH has nowhere to continue, so
    it is not a linkage, which is what keeps the chain ends intact.
    """
    out: List[Tuple[int, int, bool]] = []
    for c in S:
        is_co, o = g.is_carbonyl_carbon(c)
        if not is_co:
            continue
        for x in g.neighbours(c):
            if x not in S or x == o or g.atoms[x].element not in ("N", "O"):
                continue
            heavy = [j for j in g.neighbours(x)
                     if j in S and g.atoms[j].element != "H"]
            if len(heavy) < 2:
                continue
            if not any(g.atoms[j].element == "C" and j != c for j in heavy):
                continue
            out.append((c, x, g.atoms[x].element == "O"))
    return out


def _components(g: MolGraph, S: Set[int], cut: Set[Tuple[int, int]]) -> List[Set[int]]:
    seen: Set[int] = set()
    comps: List[Set[int]] = []
    for start in sorted(S):
        if start in seen:
            continue
        comp: Set[int] = set()
        stack = [start]
        while stack:
            a = stack.pop()
            if a in comp:
                continue
            comp.add(a)
            for b in g.neighbours(a):
                if b in S and b not in comp and (min(a, b), max(a, b)) not in cut:
                    stack.append(b)
        seen |= comp
        comps.append(comp)
    return comps


def _path_between(g: MolGraph, scope: Set[int], a: int, b: int) -> List[int]:
    if a is None or b is None or a not in scope or b not in scope:
        return []
    prev = {a: None}
    q = [a]
    while q:
        nxt = []
        for u in q:
            if u == b:
                path = []
                while u is not None:
                    path.append(u)
                    u = prev[u]
                return path[::-1]
            for v in g.neighbours(u):
                if v in scope and v not in prev:
                    prev[v] = u
                    nxt.append(v)
        q = nxt
    return []


def trace_peptide(g: MolGraph, sub_idxs: Sequence[int]) -> PeptideChain:
    """Decompose the ligand into residues by cutting its amide and ester linkages.

    Walking assumed N-CA-C(=O) units breaks on anything unusual: a beta or gamma
    residue, a non-amino-acid linker, a branch point. Cutting at the linkage bonds
    instead makes no assumption about what sits between two linkages, so every atom
    lands in exactly one residue and coverage is complete by construction.

    Residues are numbered from both ends. The C-terminal index is the useful one
    for a truncation or substitution series: it is anchored on the acyl carbon that
    the TE domain attacks, so shortening the peptide does not renumber the residues
    that remain.
    """
    S = set(sub_idxs)
    if not S:
        return PeptideChain(note="no substrate atoms")
    links = _linkage_bonds(g, S)
    cut = {(min(c, x), max(c, x)) for c, x, _ in links}
    comps = _components(g, S, cut)
    residues = [LigandResidue(idx=i, atoms=c) for i, c in enumerate(comps)]
    owner: Dict[int, int] = {a: r.idx for r in residues for a in r.atoms}

    succ: Dict[int, List[int]] = defaultdict(list)
    pred: Dict[int, List[int]] = defaultdict(list)
    for c, x, _ in links:
        i, j = owner[c], owner[x]
        if i == j:
            continue
        succ[i].append(j)
        pred[j].append(i)
        residues[i].out_link = c
        residues[j].in_link = x

    def longest(i: int, blocked: Set[int]) -> List[int]:
        best = [i]
        for j in succ[i]:
            if j in blocked:
                continue
            cand = [i] + longest(j, blocked | {i})
            if len(cand) > len(best):
                best = cand
        return best

    starts = [r.idx for r in residues if not pred[r.idx]]
    cyclic = bool(links) and not starts
    main = max((longest(i, set()) for i in (starts or [r.idx for r in residues])),
               key=len)

    # How many separate chains the ligand falls into. One is normal; more than one
    # means a linkage the tracer could not recognise, which is worth a warning
    # because it makes residue numbering incomparable between variants.
    undirected: Dict[int, Set[int]] = defaultdict(set)
    for i, js in succ.items():
        for j in js:
            undirected[i].add(j)
            undirected[j].add(i)
    seen_r: Set[int] = set()
    n_chains = 0
    for r in residues:
        if r.idx in seen_r:
            continue
        n_chains += 1
        stack = [r.idx]
        while stack:
            u = stack.pop()
            if u in seen_r:
                continue
            seen_r.add(u)
            stack.extend(undirected[u] - seen_r)
    chain = PeptideChain(
        n_units=len(residues), n_fragments=n_chains, cyclic=cyclic,
        n_ester_links=sum(1 for _, _, e in links if e),
        n_branches=len(residues) - len(main))

    # C-terminal acyl carbon: a carbonyl carbon in the last residue that is not
    # itself a linkage, i.e. it still carries a leaving group.
    last = residues[main[-1]]
    best_c: Optional[Tuple[int, int, Optional[int]]] = None
    link_cs = {c for c, _, _ in links}
    for a in sorted(last.atoms):
        if g.atoms[a].element != "C" or a in link_cs:
            continue
        is_co, o = g.is_carbonyl_carbon(a)
        if not is_co:
            continue
        n_o = len([j for j in g.neighbours(a)
                   if g.atoms[j].element == "O" and j in S])
        has_s = any(g.atoms[j].element == "S" for j in g.neighbours(a))
        rank = 0 if (n_o >= 2 or has_s) else 1
        if best_c is None or rank < best_c[0]:
            best_c = (rank, a, o)
    if best_c is not None and not cyclic:
        chain.c_term_acyl, chain.c_term_o = best_c[1], best_c[2]
    elif cyclic:
        chain.note = "ligand traced as cyclic: no free C-terminus"

    first = residues[main[0]]
    for a in sorted(first.atoms):
        if g.atoms[a].element != "N":
            continue
        heavy = [j for j in g.neighbours(a) if g.atoms[j].element != "H"]
        if len(heavy) <= 2 and not any(g.is_carbonyl_carbon(j)[0] for j in heavy):
            chain.n_term_amine = a
            break

    L = len(main)
    pos_n = {r: k + 1 for k, r in enumerate(main)}
    pos_c = {r: L - k for k, r in enumerate(main)}
    for r in residues:
        if r.idx in pos_n:
            continue
        # a branch hangs off whichever main-chain residue it links to
        anchor = next((j for j in pred[r.idx] if j in pos_n),
                      next((j for j in succ[r.idx] if j in pos_n), None))
        pos_n[r.idx] = f"{pos_n[anchor]}b" if anchor is not None else "b"
        pos_c[r.idx] = f"{pos_c[anchor]}b" if anchor is not None else "b"

    for r in residues:
        ends = [x for x in (r.in_link, r.out_link) if x is not None]
        if r.idx == main[-1] and chain.c_term_acyl is not None:
            ends.append(chain.c_term_acyl)
        if r.idx == main[0] and chain.n_term_amine is not None:
            ends.append(chain.n_term_amine)
        spine: Set[int] = set()
        for u, v in zip(ends, ends[1:]):
            spine |= set(_path_between(g, r.atoms, u, v))
        if not spine:
            spine = set(ends)
        for a in list(spine):
            is_co, o = g.is_carbonyl_carbon(a)
            if is_co and o is not None and o in r.atoms:
                spine.add(o)
        for a in r.atoms:
            chain.atom_site[a] = (pos_n[r.idx],
                                  "backbone" if a in spine else "side_chain")
            chain.atom_pos_c[a] = pos_c[r.idx]
    if chain.n_term_amine is not None:
        chain.atom_site[chain.n_term_amine] = (pos_n[main[0]], "n_terminal_amine")
    if chain.c_term_acyl is not None:
        chain.atom_site[chain.c_term_acyl] = (pos_n[main[-1]], "c_terminal_acyl")
    chain.coverage = round(len(chain.atom_site) / max(1, len(S)), 3)
    notes = []
    if chain.n_branches:
        notes.append(f"{chain.n_branches} residue(s) off the main chain, "
                     "positions suffixed b")
    if n_chains > 1:
        notes.append(f"ligand falls into {n_chains} separate chains; a linkage was "
                     "not recognised and residue numbering will not match other variants")
    if notes:
        chain.note = "; ".join(notes)
    return chain


def find_electrophile(g: MolGraph, substrate_idxs: Sequence[int],
                      reference: Optional["gemmi.Position"],
                      chain: Optional[PeptideChain] = None,
                      forced_atom: Optional[str] = None
                      ) -> Tuple[Optional[int], Optional[int], str]:
    """Locate the electrophilic (acyl) carbon and its carbonyl oxygen.

    Order: real thioester -> carbon tethered to the enzyme -> C-terminal acyl
    carbon from the traced backbone -> free carboxyl carbon nearest the active
    site. The backbone route is the reliable one: it identifies the acyl carbon
    from its position in the peptide, so it does not care whether the leaving
    group is a thioester, an acid, an ester or an amide.
    """
    sset = set(substrate_idxs)
    # 1. thioester: S bonded to a carbonyl carbon
    for i in substrate_idxs:
        if g.atoms[i].element != "S":
            continue
        for j in g.neighbours(i):
            is_co, o = g.is_carbonyl_carbon(j)
            if is_co:
                return j, o, "thioester_carbonyl"
    # 2. acyl carbon covalently attached to the enzyme (tethered construct)
    for i in substrate_idxs:
        if g.atoms[i].element != "C":
            continue
        linked = [j for j in g.neighbours(i) if g.atoms[j].is_polymer]
        if linked:
            is_co, o = g.is_carbonyl_carbon(i)
            if is_co:
                return i, o, "enzyme_tethered_acyl_carbon"
    # 3. score every acyl candidate. Position in the peptide dominates, but only
    #    if the terminal carbon actually carries a leaving group: an amide carbonyl
    #    with nothing but a continuing peptide nitrogen is not an acyl donor, it is
    #    a sign the trace stopped early.
    term_ok = False
    if chain is not None and chain.c_term_acyl is not None:
        t = chain.c_term_acyl
        hetero = [j for j in g.neighbours(t) if g.atoms[j].element in ("O", "S")]
        n_o = len([j for j in hetero if g.atoms[j].element == "O" and j in sset])
        has_s = any(g.atoms[j].element == "S" for j in hetero)
        linked = any(g.atoms[j].is_polymer for j in g.neighbours(t))
        term_ok = (n_o >= 2) or has_s or linked
    scored: List[Tuple[float, float, int, Optional[int], str]] = []
    for i in substrate_idxs:
        if g.atoms[i].element != "C":
            continue
        is_co, o = g.is_carbonyl_carbon(i)
        oxys = [j for j in g.neighbours(i) if g.atoms[j].element == "O" and j in sset]
        n_ox = len(oxys)
        # Count only intra-substrate carbons: a close packing contact with the
        # enzyme must never make a carboxylic acid read as an ester.
        bridging = [j for j in oxys
                    if len([k for k in g.neighbours(j)
                            if g.atoms[k].element == "C" and k in sset]) > 1]
        at_term = (chain is not None and term_ok
                   and chain.atom_site.get(i, (None, ""))[1] == "c_terminal_acyl")
        score, kind = 0.0, ""
        if n_ox >= 2 and not bridging:
            score, kind = 2.0, "free_carboxyl_C"
        elif n_ox >= 2 and bridging:
            score, kind = 1.5, "ester_carbonyl_C"
        elif is_co and at_term:
            score, kind = 1.0, "peptide_c_terminal_acyl"
        if not score:
            continue
        if at_term:
            score += 3.0
            kind += "_at_backbone_terminus"
        d = g.atoms[i].pos.dist(reference) if reference is not None else 0.0
        scored.append((-score, d, i, o, kind))
    if forced_atom:
        # The ensemble already voted on which carbon is the acyl centre. Keep the
        # chemical label from this pose so the reason is still visible.
        hit = next((t for t in scored if g.atoms[t[2]].name == forced_atom), None)
        if hit is not None:
            _, _, i, o, kind = hit
            if o is None:
                _, o = g.is_carbonyl_carbon(i)
            return i, o, kind + "_consensus"
        for i in substrate_idxs:
            if g.atoms[i].name == forced_atom:
                _, o = g.is_carbonyl_carbon(i)
                if o is None and chain is not None and i == chain.c_term_acyl:
                    o = chain.c_term_o
                return i, o, "consensus_acyl_C_not_scored_in_this_pose"
    if scored:
        scored.sort()
        _, _, i, o, kind = scored[0]
        if o is None:
            _, o = g.is_carbonyl_carbon(i)
        if len(scored) > 1 and scored[1][0] == scored[0][0]:
            kind += "_ambiguous"
        return i, o, kind
    if chain is not None and chain.cyclic:
        return None, None, "substrate_already_cyclic"
    return None, None, "not_found"


# ---------------------------------------------------------------------------
# Catalytic triad detection
# ---------------------------------------------------------------------------


@dataclass
class Triad:
    found: bool
    his_chain: str = ""
    his_label_chain: str = ""
    his_seqid: int = -1
    his_resname: str = ""
    elbow_chain: str = ""
    elbow_seqid: int = -1
    elbow_resname: str = ""
    elbow_atom: str = ""
    acid_chain: str = ""
    acid_seqid: int = -1
    acid_resname: str = ""
    acid_atom: str = ""
    d_his_elbow: Optional[float] = None
    d_his_acid: Optional[float] = None
    orientation: str = ""          # canonical (NE2->elbow) or flipped
    motif: str = ""                # 5-residue window around the elbow
    motif_ok: bool = False
    mutant: bool = False
    cutoff_used: Optional[float] = None
    method: str = ""
    topology_ok: bool = True
    residue_key: str = ""
    forced: bool = False
    consensus_frac: Optional[float] = None
    plddt_mean: Optional[float] = None
    n_candidates: int = 0
    warnings: List[str] = field(default_factory=list)
    # atom handles (not serialised)
    ce1: Optional["gemmi.Position"] = None
    ne2: Optional["gemmi.Position"] = None
    nd1: Optional["gemmi.Position"] = None
    elbow_pos: Optional["gemmi.Position"] = None

    def to_row(self) -> Dict[str, object]:
        d = {k: v for k, v in asdict(self).items()
             if k not in {"ce1", "ne2", "nd1", "elbow_pos", "warnings"}}
        d["triad_warnings"] = "; ".join(self.warnings)
        return d


def _atom(res, name: str):
    for a in res:
        if a.name == name and a.altloc in ("", "A", "\x00"):
            return a
    for a in res:
        if a.name == name:
            return a
    return None


def _seq_window(chain_residues: List, index: int, width: int = 2) -> str:
    lo = max(0, index - width)
    hi = min(len(chain_residues), index + width + 1)
    out = []
    for r in chain_residues[lo:hi]:
        info = gemmi.find_tabulated_residue(r.name)
        code = info.one_letter_code.upper() if info and info.one_letter_code else "X"
        out.append(code if code.strip() else "X")
    return "".join(out)


@dataclass
class _SiteAtom:
    chain: str
    subchain: str
    res_index: int
    res: object
    rlist: List
    atom: object


def collect_site_atoms(st: "gemmi.Structure", enzyme_chain_names: Set[str],
                       exclude: Optional[Set[Tuple[str, int]]] = None
                       ) -> Tuple[List[_SiteAtom], List[_SiteAtom], List[_SiteAtom]]:
    """Pre-filter the enzyme to His / elbow / acid atoms only.

    Doing this once per model keeps the triad search at a few thousand
    distance evaluations instead of a full combinatorial sweep, which matters
    when scanning tens of thousands of poses.
    """
    model = st[0]
    his: List[_SiteAtom] = []
    elbow: List[_SiteAtom] = []
    acid: List[_SiteAtom] = []
    for chain in model:
        if enzyme_chain_names and chain.name not in enzyme_chain_names:
            continue
        rlist = [r for r in chain]
        for idx, res in enumerate(rlist):
            if exclude and (chain.name, res.seqid.num) in exclude:
                continue
            rn = res.name.upper()
            if rn in HIS_ALIASES:
                ce1 = _atom(res, "CE1")
                if ce1 is not None:
                    his.append(_SiteAtom(chain.name, res.subchain, idx, res, rlist, ce1))
            if rn in ELBOW_ATOMS:
                for name in ELBOW_ATOMS[rn]:
                    a = _atom(res, name)
                    if a is not None:
                        elbow.append(_SiteAtom(chain.name, res.subchain, idx, res, rlist, a))
            if rn in ACID_RESIDUES:
                for name in ACID_RESIDUES[rn]:
                    a = _atom(res, name)
                    if a is not None:
                        acid.append(_SiteAtom(chain.name, res.subchain, idx, res, rlist, a))
    return his, elbow, acid


def _residue_key(elbow_res: int, acid_res: int, his_res: int) -> str:
    return f"{elbow_res}|{acid_res}|{his_res}"


def build_forced_triad(st: "gemmi.Structure", enzyme_chain_names: Set[str],
                       want: Tuple[int, int, int],
                       exclude: Optional[Set[Tuple[str, int]]] = None) -> Optional[Triad]:
    """Rebuild a specific triad by residue number, for ensemble consensus.

    Once the ensemble has voted on which three residues form the triad, every
    pose is measured against those same residues, so a single distorted model
    cannot silently relocate the active site.
    """
    e_want, a_want, h_want = want
    his_atoms, elbow_atoms, acid_atoms = collect_site_atoms(st, enzyme_chain_names, exclude)
    h = next((x for x in his_atoms if x.res.seqid.num == h_want), None)
    if h is None:
        return None
    nd1, ne2 = _atom(h.res, "ND1"), _atom(h.res, "NE2")
    if nd1 is None or ne2 is None:
        return None
    es = [x for x in elbow_atoms if x.res.seqid.num == e_want]
    as_ = [x for x in acid_atoms if x.res.seqid.num == a_want]
    if not es or not as_:
        return None
    best = None
    for imid_e, imid_a, orient in ((ne2, nd1, "canonical"), (nd1, ne2, "flipped")):
        e = min(es, key=lambda x: imid_e.pos.dist(x.atom.pos))
        a = min(as_, key=lambda x: imid_a.pos.dist(x.atom.pos))
        d_e = imid_e.pos.dist(e.atom.pos)
        d_a = imid_a.pos.dist(a.atom.pos)
        if best is None or (d_e + d_a) < best[0]:
            best = (d_e + d_a, d_e, d_a, e, a, orient)
    _, d_e, d_a, e, a, orient = best
    plddts = [x.b_iso for x in (h.atom, nd1, ne2, e.atom, a.atom)]
    motif = _seq_window(e.rlist, e.res_index)
    return Triad(
        found=True, his_chain=h.chain, his_label_chain=h.subchain,
        his_seqid=h.res.seqid.num, his_resname=h.res.name,
        elbow_chain=e.chain, elbow_seqid=e.res.seqid.num,
        elbow_resname=e.res.name, elbow_atom=e.atom.name,
        acid_chain=a.chain, acid_seqid=a.res.seqid.num,
        acid_resname=a.res.name, acid_atom=a.atom.name,
        d_his_elbow=round(d_e, 3), d_his_acid=round(d_a, 3),
        orientation=orient, motif=motif,
        motif_ok=bool(len(motif) == 5 and motif[0] == "G" and motif[4] == "G"),
        mutant=e.res.name.upper() in {"ALA", "GLY"}, method="ensemble_consensus",
        topology_ok=(e.res.seqid.num < a.res.seqid.num < h.res.seqid.num),
        residue_key=_residue_key(e.res.seqid.num, a.res.seqid.num, h.res.seqid.num),
        forced=True, plddt_mean=round(sum(plddts) / len(plddts), 2),
        ce1=h.atom.pos, ne2=ne2.pos, nd1=nd1.pos, elbow_pos=e.atom.pos)


def detect_triad(st: "gemmi.Structure", enzyme_chain_names: Set[str],
                 substrate_center: Optional["gemmi.Position"],
                 his_resid_hint: Optional[Tuple[int, int]] = None,
                 forced_his_seqid: Optional[int] = None,
                 exclude: Optional[Set[Tuple[str, int]]] = None,
                 require_topology: bool = True,
                 forced_triad: Optional[Tuple[int, int, int]] = None) -> Triad:
    """Geometric Ser/Cys/Ala-His-Asp/Glu triad search, ring-flip tolerant.

    Canonical alpha/beta-hydrolase geometry is His NE2 to the nucleophile
    elbow and His ND1 to the buried acid.  AF3 sometimes writes the imidazole
    flipped, so both assignments are tested and the winner is reported in the
    `orientation` column.  A Ser->Ala mutant is handled by treating Ala CB as
    the position formerly occupied by OG.
    """
    if forced_triad is not None:
        t = build_forced_triad(st, enzyme_chain_names, forced_triad, exclude)
        if t is not None:
            return t
    his_atoms, elbow_atoms, acid_atoms = collect_site_atoms(
        st, enzyme_chain_names, exclude)
    if not his_atoms:
        return Triad(found=False, method="no_histidine_in_enzyme_chain",
                     warnings=["no histidine with a CE1 atom found in the enzyme chain(s)"])
    if forced_his_seqid is not None:
        his_atoms = [h for h in his_atoms if h.res.seqid.num == forced_his_seqid] or his_atoms

    cutoffs = (4.0, 4.5, 5.5, 7.0)
    result = _triad_search(his_atoms, elbow_atoms, acid_atoms, cutoffs,
                           substrate_center, his_resid_hint, require_topology)
    if result is None and require_topology:
        result = _triad_search(his_atoms, elbow_atoms, acid_atoms, cutoffs,
                               substrate_center, his_resid_hint, False)
        if result is not None:
            result.warnings.append(
                "triad violates the canonical elbow < acid < His sequence order; "
                "accepted only because no canonical arrangement was found")
    if result is not None:
        return result
    return _triad_fallback(his_atoms, substrate_center, his_resid_hint)


def _triad_search(his_atoms, elbow_atoms, acid_atoms, cutoffs, substrate_center,
                  his_resid_hint, require_topology) -> Optional[Triad]:
    for cutoff in cutoffs:
        scored: List[Tuple[float, Triad]] = []
        for h in his_atoms:
            nd1, ne2 = _atom(h.res, "ND1"), _atom(h.res, "NE2")
            ce1 = h.atom
            if nd1 is None or ne2 is None:
                continue
            for (imid_e, imid_a, orient) in ((ne2, nd1, "canonical"),
                                             (nd1, ne2, "flipped")):
                near_elbow = [(imid_e.pos.dist(e.atom.pos), e) for e in elbow_atoms
                              if imid_e.pos.dist(e.atom.pos) <= cutoff
                              and not (e.chain == h.chain and e.res.seqid.num == h.res.seqid.num)]
                if not near_elbow:
                    continue
                near_acid = [(imid_a.pos.dist(a.atom.pos), a) for a in acid_atoms
                             if imid_a.pos.dist(a.atom.pos) <= cutoff]
                if not near_acid:
                    continue
                for d_e, e in near_elbow:
                    for d_a, a in near_acid:
                        same_chain = (e.chain == h.chain == a.chain)
                        topo_ok = ((e.res.seqid.num < a.res.seqid.num < h.res.seqid.num)
                                   if same_chain else True)
                        if require_topology and not topo_ok:
                            continue
                        motif = _seq_window(e.rlist, e.res_index)
                        motif_ok = bool(len(motif) == 5 and motif[0] == "G" and motif[4] == "G")
                        score = d_e + d_a
                        score += 0.35 * ELBOW_PRIORITY.get(e.res.name.upper(), 5)
                        if not motif_ok:
                            score += 0.8
                        if orient == "flipped":
                            score += 0.15
                        if substrate_center is not None:
                            score += 0.05 * ce1.pos.dist(substrate_center)
                        if his_resid_hint and not (
                                his_resid_hint[0] <= h.res.seqid.num <= his_resid_hint[1]):
                            score += 0.5
                        plddts = [x.b_iso for x in (ce1, imid_e, imid_a, e.atom, a.atom)]
                        scored.append((score, Triad(
                            found=True, his_chain=h.chain, his_label_chain=h.subchain,
                            his_seqid=h.res.seqid.num, his_resname=h.res.name,
                            elbow_chain=e.chain, elbow_seqid=e.res.seqid.num,
                            elbow_resname=e.res.name, elbow_atom=e.atom.name,
                            acid_chain=a.chain, acid_seqid=a.res.seqid.num,
                            acid_resname=a.res.name, acid_atom=a.atom.name,
                            d_his_elbow=round(d_e, 3), d_his_acid=round(d_a, 3),
                            orientation=orient, motif=motif, motif_ok=motif_ok,
                            mutant=e.res.name.upper() in {"ALA", "GLY"},
                            cutoff_used=cutoff, method="geometric",
                            topology_ok=topo_ok,
                            residue_key=_residue_key(e.res.seqid.num, a.res.seqid.num,
                                                     h.res.seqid.num),
                            plddt_mean=round(sum(plddts) / len(plddts), 2),
                            ce1=ce1.pos, ne2=ne2.pos, nd1=nd1.pos,
                            elbow_pos=e.atom.pos)))
        if scored:
            scored.sort(key=lambda x: x[0])
            best = scored[0][1]
            best.n_candidates = len(scored)
            if cutoff > cutoffs[0]:
                best.warnings.append(
                    f"triad only resolved after widening the H-bond cutoff to {cutoff} A; "
                    "active-site geometry may be distorted in this model")
            if best.mutant:
                best.warnings.append(
                    f"elbow is {best.elbow_resname}{best.elbow_seqid}: treated as a "
                    "Ser->Ala/Gly active-site mutant, CB/CA used as the Ser OG proxy")
            if not best.motif_ok:
                best.warnings.append(
                    f"elbow sequence window '{best.motif}' does not match G-x-S-x-G")
            if len({(t.his_chain, t.his_seqid) for _, t in scored}) > 1:
                best.warnings.append(
                    "more than one histidine satisfied the triad criteria; kept the "
                    "best-scoring one")
            return best

    return None


def _triad_fallback(his_atoms, substrate_center, his_resid_hint) -> Triad:
    # No triad geometry anywhere. Use the His nearest the substrate.
    pool = his_atoms
    if his_resid_hint:
        pool = [h for h in his_atoms
                if his_resid_hint[0] <= h.res.seqid.num <= his_resid_hint[1]] or his_atoms
    if substrate_center is not None:
        pool = sorted(pool, key=lambda h: h.atom.pos.dist(substrate_center))
        h = pool[0]
        nd1, ne2 = _atom(h.res, "ND1"), _atom(h.res, "NE2")
        return Triad(
            found=True, his_chain=h.chain, his_label_chain=h.subchain,
            his_seqid=h.res.seqid.num, his_resname=h.res.name,
            method="fallback_nearest_his_to_substrate", ce1=h.atom.pos,
            ne2=ne2.pos if ne2 else None, nd1=nd1.pos if nd1 else None,
            plddt_mean=round(h.atom.b_iso, 2),
            warnings=["NO Ser/Ala-His-Asp geometry found; fell back to the histidine "
                      "closest to the substrate. Distances from this model are "
                      "unvalidated and should be inspected manually."])
    return Triad(found=False, method="failed",
                 warnings=["no triad and no fallback histidine could be assigned"])


# ---------------------------------------------------------------------------
# Structure loading and entity classification
# ---------------------------------------------------------------------------


@dataclass
class SubstrateUnit:
    kind: str                      # ligand | short_polymer | mixed
    chain: str
    label_chain: str
    resnames: str
    seqid_first: int
    n_atoms: int
    seqids: List[int] = field(default_factory=list)
    atom_idxs: List[int] = field(default_factory=list)


def _pick_block(doc):
    for block in doc:
        if block.find_loop("_atom_site.id") or block.find_values("_atom_site.id"):
            return block
    return doc.sole_block()


def parse_chem_comp_bonds(block) -> Dict[str, List[Tuple[str, str, str, bool]]]:
    """Read authoritative ligand connectivity if the CIF carries it."""
    out: Dict[str, List[Tuple[str, str, str, bool]]] = defaultdict(list)
    try:
        table = block.find("_chem_comp_bond.",
                           ["comp_id", "atom_id_1", "atom_id_2",
                            "?value_order", "?pdbx_aromatic_flag"])
    except Exception:
        return {}
    for row in table:
        try:
            comp = row.str(0).upper()
            a1, a2 = row.str(1), row.str(2)
            order = row.str(3).upper() if row.has(3) else "SING"
            arom = (row.str(4).upper().startswith("Y")) if row.has(4) else False
        except Exception:
            continue
        out[comp].append((a1, a2, order or "SING", arom))
    return dict(out)


def read_model(path: str) -> Tuple["gemmi.Structure", Dict[str, List[Tuple[str, str, str, bool]]]]:
    lower = path.lower()
    if lower.endswith((".cif", ".mmcif", ".cif.gz", ".mmcif.gz")):
        doc = gemmi.cif.read(path)
        block = _pick_block(doc)
        st = gemmi.make_structure_from_block(block)
        comp_bonds = parse_chem_comp_bonds(block)
    else:
        st = gemmi.read_structure(path)
        comp_bonds = {}
    st.setup_entities()
    # AF3 models are not crystals; make sure no periodic images are considered.
    st.spacegroup_hm = "P 1"
    return st, comp_bonds


def clone_atom(a: GraphAtom) -> GraphAtom:
    return GraphAtom(key=a.key, chain=a.chain, label_chain=a.label_chain,
                     seqid=a.seqid, resname=a.resname, name=a.name,
                     element=a.element, pos=a.pos, bfactor=a.bfactor,
                     is_polymer=a.is_polymer)


def build_atom_list(st: "gemmi.Structure", enzyme_chain_names: Set[str],
                    substrate_keys: Optional[Set[Tuple[str, int]]] = None
                    ) -> List[GraphAtom]:
    """Flatten the model to graph atoms.

    is_polymer is decided per residue, not per chain, so a ligand that shares
    the enzyme chain id is still treated as a separate entity.
    """
    substrate_keys = substrate_keys or set()
    atoms: List[GraphAtom] = []
    for chain in st[0]:
        chain_is_enzyme = chain.name in enzyme_chain_names
        for res in chain:
            is_poly = (chain_is_enzyme
                       and (chain.name, res.seqid.num) not in substrate_keys)
            for at in res:
                if at.element == gemmi.Element("H"):
                    continue
                atoms.append(GraphAtom(
                    key=(chain.name, res.seqid.num, res.name, at.name),
                    chain=chain.name, label_chain=res.subchain,
                    seqid=res.seqid.num, resname=res.name, name=at.name,
                    element=elem(at), pos=at.pos, bfactor=at.b_iso,
                    is_polymer=is_poly))
    for i, a in enumerate(atoms):
        a.idx = i
    return atoms


def _is_solvent_or_ion(res) -> bool:
    if res.name.upper() in NON_SUBSTRATE_RESNAMES:
        return True
    heavy = [a for a in res if a.element != gemmi.Element("H")]
    if not heavy:
        return True
    if all(elem(a) in METAL_OR_ION_ELEMENTS for a in heavy):
        return True
    return not any(elem(a) == "C" for a in heavy)


def _unit_from_residues(chain_name: str, residues: List, kind_hint: str
                        ) -> Optional[SubstrateUnit]:
    if not residues:
        return None
    heavy = sum(len([a for a in r if a.element != gemmi.Element("H")]) for r in residues)
    names = [r.name for r in residues]
    if len(names) == 1:
        kind = "ligand"
    elif all(n.upper() not in STANDARD_AA for n in names):
        kind = kind_hint or "ligand"
    elif all(n.upper() in STANDARD_AA for n in names):
        kind = "short_polymer"
    else:
        kind = "mixed"
    label = ",".join(names)
    return SubstrateUnit(
        kind=kind, chain=chain_name, label_chain=residues[0].subchain,
        resnames=(label[:117] + "...") if len(label) > 120 else label,
        seqid_first=residues[0].seqid.num, n_atoms=heavy,
        seqids=[r.seqid.num for r in residues])


def classify_entities(st: "gemmi.Structure", min_enzyme_len: int = 40,
                      min_ligand_atoms: int = 8
                      ) -> Tuple[Set[str], List[SubstrateUnit], List[str]]:
    """Decide which chains are enzyme and which residues are the substrate.

    A whole non-enzyme chain is taken as one substrate unit, so a peptide built
    from several non-standard CCD components does not get split into fragments.
    Non-solvent HETATM groups sitting inside an enzyme chain are also collected,
    which covers AF3 runs where the ligand shares the protein chain id.
    """
    warnings: List[str] = []
    enzyme: Set[str] = set()
    units: List[SubstrateUnit] = []
    model = st[0]
    for chain in model:
        poly = chain.get_polymer()
        plen = len(poly) if poly else 0
        keep = [r for r in chain if not _is_solvent_or_ion(r)]
        if plen >= min_enzyme_len:
            enzyme.add(chain.name)
            het = [r for r in keep if r.name.upper() not in STANDARD_AA]
            u = _unit_from_residues(chain.name, het, "ligand")
            if u and u.n_atoms >= min_ligand_atoms:
                u.kind = "ligand_in_enzyme_chain"
                units.append(u)
                warnings.append(
                    f"substrate candidate {u.resnames} shares chain {chain.name} "
                    "with the enzyme")
        else:
            u = _unit_from_residues(chain.name, keep, "ligand")
            if u and u.n_atoms >= min_ligand_atoms:
                units.append(u)
    if not enzyme:
        best, blen = None, -1
        for chain in model:
            poly = chain.get_polymer()
            L = len(poly) if poly else 0
            if L > blen:
                best, blen = chain.name, L
        if best is not None and blen >= 2:
            enzyme.add(best)
            warnings.append(
                f"no chain reached {min_enzyme_len} residues; treating chain {best} "
                f"({blen} residues) as the enzyme")
            units = [u for u in units if u.chain != best]
    if not units:
        warnings.append("no substrate ligand or short peptide chain identified")
    return enzyme, units, warnings


def neighbours_within(sub_pos: Sequence["gemmi.Position"],
                      other_pos: Sequence["gemmi.Position"],
                      cutoff: float) -> List[int]:
    """Indices of `other_pos` within `cutoff` of any point in `sub_pos`."""
    if not sub_pos or not other_pos:
        return []
    if np is not None:
        a = np.array([[p.x, p.y, p.z] for p in sub_pos])
        b = np.array([[p.x, p.y, p.z] for p in other_pos])
        lo = a.min(axis=0) - cutoff
        hi = a.max(axis=0) + cutoff
        mask = np.all((b >= lo) & (b <= hi), axis=1)
        cand = np.flatnonzero(mask)
        if cand.size == 0:
            return []
        d2 = ((b[cand][:, None, :] - a[None, :, :]) ** 2).sum(axis=2)
        keep = cand[(d2 <= cutoff * cutoff).any(axis=1)]
        return keep.tolist()
    cell = max(cutoff, 1e-3)
    grid: Dict[Tuple[int, int, int], List[int]] = defaultdict(list)
    for i, p in enumerate(other_pos):
        grid[(int(p.x // cell), int(p.y // cell), int(p.z // cell))].append(i)
    hits: Set[int] = set()
    for p in sub_pos:
        cx, cy, cz = int(p.x // cell), int(p.y // cell), int(p.z // cell)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for i in grid.get((cx + dx, cy + dy, cz + dz), ()):
                        if p.dist(other_pos[i]) <= cutoff:
                            hits.add(i)
    return sorted(hits)


# ---------------------------------------------------------------------------
# AF3 output discovery and confidence extraction
# ---------------------------------------------------------------------------

SEED_SAMPLE_RE = re.compile(r"seed[-_]?(\d+)[-_]sample[-_]?(\d+)", re.I)
SEED_ONLY_RE = re.compile(r"seed[-_]?(\d+)", re.I)
SAMPLE_ONLY_RE = re.compile(r"sample[-_]?(\d+)", re.I)
MODEL_N_RE = re.compile(r"model[-_](\d+)", re.I)
# Folders AF3 uses for bookkeeping rather than for identity. A model inside one of
# these is usually a second copy of a pose that also lives in its own seed folder.
BOOKKEEPING_RE = re.compile(
    r"^(collected|ranked|top|best|all)[-_]?models?$|^predictions?$", re.I)


def looks_like_coordinate_cif(path: str) -> bool:
    try:
        with open(path, "r", errors="ignore") as fh:
            head = fh.read(200000)
    except OSError:
        return False
    return "_atom_site." in head


def discover_models(roots: Sequence[str], patterns: Sequence[str]) -> List[str]:
    found: List[str] = []
    for root in roots:
        if os.path.isfile(root):
            found.append(root)
            continue
        for pat in patterns:
            found.extend(glob.glob(os.path.join(root, "**", pat), recursive=True))
    uniq = sorted({os.path.abspath(p) for p in found})
    return [p for p in uniq if looks_like_coordinate_cif(p) or p.lower().endswith(".pdb")]


def _num(text: str, rx) -> Optional[int]:
    m = rx.search(text)
    return int(m.group(1)) if m else None


def parse_identity(path: str) -> Dict[str, object]:
    """Recover job, seed, sample and the substrate folder from an AF3 output path.

    Layouts vary between AF3 versions and between run scripts: the seed can live
    in a directory called seed-1_sample-0, in one called <job>_seed0, or only in
    the filename as ..._seed000_model.cif. All three are read here, because a
    missing seed silently disables deduplication.

    The substrate folder matters too. When a series is organised as
    <substrate>/<job>/... the substrate label exists only in the path, and the
    residue name inside the CIF is whatever CCD code was assigned, which need not
    match. Both are reported so they can be compared.
    """
    parent = os.path.dirname(path)
    base = os.path.basename(path)
    segs = [x for x in parent.split(os.sep) if x]

    seed = sample = None
    for text in [base] + list(reversed(segs)):
        ms = SEED_SAMPLE_RE.search(text)
        if ms:
            seed, sample = int(ms.group(1)), int(ms.group(2))
            break
        only = _num(text, SEED_ONLY_RE)
        if only is not None:
            seed, sample = only, _num(text, SAMPLE_ONLY_RE)
            break
    if sample is None:
        sample = _num(base, SAMPLE_ONLY_RE)
        if sample is None:
            sample = _num(base, MODEL_N_RE)

    # The job directory is the model's own folder, unless that folder only
    # identifies a seed, a sample or a bookkeeping bucket, in which case climb.
    job_dir = parent
    while True:
        b = os.path.basename(job_dir.rstrip(os.sep))
        up = os.path.dirname(job_dir)
        if b and up and up != job_dir and (
                BOOKKEEPING_RE.match(b) or SEED_ONLY_RE.search(b)
                or SAMPLE_ONLY_RE.search(b)):
            job_dir = up
            continue
        break
    job = os.path.basename(job_dir.rstrip(os.sep)) or "job"

    collected = any(BOOKKEEPING_RE.match(x) for x in segs)
    # nearest ancestor that identifies the substrate rather than the run
    substrate_dir = ""
    for x in reversed(segs):
        if x == job or BOOKKEEPING_RE.match(x) or SEED_ONLY_RE.search(x):
            continue
        if x.startswith(job) or job.startswith(x):
            continue
        substrate_dir = x
        break
    return {"job": job, "job_dir": job_dir, "seed": seed, "sample": sample,
            "model_dir": parent, "model_file": base,
            "substrate_dir": substrate_dir, "collected_copy": collected}


def load_confidences(ident: Dict[str, object]) -> Dict[str, object]:
    """Pull AF3 confidence metrics from whatever sidecars exist."""
    out: Dict[str, object] = {}
    cands: List[str] = []
    for d in (ident["model_dir"], ident["job_dir"]):
        cands.extend(sorted(glob.glob(os.path.join(str(d), "*summary_confidences.json"))))
    for p in cands:
        try:
            with open(p) as fh:
                js = json.load(fh)
        except Exception:
            continue
        for key in ("ranking_score", "iptm", "ptm", "fraction_disordered", "has_clash"):
            if key in js and key not in out:
                out[key] = js[key]
        cp = js.get("chain_pair_pae_min")
        if cp and "pae_min_interchain" not in out:
            vals = [cp[i][j] for i in range(len(cp)) for j in range(len(cp[i])) if i != j]
            if vals:
                out["pae_min_interchain"] = round(min(vals), 3)
        ci = js.get("chain_iptm")
        if ci and "chain_iptm_min" not in out:
            try:
                out["chain_iptm_min"] = round(min(ci), 4)
            except Exception:
                pass
        if out.get("ranking_score") is not None:
            break
    if "ranking_score" not in out:
        rs = os.path.join(str(ident["job_dir"]), "ranking_scores.csv")
        if os.path.exists(rs) and ident.get("seed") is not None:
            try:
                with open(rs) as fh:
                    for row in csv.DictReader(fh):
                        if (int(row.get("seed", -1)) == ident["seed"]
                                and int(row.get("sample", -1)) == ident["sample"]):
                            out["ranking_score"] = float(row["ranking_score"])
                            break
            except Exception:
                pass
    return out


# ---------------------------------------------------------------------------
# Per-model scan
# ---------------------------------------------------------------------------

RANK_ATOM_FIELDS = {"CE1": "d_ce1", "NE2": "d_ne2", "ND1": "d_nd1",
                    "IMID_MIN": "d_imid_min", "ELBOW": "d_elbow",
                    "ATTACK": "d_attack"}


@dataclass
class ModelResult:
    path: str
    ok: bool
    ident: Dict[str, object]
    triad: Triad
    candidates: List[NucCandidate]
    substrate: Optional[SubstrateUnit]
    enzyme_chains: str
    electrophile: str
    electrophile_atom: str
    bond_source: str
    confidences: Dict[str, object]
    warnings: List[str]
    error: str = ""
    group_key: str = ""
    series_key: str = ""
    n_enzyme_residues: int = 0
    bond_pairs: List[Tuple[str, str]] = field(default_factory=list)
    ring_votes: List[Tuple[str, bool]] = field(default_factory=list)
    bonds_dropped: int = 0
    backbone_units: int = 0
    backbone_cyclic: bool = False
    backbone_note: str = ""
    backbone_fragments: int = 0
    backbone_coverage: float = 0.0
    backbone_ester_links: int = 0
    backbone_branches: int = 0


def group_key_for(ident: Dict[str, object], unit: "SubstrateUnit") -> str:
    """Key under which ensemble consensus is pooled: same job, same molecule."""
    return f"{ident.get('job')}|{unit.resnames}|{unit.n_atoms}"


def scan_model(path: str, opts: Dict[str, object],
               consensus: Optional[Dict[str, object]] = None) -> ModelResult:
    ident = parse_identity(path)
    warnings: List[str] = []
    try:
        st, comp_bonds = read_model(path)
        enzyme, subs, w = classify_entities(
            st, int(opts["min_enzyme_len"]), int(opts["min_ligand_atoms"]))
        warnings.extend(w)
        if opts.get("enzyme_chain"):
            enzyme = {str(opts["enzyme_chain"])}
        if opts.get("substrate_chain"):
            subs = [s for s in subs if s.chain == str(opts["substrate_chain"])] or subs
        if opts.get("substrate_resname"):
            want = str(opts["substrate_resname"]).upper()
            subs = [s for s in subs if want in s.resnames.upper()] or subs
        if not subs:
            return ModelResult(path, False, ident, Triad(found=False), [], None,
                               ",".join(sorted(enzyme)), "not_found", "", "none", {},
                               warnings + ["no substrate identified"], "no_substrate")
        # if several substrate units, keep the one with the most atoms
        subs.sort(key=lambda s: -s.n_atoms)
        if len(subs) > 1:
            warnings.append(
                "multiple substrate candidates ("
                + "; ".join(f"{s.chain}/{s.resnames}:{s.n_atoms}at" for s in subs[:4])
                + "); used the largest")
        unit = subs[0]
        sub_keys = {(unit.chain, sq) for sq in (unit.seqids or [unit.seqid_first])}
        atoms = build_atom_list(st, enzyme, sub_keys)
        unit.atom_idxs = [a.idx for a in atoms if (a.chain, a.seqid) in sub_keys]
        if not unit.atom_idxs:
            return ModelResult(path, False, ident, Triad(found=False), [], unit,
                               ",".join(sorted(enzyme)), "not_found", "", "none", {},
                               warnings + ["substrate atom lookup returned nothing"],
                               "substrate_atom_lookup_failed")
        sub_pos = [atoms[i].pos for i in unit.atom_idxs]
        cx = sum(p.x for p in sub_pos) / len(sub_pos)
        cy = sum(p.y for p in sub_pos) / len(sub_pos)
        cz = sum(p.z for p in sub_pos) / len(sub_pos)
        centroid = gemmi.Position(cx, cy, cz)

        gkey = group_key_for(ident, unit)
        n_enz_res = sum(1 for ch in st[0] if ch.name in enzyme for _ in ch)
        cons_triad = None
        cons_frac = None
        if consensus:
            ct = (consensus.get("triad") or {}).get(str(ident.get("job")))
            if ct:
                cons_triad = tuple(ct["residues"])
                cons_frac = ct.get("frac")
        hint = opts.get("his_window") or None
        triad = detect_triad(st, enzyme, centroid,
                             tuple(hint) if hint else None,
                             opts.get("force_his"), exclude=sub_keys,
                             require_topology=not opts.get("no_topology_filter"),
                             forced_triad=cons_triad)
        if cons_triad and not triad.forced:
            warnings.append(
                f"ensemble consensus triad {cons_triad} could not be built in this "
                "model; fell back to per-pose detection")
        triad.consensus_frac = cons_frac
        warnings.extend(triad.warnings)

        # graph over substrate + a thin enzyme shell (for tether detection)
        enz_idx = [a.idx for a in atoms if a.is_polymer]
        shell_local = neighbours_within(sub_pos, [atoms[i].pos for i in enz_idx], 2.6)
        shell = [enz_idx[i] for i in shell_local]
        graph_atoms = [atoms[i] for i in list(unit.atom_idxs) + shell]
        remap = {a.idx: k for k, a in enumerate(graph_atoms)}
        g = MolGraph([clone_atom(a) for a in graph_atoms])
        dropped = g.build_by_distance()
        if comp_bonds:
            g.apply_chem_comp_bonds(comp_bonds)
        sub_local = [remap[i] for i in unit.atom_idxs]
        sub_scope = set(sub_local)
        perceived_pairs = g.bond_name_pairs(sub_scope)
        cons_pairs = (consensus.get("bonds") or {}).get(gkey) if consensus else None
        if cons_pairs:
            g.apply_consensus_bonds([tuple(x) for x in cons_pairs], scope=sub_scope)
            changed = set(map(tuple, cons_pairs)) ^ set(perceived_pairs)
            if changed:
                warnings.append(
                    f"consensus connectivity differed from this pose's perception on "
                    f"{len(changed)} bond(s): " +
                    ", ".join(f"{a}-{b}" for a, b in sorted(changed)[:6]))
        if dropped:
            warnings.append(
                "valence cap dropped " + ", ".join(f"{a}-{b}@{d}A" for a, b, d in dropped[:6]))
        g.find_rings()
        if consensus and consensus.get("aromatic"):
            n_arom = g.apply_consensus_aromaticity(dict(consensus["aromatic"]))
            if n_arom:
                warnings.append(
                    f"{n_arom} ring(s) had their aromaticity set by the run-wide "
                    "consensus rather than by this pose's geometry")
        chain = trace_peptide(g, sub_local)

        if chain.n_fragments > 1 or chain.n_branches:
            warnings.append(chain.note or "incomplete backbone trace")
        cands = enumerate_nucleophiles(g, sub_local, chain)
        cons_nuc = (consensus.get("nucleophiles") or {}).get(gkey) if consensus else None
        if cons_nuc:
            missing = set(cons_nuc) - {c.atom_name for c in cands}
            flipped = 0
            for c in cands:
                v = cons_nuc.get(c.atom_name)
                if not v:
                    continue
                c.pose_class = c.nuc_class
                if c.nuc_class != v["class"] or c.accepted != v["accepted"]:
                    flipped += 1
                    c.nuc_class = str(v["class"])
                    c.accepted = bool(v["accepted"])
                    c.priority = CLASS_PRIORITY.get(c.nuc_class, 99)
                    c.class_group = CLASS_GROUP_OF.get(c.nuc_class, "")
                    c.reason = (f"{c.reason} [ensemble consensus: {v['class']}, "
                                f"this pose read it as {c.pose_class}]")
                c.class_source = "consensus"
                c.class_agreement = v["agreement"]
                c.presence = v["presence"]
                if not c.nuc_residue and v.get("residue") != "":
                    c.nuc_residue = v["residue"]
                if not c.nuc_residue_c and v.get("residue_c") != "":
                    c.nuc_residue_c = v.get("residue_c", "")
                if c.nuc_site in ("", "unassigned") and v.get("site"):
                    c.nuc_site = str(v["site"])
            if flipped:
                warnings.append(
                    f"{flipped} heteroatom classification(s) in this pose disagreed with "
                    "the ensemble; the consensus classification was used")
            if missing:
                warnings.append(
                    "ligand atoms in the ensemble inventory absent from this pose: "
                    + ", ".join(sorted(missing)[:8]))
        ref = triad.elbow_pos or triad.ce1 or centroid
        forced_e = None
        if consensus:
            ce = (consensus.get("electrophile") or {}).get(gkey)
            if ce:
                forced_e = str(ce["atom"])
        e_c, e_o, e_kind = find_electrophile(g, sub_local, ref, chain, forced_e)
        e_name = g.atoms[e_c].name if e_c is not None else ""

        for c in cands:
            p = g.atoms[c.atom_idx].pos
            if triad.ce1 is not None:
                c.d_ce1 = round(p.dist(triad.ce1), 3)
            if triad.ne2 is not None:
                c.d_ne2 = round(p.dist(triad.ne2), 3)
            if triad.nd1 is not None:
                c.d_nd1 = round(p.dist(triad.nd1), 3)
            ds = [d for d in (c.d_ne2, c.d_nd1) if d is not None]
            c.d_imid_min = min(ds) if ds else None
            if triad.elbow_pos is not None:
                c.d_elbow = round(p.dist(triad.elbow_pos), 3)
            if e_c is not None and e_c != c.atom_idx and e_c not in g.adj[c.atom_idx]:
                c.d_attack = round(p.dist(g.atoms[e_c].pos), 3)
                if e_o is not None:
                    c.bd_angle = round(math.degrees(gemmi.calculate_angle(
                        p, g.atoms[e_c].pos, g.atoms[e_o].pos)), 1)
            field_name = RANK_ATOM_FIELDS[str(opts["rank_atom"]).upper()]
            c.band = band_for(getattr(c, field_name))

        return ModelResult(path, True, ident, triad, cands, unit,
                           ",".join(sorted(enzyme)), e_kind, e_name,
                           g.bond_source, load_confidences(ident), warnings,
                           group_key=gkey,
                           series_key=(triad.residue_key or "?") + f"@{n_enz_res}",
                           n_enzyme_residues=n_enz_res,
                           bond_pairs=(perceived_pairs if opts.get("collect_bonds") else []),
                           ring_votes=(g.aromaticity_votes()
                                       if opts.get("collect_bonds") else []),
                           bonds_dropped=len(dropped),
                           backbone_units=chain.n_units,
                           backbone_cyclic=chain.cyclic,
                           backbone_note=chain.note,
                           backbone_fragments=chain.n_fragments,
                           backbone_coverage=chain.coverage,
                           backbone_ester_links=chain.n_ester_links,
                           backbone_branches=chain.n_branches)
    except Exception as exc:  # keep the sweep alive on a single bad file
        return ModelResult(path, False, ident, Triad(found=False), [], None, "",
                           "error", "", "none", {},
                           warnings + [f"exception: {exc}"],
                           traceback.format_exc(limit=3))


# ---------------------------------------------------------------------------
# Row assembly
# ---------------------------------------------------------------------------

def _conf(res: ModelResult, key: str):
    v = res.confidences.get(key)
    return v if v is not None else ""


def base_row(res: ModelResult, opts: Dict[str, object]) -> Dict[str, object]:
    t = res.triad
    s = res.substrate
    return {
        "job": res.ident.get("job"),
        "substrate_dir": res.ident.get("substrate_dir"),
        "collected_copy": res.ident.get("collected_copy"),
        "seed": res.ident.get("seed"),
        "sample": res.ident.get("sample"),
        "model_file": res.ident.get("model_file"),
        "enzyme_chains": res.enzyme_chains,
        "substrate_kind": s.kind if s else "",
        "substrate_chain": s.chain if s else "",
        "substrate_label_chain": s.label_chain if s else "",
        "substrate_resname": s.resnames if s else "",
        "substrate_seqid": s.seqid_first if s else "",
        "n_substrate_atoms": s.n_atoms if s else "",
        "bond_source": res.bond_source,
        "bonds_dropped_by_valence": res.bonds_dropped,
        "backbone_units": res.backbone_units,
        "backbone_cyclic": res.backbone_cyclic,
        "backbone_fragments": res.backbone_fragments,
        "backbone_coverage": res.backbone_coverage,
        "backbone_ester_links": res.backbone_ester_links,
        "backbone_branches": res.backbone_branches,
        "his_chain": t.his_chain, "his_label_chain": t.his_label_chain,
        "his_resname": t.his_resname, "his_seqid": t.his_seqid,
        "elbow": f"{t.elbow_resname}{t.elbow_seqid}.{t.elbow_atom}" if t.elbow_atom else "",
        "acid": f"{t.acid_resname}{t.acid_seqid}.{t.acid_atom}" if t.acid_atom else "",
        "d_his_elbow": t.d_his_elbow, "d_his_acid": t.d_his_acid,
        "triad_orientation": t.orientation, "triad_motif": t.motif,
        "triad_motif_ok": t.motif_ok, "triad_mutant": t.mutant,
        "triad_method": t.method, "triad_cutoff": t.cutoff_used,
        "triad_residue_key": t.residue_key, "triad_topology_ok": t.topology_ok,
        "triad_forced": t.forced, "triad_consensus_frac": t.consensus_frac,
        "triad_plddt_mean": t.plddt_mean,
        "electrophile_kind": res.electrophile,
        "electrophile_atom": res.electrophile_atom,
        "ranking_score": _conf(res, "ranking_score"),
        "iptm": _conf(res, "iptm"), "ptm": _conf(res, "ptm"),
        "pae_min_interchain": _conf(res, "pae_min_interchain"),
        "chain_iptm_min": _conf(res, "chain_iptm_min"),
        "has_clash": _conf(res, "has_clash"),
        "fraction_disordered": _conf(res, "fraction_disordered"),
        "rank_atom": str(opts["rank_atom"]).upper(),
        "model_path": res.path,
    }


def allowed_classes(spec: str) -> Set[str]:
    if not spec or spec.lower() in ("default", "auto"):
        return set(DEFAULT_CLASSES)
    if spec.lower() == "all":
        return set(CLASS_PRIORITY)
    out: Set[str] = set()
    for tok in spec.split(","):
        tok = tok.strip().lower()
        if not tok:
            continue
        if tok in CLASS_GROUPS:
            out |= CLASS_GROUPS[tok]
        else:
            out.add(tok)
    return out


def candidate_rows(res: ModelResult, opts: Dict[str, object]
                   ) -> Tuple[List[Dict[str, object]], Optional[Dict[str, object]]]:
    field_name = RANK_ATOM_FIELDS[str(opts["rank_atom"]).upper()]
    keep = allowed_classes(str(opts["classes"]))
    rows: List[Dict[str, object]] = []
    base = base_row(res, opts)
    usable = [c for c in res.candidates
              if c.accepted and c.nuc_class in keep and getattr(c, field_name) is not None]
    if str(opts["rank_by"]) == "priority":
        usable.sort(key=lambda c: (c.priority, getattr(c, field_name)))
    else:
        usable.sort(key=lambda c: (getattr(c, field_name), c.priority))
    best_key = (usable[0].atom_name, usable[0].seqid) if usable else None
    for c in res.candidates:
        if not (c.accepted and c.nuc_class in keep) and not opts["include_rejected"]:
            continue
        row = dict(base)
        row.update({
            "nuc_atom": c.atom_name, "nuc_element": c.element,
            "nuc_resname": c.resname, "nuc_seqid": c.seqid,
            "nuc_class": c.nuc_class, "nuc_accepted": c.accepted,
            "nuc_reason": c.reason, "nuc_priority": c.priority,
            "nuc_neighbours": c.neighbours, "nuc_plddt": round(c.plddt, 2),
            "nuc_residue": c.nuc_residue, "nuc_site": c.nuc_site,
            "nuc_residue_from_c": c.nuc_residue_c, "nuc_class_group": c.class_group,
            "nuc_class_source": c.class_source,
            "nuc_class_agreement": c.class_agreement,
            "nuc_pose_class": c.pose_class or c.nuc_class,
            "nuc_presence": c.presence,
            "d_ce1": c.d_ce1, "d_ne2": c.d_ne2, "d_nd1": c.d_nd1,
            "d_imid_min": c.d_imid_min, "d_elbow": c.d_elbow,
            "d_attack": c.d_attack, "burgi_dunitz_angle": c.bd_angle,
            "band": c.band,
            "is_best": bool(best_key and (c.atom_name, c.seqid) == best_key and c.accepted),
            "warnings": "; ".join(res.warnings),
        })
        rows.append(row)

    summary: Optional[Dict[str, object]] = None
    if usable:
        b = usable[0]
        summary = dict(base)
        summary.update({
            "best_nuc_atom": b.atom_name, "best_nuc_element": b.element,
            "best_nuc_class": b.nuc_class, "best_nuc_resname": b.resname,
            "best_nuc_seqid": b.seqid, "best_nuc_plddt": round(b.plddt, 2),
            "best_nuc_residue": b.nuc_residue, "best_nuc_site": b.nuc_site,
            "best_d_ce1": b.d_ce1, "best_d_ne2": b.d_ne2, "best_d_nd1": b.d_nd1,
            "best_d_imid_min": b.d_imid_min, "best_d_elbow": b.d_elbow,
            "best_d_attack": b.d_attack, "best_bd_angle": b.bd_angle,
            "best_band": b.band,
            "n_candidates_accepted": len(usable),
            "n_candidates_total": len(res.candidates),
            "n_near_attack": sum(1 for c in usable if c.band == "near_attack"),
            "runner_up_atom": usable[1].atom_name if len(usable) > 1 else "",
            "runner_up_d": getattr(usable[1], field_name) if len(usable) > 1 else "",
            "warnings": "; ".join(res.warnings),
        })
        for group, members in CLASS_GROUPS.items():
            sub = [c for c in usable if c.nuc_class in members]
            summary[f"best_{group}_atom"] = sub[0].atom_name if sub else ""
            summary[f"best_{group}_d"] = getattr(sub[0], field_name) if sub else ""
    else:
        summary = dict(base)
        summary.update({"best_nuc_atom": "", "best_band": "no_candidate",
                        "n_candidates_accepted": 0,
                        "n_candidates_total": len(res.candidates),
                        "warnings": "; ".join(res.warnings + ["no accepted nucleophile"])})
    return rows, summary


def inventory_rows(res: "ModelResult", opts: Dict[str, object]) -> List[Dict[str, object]]:
    """One row per ligand heteroatom per pose, regardless of --classes.

    The point is that nothing is invisible: an atom rejected as an amide, an
    ester oxygen or a tether still appears, with the reason. Without this, a
    nucleophile that a substitution introduced but the classifier declined is
    indistinguishable from one that was never there.
    """
    out = []
    for c in res.candidates:
        out.append({
            "job": res.ident.get("job"), "series_key": res.series_key,
            "substrate_dir": res.ident.get("substrate_dir", ""),
            "substrate": res.substrate.resnames if res.substrate else "",
            "n_substrate_atoms": res.substrate.n_atoms if res.substrate else "",
            "model_path": res.path,
            "atom": c.atom_name, "element": c.element, "nuc_class": c.nuc_class,
            "accepted": c.accepted, "residue": c.nuc_residue, "site": c.nuc_site,
            "residue_from_c": c.nuc_residue_c, "class_group": c.class_group,
            "class_source": c.class_source, "class_agreement": c.class_agreement,
            "presence": c.presence, "pose_class": c.pose_class or c.nuc_class,
            "reason": c.reason, "d_ce1": c.d_ce1, "d_ne2": c.d_ne2,
            "d_imid_min": c.d_imid_min, "d_attack": c.d_attack, "plddt": round(c.plddt, 2),
        })
    return out


def aggregate_inventory(rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    """Collapse per-pose inventory into one row per (job, substrate, atom)."""
    # Keyed on ligand size as well as name: a truncation series can reuse one CCD
    # name across variants of different length, and collapsing those would hide
    # exactly the added/removed nucleophiles this table exists to show.
    by: Dict[Tuple[str, str, str, str], List[Dict[str, object]]] = defaultdict(list)
    for r in rows:
        by[(str(r["job"]), str(r["substrate"]), str(r["n_substrate_atoms"]),
            str(r["atom"]))].append(r)
    out = []
    for (job, sub, size, atom), rs in sorted(by.items()):
        ds = [d for d in (_as_float(r.get("d_ce1")) for r in rs) if d is not None]
        cls = Counter(str(r["nuc_class"]) for r in rs).most_common(1)[0][0]
        pose_cls = Counter(str(r["pose_class"]) for r in rs)
        out.append({
            "job": job, "series_key": rs[0]["series_key"], "substrate": sub,
            "substrate_dir": rs[0].get("substrate_dir", ""),
            "n_substrate_atoms": rs[0]["n_substrate_atoms"],
            "variant": f"{sub}:{size}at",
            "atom": atom, "element": rs[0]["element"], "nuc_class": cls,
            "accepted": _as_bool(rs[0]["accepted"]),
            "residue": rs[0]["residue"], "site": rs[0]["site"],
            "residue_from_c": rs[0]["residue_from_c"],
            "class_group": rs[0].get("class_group", ""),
            "n_poses": len(rs),
            "class_agreement": _as_float(rs[0]["class_agreement"]),
            "presence": _as_float(rs[0]["presence"]),
            "pose_class_variants": len(pose_cls),
            "pose_class_breakdown": "; ".join(f"{k}:{v}" for k, v in pose_cls.most_common(3)),
            "min_d_ce1": round(min(ds), 3) if ds else "",
            "median_d_ce1": median(ds) if ds else "",
            "n_near_attack": sum(1 for d in ds if d < BAND_NEAR_ATTACK),
            "reason": rs[0]["reason"],
        })
    return out


def _as_float(v) -> Optional[float]:
    """CSV round-trips turn every value into a string; keep aggregation numeric."""
    if v in ("", None):
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _as_bool(v) -> bool:
    if isinstance(v, bool):
        return v
    return str(v).strip().lower() in ("true", "1", "yes")


def nucleophile_signature(rows: List[Dict[str, object]], anchor: str = "c",
                          detail: str = "group") -> List[str]:
    """Naming-independent description of a substrate's nucleophile set.

    Atom names change between CCD components, so a series cannot be compared on
    them. Position along the traced backbone plus site plus class survives
    renaming. The position is counted from the C-terminal acyl carbon by default,
    because that is the atom the TE attacks and the one a truncation series keeps
    fixed; numbering from the N-terminus shifts every index when the peptide is
    shortened, which manufactures spurious differences.
    """
    key = "residue_from_c" if anchor == "c" else "residue"
    ckey = "nuc_class" if detail == "fine" else "class_group"
    sig = []
    for r in rows:
        if not _as_bool(r["accepted"]):
            continue
        resid = r.get(key)
        resid = resid if resid not in ("", None) and str(resid) != "nan" else "?"
        cls = r.get(ckey) or CLASS_GROUP_OF.get(str(r["nuc_class"]), str(r["nuc_class"]))
        sig.append(f"res{resid}:{r['site']}:{cls}")
    return sorted(sig)


def class_signature(rows: List[Dict[str, object]], detail: str = "group") -> Counter:
    """Nucleophile counts by class: order-free, numbering-free, always comparable.

    `group` merges phenol with aliphatic hydroxyl and the amine subtypes, which
    removes the aromaticity judgement from the comparison. `fine` keeps them apart.
    """
    key = "nuc_class" if detail == "fine" else "class_group"
    out: Counter = Counter()
    for r in rows:
        if not _as_bool(r["accepted"]):
            continue
        v = r.get(key) or CLASS_GROUP_OF.get(str(r["nuc_class"]), str(r["nuc_class"]))
        out[str(v)] += 1
    return out


def _pick_anchor(entries, detail: str) -> str:
    """Choose the residue anchor that actually lines the series up.

    Which end is truncated is a property of the experiment, not something the
    tool can assume. Numbering from the C-terminal acyl carbon aligns a series
    shortened from the N-terminus; numbering from the N-terminus aligns one
    shortened from the C-terminal end. Both are computed and the anchor with more
    matching nucleophiles against the reference wins, so neither choice has to be
    made in advance or supplied by the user.
    """
    scores: Dict[str, int] = {}
    for anchor in ("n", "c"):
        ref = Counter(nucleophile_signature(entries[0][2], anchor, detail))
        scores[anchor] = sum(
            sum((Counter(nucleophile_signature(rs, anchor, detail)) & ref).values())
            for _, _, rs in entries)
    return "n" if scores["n"] > scores["c"] else "c"


def series_rows(inv: List[Dict[str, object]], anchor: str = "auto",
                detail: str = "group") -> List[Dict[str, object]]:
    """Compare substrates that were run against the same enzyme.

    Rows are ordered by ligand size so a truncation series reads top to bottom.
    Two diffs are reported. The class-count diff is unconditionally reliable: it
    cannot show a gain for a ligand whose nucleophile set is a subset. The
    positional diff is more informative when it works, but depends on residue
    numbering lining up between variants, so `positional_diff_agrees` flags the
    cases where the two disagree and the positional columns should be ignored.
    """
    by_sub: Dict[Tuple[str, str, str, str], List[Dict[str, object]]] = defaultdict(list)
    for r in inv:
        by_sub[(str(r["series_key"]), str(r["substrate"]),
                str(r["n_substrate_atoms"]), str(r["job"]))].append(r)
    by_series: Dict[str, List[Tuple[str, str, List[Dict[str, object]]]]] = defaultdict(list)
    for (sk, sub, size, job), rs in by_sub.items():
        by_series[sk].append((sub, job, rs))
    out = []
    for sk, entries in sorted(by_series.items()):
        entries.sort(key=lambda e: -(_as_float(e[2][0]["n_substrate_atoms"]) or 0))
        use = _pick_anchor(entries, detail) if anchor == "auto" else anchor
        ref_sub, _, ref_rows = entries[0]
        ref_pos = Counter(nucleophile_signature(ref_rows, use, detail))
        ref_cls = class_signature(ref_rows, detail)
        for sub, job, rs in entries:
            pos = Counter(nucleophile_signature(rs, use, detail))
            cls = class_signature(rs, detail)
            gained_pos, lost_pos = pos - ref_pos, ref_pos - pos
            gained_cls, lost_cls = cls - ref_cls, ref_cls - cls
            acc = [r for r in rs if _as_bool(r["accepted"])]
            best = min((d for d in (_as_float(r["min_d_ce1"]) for r in acc)
                        if d is not None), default="")
            agrees = (sum(gained_pos.values()) == sum(gained_cls.values())
                      and sum(lost_pos.values()) == sum(lost_cls.values()))
            out.append({
                "series_key": sk, "job": job, "substrate": sub,
                "n_substrate_atoms": rs[0]["n_substrate_atoms"],
                "variant": rs[0].get("variant", sub),
                "substrate_dir": rs[0].get("substrate_dir", ""),
                "reference_substrate": ref_sub,
                "n_heteroatoms": len(rs), "n_nucleophiles": len(acc),
                "n_rejected": len(rs) - len(acc),
                "classes": "; ".join(f"{k}:{v}" for k, v in sorted(cls.items())),
                "gained_classes": "; ".join(
                    f"{k}:{v}" for k, v in sorted(gained_cls.items())),
                "lost_classes": "; ".join(
                    f"{k}:{v}" for k, v in sorted(lost_cls.items())),
                "n_gained": sum(gained_cls.values()),
                "n_lost": sum(lost_cls.values()),
                "positional_diff_agrees": agrees,
                "gained_positional": "; ".join(sorted(gained_pos.elements())),
                "lost_positional": "; ".join(sorted(lost_pos.elements())),
                "anchor": f"residue_from_{use}",
                "anchor_selected": "auto" if anchor == "auto" else "fixed",
                "class_detail": detail,
                "chains": rs[0].get("chains", ""),
                "branches": rs[0].get("branches", ""),
                "unstable_atoms": "; ".join(
                    f"{r['atom']}({(_as_float(r['class_agreement']) or 1.0):.0%})"
                    for r in rs
                    if (_as_float(r["class_agreement"]) or 1.0) < 1.0),
                "best_d_ce1": best,
                "signature": " | ".join(sorted(pos.elements())),
            })
    return out


def median(values: Sequence[float]) -> Optional[float]:
    v = sorted(values)
    if not v:
        return None
    n = len(v)
    return v[n // 2] if n % 2 else round((v[n // 2 - 1] + v[n // 2]) / 2, 3)


def aggregate_jobs(model_summaries: List[Dict[str, object]],
                   opts: Dict[str, object]) -> List[Dict[str, object]]:
    by_job: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for r in model_summaries:
        by_job[str(r.get("job"))].append(r)
    metric = {"CE1": "best_d_ce1", "NE2": "best_d_ne2", "ND1": "best_d_nd1",
              "IMID_MIN": "best_d_imid_min", "ELBOW": "best_d_elbow",
              "ATTACK": "best_d_attack"}[str(opts["rank_atom"]).upper()]
    out: List[Dict[str, object]] = []
    for job, rows in sorted(by_job.items()):
        ds = [float(r[metric]) for r in rows
              if r.get(metric) not in ("", None)]
        atoms = [str(r.get("best_nuc_atom")) for r in rows if r.get("best_nuc_atom")]
        counts = Counter(atoms)
        modal, modal_n = (counts.most_common(1)[0] if counts else ("", 0))
        bands = Counter(str(r.get("best_band")) for r in rows)
        rank = [float(r["ranking_score"]) for r in rows
                if r.get("ranking_score") not in ("", None)]
        best_row = min((r for r in rows if r.get(metric) not in ("", None)),
                       key=lambda r: float(r[metric]), default=None)
        out.append({
            "job": job,
            "n_models": len(rows),
            "n_with_distance": len(ds),
            "min_d": min(ds) if ds else "",
            "median_d": median(ds) if ds else "",
            "mean_d": round(sum(ds) / len(ds), 3) if ds else "",
            "max_d": max(ds) if ds else "",
            "n_near_attack": bands.get("near_attack", 0),
            "frac_near_attack": round(bands.get("near_attack", 0) / len(rows), 4) if rows else "",
            "n_intermediate": bands.get("intermediate", 0),
            "n_too_far": bands.get("too_far", 0),
            "n_no_candidate": bands.get("no_candidate", 0),
            "modal_nucleophile": modal,
            "modal_nucleophile_class": next(
                (str(r.get("best_nuc_class") or "") for r in rows
                 if str(r.get("best_nuc_atom") or "") == modal), ""),
            "nucleophile_agreement": round(modal_n / len(atoms), 4) if atoms else "",
            "mean_ranking_score": round(sum(rank) / len(rank), 4) if rank else "",
            "modal_his": Counter(str(r.get("his_seqid")) for r in rows).most_common(1)[0][0],
            "modal_triad_method": Counter(
                str(r.get("triad_method")) for r in rows).most_common(1)[0][0],
            "modal_elbow": Counter(str(r.get("elbow")) for r in rows).most_common(1)[0][0],
            "best_pose_seed": best_row.get("seed") if best_row else "",
            "best_pose_sample": best_row.get("sample") if best_row else "",
            "best_pose_path": best_row.get("model_path") if best_row else "",
            "rank_atom": str(opts["rank_atom"]).upper(),
        })
    return out


def write_csv(path: str, rows: List[Dict[str, object]], preferred: Sequence[str] = ()) -> None:
    if not rows:
        with open(path, "w") as fh:
            fh.write("")
        return
    keys: List[str] = [k for k in preferred if any(k in r for r in rows)]
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


LONG_ORDER = ["job", "seed", "sample", "nuc_atom", "nuc_element", "nuc_class",
              "nuc_residue", "nuc_site",
              "d_ce1", "band", "is_best", "d_ne2", "d_nd1", "d_imid_min",
              "d_elbow", "d_attack", "burgi_dunitz_angle", "nuc_plddt",
              "his_resname", "his_seqid", "elbow", "acid", "triad_method",
              "triad_orientation", "ranking_score"]
INVENTORY_ORDER = ["job", "substrate_dir", "substrate", "variant", "n_substrate_atoms", "atom", "element",
                   "nuc_class", "class_group", "accepted", "residue", "residue_from_c",
                   "site", "presence",
                   "class_agreement", "pose_class_variants", "min_d_ce1",
                   "median_d_ce1", "n_near_attack", "n_poses", "reason"]
SERIES_ORDER = ["series_key", "job", "substrate_dir", "substrate", "variant",
                "n_substrate_atoms",
                "n_nucleophiles", "n_heteroatoms", "n_rejected", "classes",
                "gained_classes", "lost_classes", "n_gained", "n_lost",
                "positional_diff_agrees", "gained_positional", "lost_positional",
                "best_d_ce1", "unstable_atoms", "anchor", "anchor_selected",
                "class_detail",
                "reference_substrate",
                "signature"]
SUMMARY_ORDER = ["job", "seed", "sample", "best_nuc_atom", "best_nuc_class",
                 "best_nuc_residue", "best_nuc_site",
                 "best_d_ce1", "best_band", "best_d_attack", "best_bd_angle",
                 "best_nuc_plddt", "n_candidates_accepted", "runner_up_atom",
                 "runner_up_d", "his_resname", "his_seqid", "elbow", "acid",
                 "triad_method", "ranking_score", "iptm"]


def write_pml(res: ModelResult, best_atom: str, out_path: str) -> None:
    t = res.triad
    s = res.substrate
    if not (s and t.found and t.his_seqid > 0):
        return
    lines = [
        f"# generated by te_autoscan {__version__}",
        f"load {res.path}, model",
        "hide everything", "bg_color white",
        "show cartoon, polymer", "color grey80, polymer",
        f"select nuc, (chain {s.chain} and resi {s.seqid_first} and name {best_atom})"
        if len(s.seqids) <= 1 else
        f"select nuc, (chain {s.chain} and name {best_atom})",
        f"select cathis, (chain {t.his_chain} and resi {t.his_seqid})",
        f"select elbow, (chain {t.elbow_chain} and resi {t.elbow_seqid})" if t.elbow_seqid > 0 else "",
        f"select acid, (chain {t.acid_chain} and resi {t.acid_seqid})" if t.acid_seqid > 0 else "",
        "select substrate, (chain %s and resi %s)" % (
            s.chain, "+".join(str(x) for x in (s.seqids or [s.seqid_first]))),
        "show sticks, cathis or elbow or acid or substrate",
        "color cyan, substrate", "color orange, cathis",
        f"distance d_ce1, nuc, (cathis and name CE1)",
        f"distance d_ne2, nuc, (cathis and name NE2)",
        "set dash_color, red", "set label_size, 16",
        "orient substrate or cathis", "zoom substrate or cathis, 4",
    ]
    with open(out_path, "w") as fh:
        fh.write("\n".join([L for L in lines if L]) + "\n")


# ---------------------------------------------------------------------------
# Worker + driver
# ---------------------------------------------------------------------------

def _survey_worker(payload):
    """Cheap pass used only to gather votes: connectivity, rings, triad, nucleophiles."""
    path, opts, pre = payload
    o = dict(opts); o["collect_bonds"] = True
    try:
        res = scan_model(path, o, consensus=pre)
        return {"ok": res.ok, "job": str(res.ident.get("job")),
                "gkey": res.group_key, "pairs": res.bond_pairs,
                "elec_atom": res.electrophile_atom, "rings": res.ring_votes,
                "nucs": [(c.atom_name, c.element, c.nuc_class, c.accepted,
                          c.nuc_residue, c.nuc_site, c.nuc_residue_c)
                         for c in res.candidates],
                "triad_key": res.triad.residue_key if res.triad.found else "",
                "triad_method": res.triad.method,
                "topology_ok": res.triad.topology_ok}
    except Exception:
        return {"ok": False, "job": str(parse_identity(path).get("job")), "gkey": "",
                "pairs": [], "elec_atom": "", "nucs": [], "rings": [],
                "triad_key": "", "triad_method": "", "topology_ok": False}


def run_survey(models: Sequence[str], opts: Dict[str, object], n_jobs: int,
               per_job: int, bond_vote: float, quiet: bool = False) -> Dict[str, object]:
    """Vote on connectivity and on the catalytic triad across the ensemble.

    The same molecule and the same protein appear in every pose of a job, so the
    parts of the analysis that must not vary with conformation are decided once,
    by majority, and then applied to every pose. Nothing is supplied by the user.
    """
    by_job: Dict[str, List[str]] = defaultdict(list)
    for m in models:
        by_job[str(parse_identity(m).get("job"))].append(m)
    sample: List[str] = []
    for job, paths in by_job.items():
        paths = sorted(paths)
        if len(paths) <= per_job:
            sample.extend(paths)
        else:
            step = len(paths) / per_job
            sample.extend(paths[int(i * step)] for i in range(per_job))
    if not quiet:
        print(f"  survey: voting on connectivity and triad from {len(sample)} of "
              f"{len(models)} poses", flush=True)

    def _sweep(pre):
        payloads = [(m, opts, pre) for m in sample]
        if n_jobs > 1:
            from concurrent.futures import ProcessPoolExecutor
            with ProcessPoolExecutor(max_workers=n_jobs) as ex:
                return list(ex.map(_survey_worker, payloads, chunksize=4))
        return [_survey_worker(pl) for pl in payloads]

    outs = _sweep(None)

    # Aromaticity is a property of the molecule, so it is decided once for the
    # whole run rather than per pose or per job. Judging it from each pose's ring
    # geometry lets the same chemical group read as a phenol in one substrate and
    # an aliphatic hydroxyl in another, which makes a series incomparable.
    ring_votes: Dict[str, Counter] = defaultdict(Counter)
    for o in outs:
        for fp, arom in o.get("rings", ()):
            ring_votes[fp][bool(arom)] += 1
    aromatic: Dict[str, bool] = {}
    ring_split: List[Tuple[str, float, int]] = []
    for fp, counts in ring_votes.items():
        verdict, n = counts.most_common(1)[0]
        total = sum(counts.values())
        aromatic[fp] = bool(verdict)
        if n < total:
            ring_split.append((fp, round(n / total, 3), total))
    if ring_split:
        outs = _sweep({"aromatic": aromatic})
        if not quiet:
            for fp, frac, total in sorted(ring_split, key=lambda x: x[1])[:6]:
                print(f"    ring {fp:<28} aromatic={aromatic[fp]} by {frac:.0%} "
                      f"of {total} observations; applied run-wide", flush=True)

    bond_counts: Dict[str, Counter] = defaultdict(Counter)
    group_n: Counter = Counter()
    triad_votes: Dict[str, Counter] = defaultdict(Counter)
    elec_votes: Dict[str, Counter] = defaultdict(Counter)
    nuc_votes: Dict[str, Dict[str, Counter]] = defaultdict(lambda: defaultdict(Counter))
    nuc_seen: Dict[str, Counter] = defaultdict(Counter)
    nuc_meta: Dict[str, Dict[str, Tuple[object, str, str]]] = defaultdict(dict)
    for o in outs:
        if not o["ok"]:
            continue
        if o["gkey"]:
            group_n[o["gkey"]] += 1
            for pr in o["pairs"]:
                bond_counts[o["gkey"]][tuple(pr)] += 1
            if o.get("elec_atom"):
                elec_votes[o["gkey"]][o["elec_atom"]] += 1
            for (name, el, cls, acc, resid, site, resid_c) in o.get("nucs", ()):
                nuc_votes[o["gkey"]][name][(cls, bool(acc))] += 1
                nuc_seen[o["gkey"]][name] += 1
                nuc_meta[o["gkey"]].setdefault(name, (resid, site, el, resid_c))
        if o["triad_key"] and o["triad_method"] == "geometric" and o["topology_ok"]:
            triad_votes[o["job"]][o["triad_key"]] += 1

    bonds: Dict[str, List[Tuple[str, str]]] = {}
    bond_stats: Dict[str, Dict[str, object]] = {}
    for gk, counts in bond_counts.items():
        n = group_n[gk]
        keep = [pr for pr, c in counts.items() if c >= max(1, math.ceil(bond_vote * n))]
        unstable = [(f"{a}-{b}", round(c / n, 3)) for (a, b), c in counts.items()
                    if c < n]
        bonds[gk] = sorted(keep)
        bond_stats[gk] = {"n_poses": n, "n_bonds": len(keep),
                          "n_unstable": len(unstable),
                          "unstable": sorted(unstable, key=lambda x: x[1])[:12]}
    triad: Dict[str, Dict[str, object]] = {}
    for job, counts in triad_votes.items():
        key, n = counts.most_common(1)[0]
        total = sum(counts.values())
        triad[job] = {"residues": [int(x) for x in key.split("|")],
                      "frac": round(n / total, 4), "n_votes": total,
                      "alternatives": {k: v for k, v in counts.items() if k != key}}
    # Every heteroatom of the ligand is classified in every surveyed pose, and the
    # majority verdict becomes the classification used everywhere. Aromaticity and
    # contact-derived decisions (a phenol pressed against the enzyme reads as
    # "tethered") are pose dependent, so without this an atom can drop in and out
    # of the nucleophile set as the conformation changes.
    nucleophiles: Dict[str, Dict[str, Dict[str, object]]] = {}
    for gk, atoms in nuc_votes.items():
        n = group_n[gk]
        inv: Dict[str, Dict[str, object]] = {}
        for name, counts in atoms.items():
            (cls, acc), top = counts.most_common(1)[0]
            total = sum(counts.values())
            resid, site, el, resid_c = nuc_meta[gk].get(name, ("", "", "", ""))
            inv[name] = {
                "element": el, "class": cls, "accepted": acc,
                "agreement": round(top / total, 4),
                "presence": round(nuc_seen[gk][name] / n, 4) if n else 0.0,
                "residue": resid, "site": site, "residue_c": resid_c,
                "alternatives": {f"{c}/{a}": v for (c, a), v in counts.items()
                                 if (c, a) != (cls, acc)},
            }
        nucleophiles[gk] = inv

    electrophile: Dict[str, Dict[str, object]] = {}
    for gk, counts in elec_votes.items():
        atom, n = counts.most_common(1)[0]
        electrophile[gk] = {"atom": atom, "frac": round(n / sum(counts.values()), 4),
                            "alternatives": {k: v for k, v in counts.items() if k != atom}}
    if not quiet:
        for job, t in sorted(triad.items()):
            e, a, h = t["residues"]
            print(f"    {job:<30} triad elbow {e} / acid {a} / His {h} "
                  f"({t['frac']:.0%} of {t['n_votes']} surveyed poses)", flush=True)
        for gk, st_ in sorted(bond_stats.items()):
            if st_["n_unstable"]:
                print(f"    {gk:<40} {st_['n_unstable']} conformation-dependent bond(s) "
                      f"resolved by vote", flush=True)
        for gk, e in sorted(electrophile.items()):
            if e["frac"] < 1.0:
                print(f"    {gk:<40} acyl carbon {e['atom']} ({e['frac']:.0%} of poses; "
                      f"outvoted {e['alternatives']})", flush=True)
        for gk, inv in sorted(nucleophiles.items()):
            unstable = {a: v for a, v in inv.items() if v["agreement"] < 1.0}
            acc = sum(1 for v in inv.values() if v["accepted"])
            print(f"    {gk:<40} {len(inv)} heteroatom(s), {acc} nucleophile(s) by vote"
                  + (f"; {len(unstable)} needed a majority: " +
                     ", ".join(f"{a} {v['class']} {v['agreement']:.0%}"
                               for a, v in sorted(unstable.items())[:6])
                     if unstable else ""), flush=True)
    return {"bonds": bonds, "triad": triad, "bond_stats": bond_stats,
            "electrophile": electrophile, "nucleophiles": nucleophiles,
            "aromatic": aromatic,
            "ring_disagreements": {fp: f for fp, f, _ in ring_split}}


def _worker(payload):
    path, opts, consensus = payload
    global BAND_NEAR_ATTACK, BAND_INTERMEDIATE
    BAND_NEAR_ATTACK, BAND_INTERMEDIATE = opts["bands"]
    try:
        res = scan_model(path, opts, consensus=consensus)
        rows, summary = candidate_rows(res, opts)
        if opts.get("write_pml") and summary and summary.get("best_nuc_atom"):
            pdir = os.path.join(str(opts["outdir"]), "pml")
            os.makedirs(pdir, exist_ok=True)
            tag = f"{res.ident['job']}_seed{res.ident.get('seed')}_sample{res.ident.get('sample')}"
            write_pml(res, str(summary["best_nuc_atom"]),
                      os.path.join(pdir, re.sub(r"[^A-Za-z0-9_.-]", "_", tag) + ".pml"))
        return {"rows": rows, "summary": summary, "ok": res.ok, "path": path,
                "error": res.error, "warnings": res.warnings,
                "inventory": inventory_rows(res, opts)}
    except Exception:
        return {"rows": [], "summary": None, "ok": False, "path": path,
                "error": traceback.format_exc(limit=4), "warnings": [],
                "inventory": []}


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="te_autoscan.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Automatic nucleophile detection and catalytic-His geometry "
                    "scanning across AlphaFold3 TE-domain models.",
        epilog="""examples:
  # scan an entire AF3 output tree, 16 processes
  python te_autoscan.py /fs04/scratch2/zf88/Anusara/af3_out -o scan_results -j 16

  # one SLURM array task out of 20
  python te_autoscan.py af3_out -o scan_results --shard $SLURM_ARRAY_TASK_ID/20
  python te_autoscan.py --merge -o scan_results

  # audit a single suspicious model, keeping rejected atoms in the output
  python te_autoscan.py job/seed-1_sample-0/model.cif -o /tmp/audit \\
      --include-rejected --classes all --write-pml

  # verify the perception logic on a built-in synthetic structure
  python te_autoscan.py --selftest
""")
    p.add_argument("inputs", nargs="*", help="AF3 output dirs, or single CIF/PDB files")
    p.add_argument("-o", "--outdir", default="te_autoscan_out")
    p.add_argument("-j", "--jobs", type=int, default=1, help="parallel processes")
    p.add_argument("--pattern", action="append", default=None,
                   help="glob for model files (repeatable). Default: *model*.cif, *.cif, *.pdb")
    p.add_argument("--rank-atom", default="CE1",
                   choices=sorted(RANK_ATOM_FIELDS), help="distance used for ranking and banding")
    p.add_argument("--rank-by", default="distance", choices=("distance", "priority"))
    p.add_argument("--classes", default="default",
                   help="nucleophile classes to keep: 'default', 'all', a group "
                        "(amine,hydroxyl,thiol,carboxyl) or explicit class names")
    p.add_argument("--bands", default=f"{BAND_NEAR_ATTACK},{BAND_INTERMEDIATE}",
                   help="near-attack and intermediate cutoffs in Angstrom")
    p.add_argument("--no-consensus", action="store_true",
                   help="disable ensemble voting; perceive bonds and the triad "
                        "independently in every pose (the old, pose-dependent behaviour)")
    p.add_argument("--survey-per-job", type=int, default=60,
                   help="poses sampled per job for the voting pass (default 60)")
    p.add_argument("--bond-vote", type=float, default=0.5,
                   help="fraction of surveyed poses a bond must appear in (default 0.5)")
    p.add_argument("--series-class", default="group", choices=("group", "fine"),
                   help="compare substrates on coarse nucleophile groups (default) or "
                        "on fine classes, which include the phenol/alcohol distinction")
    p.add_argument("--series-anchor", default="auto", choices=("auto", "c", "n"),
                   help="anchor residue numbering on the C-terminal acyl carbon, the "
                        "N-terminus, or (default) whichever lines the series up better")
    p.add_argument("--no-topology-filter", action="store_true",
                   help="allow triads that violate the elbow < acid < His sequence order")
    p.add_argument("--include-rejected", action="store_true",
                   help="also write non-nucleophilic heteroatoms with the reason")
    p.add_argument("--min-enzyme-len", type=int, default=40)
    p.add_argument("--min-ligand-atoms", type=int, default=8)
    p.add_argument("--his-window", default=None,
                   help="soft prior on the catalytic His residue number, e.g. 200,240")
    p.add_argument("--force-his", type=int, default=None,
                   help="force the catalytic His residue number")
    p.add_argument("--enzyme-chain", default=None)
    p.add_argument("--substrate-chain", default=None)
    p.add_argument("--substrate-resname", default=None)
    p.add_argument("--write-pml", action="store_true",
                   help="write a PyMOL script per model with the distances drawn")
    p.add_argument("--shard", default=None, metavar="I/N",
                   help="process only shard I of N (0- or 1-based I both accepted)")
    p.add_argument("--merge", action="store_true", help="merge shard CSVs in --outdir")
    p.add_argument("--max-models", type=int, default=None, help="stop after N models")
    p.add_argument("--selftest", action="store_true")
    p.add_argument("-q", "--quiet", action="store_true")
    p.add_argument("--version", action="version", version=f"te_autoscan {__version__}")
    return p


def merge_shards(outdir: str) -> None:
    for stem in ("nucleophile_distances_long", "per_model_summary"):
        parts = sorted(glob.glob(os.path.join(outdir, f"{stem}.shard*.csv")))
        rows: List[Dict[str, object]] = []
        for part in parts:
            with open(part) as fh:
                rows.extend(list(csv.DictReader(fh)))
        if rows:
            order = (LONG_ORDER if "long" in stem else
                     INVENTORY_ORDER if "inventory" in stem else SUMMARY_ORDER)
            write_csv(os.path.join(outdir, f"{stem}.csv"), rows, order)
            print(f"merged {len(parts)} shards -> {stem}.csv ({len(rows)} rows)")
    summ = os.path.join(outdir, "per_model_summary.csv")
    if os.path.exists(summ):
        with open(summ) as fh:
            rows = list(csv.DictReader(fh))
        opts = {"rank_atom": rows[0].get("rank_atom", "CE1") if rows else "CE1"}
        write_csv(os.path.join(outdir, "per_job_summary.csv"),
                  aggregate_jobs(rows, opts))
        print("rebuilt per_job_summary.csv")
    parts = sorted(glob.glob(os.path.join(outdir, "nucleophile_inventory_poses.shard*.csv")))
    pose_rows: List[Dict[str, object]] = []
    for part in parts:
        with open(part) as fh:
            pose_rows.extend(list(csv.DictReader(fh)))
    if pose_rows:
        inv = aggregate_inventory(pose_rows)
        write_csv(os.path.join(outdir, "nucleophile_inventory.csv"), inv, INVENTORY_ORDER)
        write_csv(os.path.join(outdir, "nucleophile_series.csv"),
                  series_rows(inv), SERIES_ORDER)
        print(f"merged {len(parts)} shards -> nucleophile_inventory.csv "
              f"({len(inv)} rows) and nucleophile_series.csv")


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.selftest:
        return run_selftest()
    os.makedirs(args.outdir, exist_ok=True)
    if args.merge:
        merge_shards(args.outdir)
        return 0
    if not args.inputs:
        build_parser().print_help()
        return 2

    bands = [float(x) for x in str(args.bands).split(",")]
    hw = tuple(int(x) for x in args.his_window.split(",")) if args.his_window else None
    opts: Dict[str, object] = {
        "rank_atom": args.rank_atom, "rank_by": args.rank_by, "classes": args.classes,
        "bands": (bands[0], bands[1]), "include_rejected": args.include_rejected,
        "min_enzyme_len": args.min_enzyme_len, "min_ligand_atoms": args.min_ligand_atoms,
        "his_window": hw, "force_his": args.force_his,
        "enzyme_chain": args.enzyme_chain, "substrate_chain": args.substrate_chain,
        "substrate_resname": args.substrate_resname, "write_pml": args.write_pml,
        "outdir": args.outdir, "no_topology_filter": args.no_topology_filter,
        "collect_bonds": False,
    }
    patterns = args.pattern or ["*model*.cif", "*.cif", "*.pdb"]
    models = discover_models(args.inputs, patterns)
    suffix = ""
    if args.shard:
        i, n = (int(x) for x in args.shard.split("/"))
        i = i - 1 if i >= n and i > 0 else i        # tolerate 1-based array ids
        models = [m for k, m in enumerate(models) if k % n == i % n]
        suffix = f".shard{i % n}"
    if args.max_models:
        models = models[: args.max_models]
    if not models:
        sys.stderr.write("no coordinate files found; check the path and --pattern\n")
        return 1
    if not args.quiet:
        print(f"te_autoscan {__version__}: {len(models)} model file(s) to scan", flush=True)

    consensus: Optional[Dict[str, object]] = None
    if not args.no_consensus:
        consensus = run_survey(models, opts, args.jobs, args.survey_per_job,
                               args.bond_vote, args.quiet)

    long_rows: List[Dict[str, object]] = []
    summaries: List[Dict[str, object]] = []
    inv_pose: List[Dict[str, object]] = []
    failures: List[Dict[str, str]] = []
    payloads = [(m, opts, consensus) for m in models]
    done = 0

    def absorb(out) -> None:
        nonlocal done
        done += 1
        long_rows.extend(out["rows"])
        inv_pose.extend(out.get("inventory", []))
        if out["summary"]:
            summaries.append(out["summary"])
        if not out["ok"]:
            failures.append({"path": out["path"], "error": out["error"][:400]})
        if not args.quiet and (done % 250 == 0 or done == len(models)):
            print(f"  {done}/{len(models)} scanned", flush=True)

    if args.jobs > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            for out in ex.map(_worker, payloads, chunksize=8):
                absorb(out)
    else:
        for pl in payloads:
            absorb(_worker(pl))

    write_csv(os.path.join(args.outdir, f"nucleophile_distances_long{suffix}.csv"),
              long_rows, LONG_ORDER)
    write_csv(os.path.join(args.outdir, f"per_model_summary{suffix}.csv"),
              summaries, SUMMARY_ORDER)
    inv = aggregate_inventory(inv_pose)
    if suffix:
        # Shards hold only part of each ensemble, so per-shard aggregates cannot be
        # concatenated. Write the per-pose rows and let --merge aggregate them once.
        write_csv(os.path.join(args.outdir, f"nucleophile_inventory_poses{suffix}.csv"),
                  inv_pose)
    else:
        write_csv(os.path.join(args.outdir, "nucleophile_inventory.csv"), inv,
                  INVENTORY_ORDER)
    jobs_rows = aggregate_jobs(summaries, opts)
    if not suffix:
        write_csv(os.path.join(args.outdir, "per_job_summary.csv"), jobs_rows)
        write_csv(os.path.join(args.outdir, "nucleophile_series.csv"),
                  series_rows(inv, args.series_anchor, args.series_class), SERIES_ORDER)
    report = {
        "version": __version__, "n_models": len(models), "n_failed": len(failures),
        "rank_atom": args.rank_atom, "bands": opts["bands"],
        "classes_kept": sorted(allowed_classes(args.classes)),
        "inputs": [os.path.abspath(x) for x in args.inputs],
        "triad_methods": dict(Counter(str(r.get("triad_method")) for r in summaries)),
        "consensus": {
            "enabled": not args.no_consensus,
            "aromaticity_disagreements": (consensus or {}).get("ring_disagreements", {}),
            "triad": (consensus or {}).get("triad", {}),
            "bond_stats": (consensus or {}).get("bond_stats", {}),
        },
        "warning_counts": dict(Counter(
            w.strip() for r in summaries for w in str(r.get("warnings", "")).split(";") if w.strip()
        ).most_common(30)),
        "failures": failures[:200],
    }
    with open(os.path.join(args.outdir, f"run_report{suffix}.json"), "w") as fh:
        json.dump(report, fh, indent=2, default=str)

    if not args.quiet:
        print(f"\n{'job':<34}{'n':>5}{'min':>8}{'med':>8}{'near%':>8}  nucleophile")
        for r in jobs_rows:
            frac = r["frac_near_attack"]
            print(f"{str(r['job'])[:33]:<34}{r['n_models']:>5}"
                  f"{str(r['min_d']):>8}{str(r['median_d']):>8}"
                  f"{(round(100*float(frac),1) if frac != '' else ''):>8}  "
                  f"{r['modal_nucleophile']} ({r['modal_nucleophile_class']}, "
                  f"agree {r['nucleophile_agreement']})")
        srows = series_rows(inv, args.series_anchor, args.series_class) if not suffix else []
        if srows and len({r["series_key"] for r in srows}) < len(srows):
            print(f"\n{'substrate':<20}{'atoms':>6}{'nucs':>6}{'rej':>5}  changes vs "
                  "largest substrate on the same enzyme")
            for r in srows:
                ch = []
                if r["n_gained"]:
                    ch.append(f"+{r['gained_classes']}")
                if r["n_lost"]:
                    ch.append(f"-{r['lost_classes']}")
                print(f"{str(r['substrate'])[:19]:<20}{r['n_substrate_atoms']:>6}"
                      f"{r['n_nucleophiles']:>6}{r['n_rejected']:>5}  "
                      + ("; ".join(ch) if ch else "(reference)"))
        if failures:
            print(f"\n{len(failures)} model(s) failed; see run_report{suffix}.json")
        print(f"\nwrote {args.outdir}/")
    return 0


# ---------------------------------------------------------------------------
# Self test: synthetic structure with a known answer
# ---------------------------------------------------------------------------

def _mkatom(name: str, element: str, x: float, y: float, z: float, b: float = 88.0):
    a = gemmi.Atom()
    a.name = name
    a.element = gemmi.Element(element)
    a.pos = gemmi.Position(x, y, z)
    a.b_iso = b
    a.occ = 1.0
    return a


def _mkres(name: str, seqid: int, atoms, het: bool = False):
    r = gemmi.Residue()
    r.name = name
    r.seqid = gemmi.SeqId(seqid, " ")
    r.het_flag = "H" if het else "A"
    for a in atoms:
        r.add_atom(a)
    return r


def build_synthetic(mutant: bool = False) -> "gemmi.Structure":
    """A minimal enzyme+ligand model whose answers are known by construction."""
    st = gemmi.Structure()
    st.name = "selftest"
    st.spacegroup_hm = "P 1"
    model = gemmi.Model("1")
    ch = gemmi.Chain("A")

    # filler polymer so the chain is recognised as the enzyme
    special = {89: "GLY", 90: "ALA", 91: ("ALA" if mutant else "SER"),
               92: "ALA", 93: "GLY", 190: "ASP", 220: "HIS"}
    for i in range(60, 271):
        name = special.get(i, "ALA")
        x, y, z = -40.0 + 0.9 * i, -30.0 + 3.0 * math.sin(i / 2.0), 3.0 * math.cos(i / 2.0)
        atoms = [_mkatom("N", "N", x - 1.2, y, z), _mkatom("CA", "C", x, y, z),
                 _mkatom("C", "C", x + 1.2, y, z), _mkatom("O", "O", x + 1.9, y + 0.9, z)]
        if name != "GLY":
            atoms.append(_mkatom("CB", "C", x, y + 1.5, z))
        ch.add_residue(_mkres(name, i, atoms))

    # real active-site geometry, overwriting the placeholders
    residues = [r for r in ch]

    his_atoms = [
        _mkatom("N", "N", -3.5, -0.6, 0.0), _mkatom("CA", "C", -2.5, 0.2, 0.0),
        _mkatom("C", "C", -3.0, 1.4, 0.0), _mkatom("O", "O", -3.6, 2.3, 0.0),
        _mkatom("CB", "C", -1.4, -0.6, 0.5),
        _mkatom("CG", "C", 0.0, 0.0, 0.0), _mkatom("ND1", "N", 1.107, -0.797, 0.0),
        _mkatom("CE1", "C", 2.245, -0.106, 0.0), _mkatom("NE2", "N", 2.010, 1.220, 0.0),
        _mkatom("CD2", "C", 0.676, 1.348, 0.0),
    ]
    if mutant:
        elbow_atoms = [
            _mkatom("N", "N", 2.0, 7.4, 1.0), _mkatom("CA", "C", 2.6, 6.4, 0.2),
            _mkatom("C", "C", 2.4, 6.8, -1.2), _mkatom("O", "O", 3.0, 7.5, -2.0),
            _mkatom("CB", "C", 2.5, 3.95, 0.3),          # CB where OG used to be
        ]
        elbow_name = "ALA"
    else:
        elbow_atoms = [
            _mkatom("N", "N", 2.0, 7.4, 1.0), _mkatom("CA", "C", 2.6, 6.4, 0.2),
            _mkatom("C", "C", 2.4, 6.8, -1.2), _mkatom("O", "O", 3.0, 7.5, -2.0),
            _mkatom("CB", "C", 2.0, 5.2, 0.6), _mkatom("OG", "O", 2.5, 3.95, 0.3),
        ]
        elbow_name = "SER"
    asp_atoms = [
        _mkatom("N", "N", -1.0, -7.8, -0.2), _mkatom("CA", "C", -0.3, -6.9, -0.5),
        _mkatom("C", "C", 0.4, -7.5, -1.7), _mkatom("O", "O", 0.6, -8.7, -1.8),
        _mkatom("CB", "C", 0.5, -5.7, -0.7), _mkatom("CG", "C", 0.0, -4.5, 0.0),
        _mkatom("OD1", "O", 0.6, -3.45, 0.3), _mkatom("OD2", "O", -1.2, -4.6, 0.3),
    ]
    new_res = {220: _mkres("HIS", 220, his_atoms),
               91: _mkres(elbow_name, 91, elbow_atoms),
               190: _mkres("ASP", 190, asp_atoms)}

    ch2 = gemmi.Chain("A")
    for r in residues:
        ch2.add_residue(new_res.get(r.seqid.num, r))

    # ---- ligand: known nucleophile 3.20 A from CE1 ----
    # Backbone is a realistic alkane zig-zag in the y-z plane; every substituent
    # is placed in that plane along the outward bisector so no 1-3 contact is
    # short enough to be mistaken for a bond.
    ce1 = gemmi.Position(2.245, -0.106, 0.0)
    u0 = (0.3008, 0.2005, 0.9326)
    bx, by, bz = (ce1.x + 3.2 * u0[0], ce1.y + 3.2 * u0[1], ce1.z + 3.2 * u0[2])
    p = [(bx, by + 0.86 * (i % 2), bz + 1.25 * i) for i in range(8)]

    def _unit(a, b):
        v = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
        n = math.sqrt(sum(c * c for c in v)) or 1.0
        return (v[0] / n, v[1] / n, v[2] / n)

    def _outward(i):
        nbrs = [j for j in (i - 1, i + 1) if 0 <= j < len(p)]
        s = [0.0, 0.0, 0.0]
        for j in nbrs:
            uu = _unit(p[i], p[j])
            s = [s[k] + uu[k] for k in range(3)]
        n = math.sqrt(sum(c * c for c in s)) or 1.0
        return (-s[0] / n, -s[1] / n, -s[2] / n)

    def _at(base, d, length):
        return (base[0] + d[0] * length, base[1] + d[1] * length,
                base[2] + d[2] * length)

    def _rot_yz(v, deg):
        a = math.radians(deg)
        return (v[0], v[1] * math.cos(a) - v[2] * math.sin(a),
                v[1] * math.sin(a) + v[2] * math.cos(a))

    o1 = _at(p[2], _outward(2), 1.23)                      # amide carbonyl O
    thio_c = _at(p[1], _outward(1), 1.52)                  # CH2 of the thiol arm
    thio_s = _at(p[1], _outward(1), 3.33)                  # SH
    hyd_c = _at(p[5], _outward(5), 1.52)                   # CH2 of the alcohol arm
    hyd_o = _at(p[5], _outward(5), 2.94)                   # OH
    v76 = _unit(p[7], p[6])
    cooh_o = _at(p[7], _rot_yz(v76, 120.0), 1.23)          # C=O
    cooh_oh = _at(p[7], _rot_yz(v76, -120.0), 1.34)        # C-OH

    lig = [
        _mkatom("N1", "N", *p[0]),          # free primary amine, the target
        _mkatom("C1", "C", *p[1]),
        _mkatom("C2", "C", *p[2]),
        _mkatom("O1", "O", *o1),            # carbonyl O -> reject
        _mkatom("N2", "N", *p[3]),          # amide N -> reject
        _mkatom("C3", "C", *p[4]),
        _mkatom("C4", "C", *p[5]),
        _mkatom("C6", "C", *p[6]),
        _mkatom("C7", "C", *p[7]),          # carboxyl C -> electrophile
        _mkatom("O2", "O", *cooh_o),
        _mkatom("O3", "O", *cooh_oh),       # carboxyl OH
        _mkatom("C5", "C", *hyd_c),
        _mkatom("O4", "O", *hyd_o),         # aliphatic hydroxyl
        _mkatom("C8", "C", *thio_c),
        _mkatom("S1", "S", *thio_s),        # thiol
    ]
    # aromatic ring with a phenol, placed away from the active site
    rc = (bx + 9.0, by - 9.0, bz + 9.0)
    for k in range(6):
        ang = math.radians(60 * k)
        lig.append(_mkatom(f"CR{k}", "C", rc[0] + 1.39 * math.cos(ang),
                           rc[1] + 1.39 * math.sin(ang), rc[2]))
    lig.append(_mkatom("O5", "O", rc[0] + 2.76, rc[1], rc[2]))       # phenol
    ligres = _mkres("LIG", 1, lig, het=True)

    ch_l = gemmi.Chain("B")
    ch_l.add_residue(ligres)
    model.add_chain(ch2)
    model.add_chain(ch_l)
    st.add_model(model)
    st.setup_entities()
    return st



def _nerf(a, b, c, bond: float, angle: float, torsion: float):
    """Place a fourth atom from three references (bond, angle, torsion)."""
    import math as _m
    ax, ay, az = a; bx, by, bz = b; cx, cy, cz = c
    v1 = (bx - ax, by - ay, bz - az)
    v2 = (cx - bx, cy - by, cz - bz)
    n2 = _m.sqrt(sum(t * t for t in v2))
    u = tuple(t / n2 for t in v2)
    cr = (v1[1] * u[2] - v1[2] * u[1], v1[2] * u[0] - v1[0] * u[2],
          v1[0] * u[1] - v1[1] * u[0])
    ncr = _m.sqrt(sum(t * t for t in cr)) or 1.0
    nvec = tuple(t / ncr for t in cr)
    m = (nvec[1] * u[2] - nvec[2] * u[1], nvec[2] * u[0] - nvec[0] * u[2],
         nvec[0] * u[1] - nvec[1] * u[0])
    ar, tr = _m.radians(angle), _m.radians(torsion)
    d = (-bond * _m.cos(ar), bond * _m.sin(ar) * _m.cos(tr), bond * _m.sin(ar) * _m.sin(tr))
    return (cx + d[0] * u[0] + d[1] * m[0] + d[2] * nvec[0],
            cy + d[0] * u[1] + d[1] * m[1] + d[2] * nvec[1],
            cz + d[0] * u[2] + d[1] * m[2] + d[2] * nvec[2])


def build_synthetic_peptide(spurious_contact: bool = False, jitter: float = 0.0,
                            depsi: bool = False, extra_acid: bool = False,
                            drop_hydroxyl: bool = False, truncate_to: int = 3,
                            clash_amine: bool = False) -> "gemmi.Structure":
    """Enzyme plus a real tripeptide ligand built from internal coordinates.

    Residue 1 has a free N-terminal amine, residue 2 a Ser-like hydroxyl,
    residue 3 a Dab-like side-chain amine, and the C-terminus is a free acid.
    With spurious_contact the terminal OXT is pushed against CB of residue 3, the
    exact artefact that makes a carboxylic acid read as an ester in one pose.
    """
    base = build_synthetic(mutant=True)
    enz = base[0]["A"]
    st = gemmi.Structure(); st.name = "selftest_peptide"; st.spacegroup_hm = "P 1"
    model = gemmi.Model("1")
    ch = gemmi.Chain("A")
    for r in enz:
        ch.add_residue(r)

    phi, psi, omega = -125.0, 135.0, 180.0
    N = [(6.0, 4.0, 6.0 + jitter)]
    CA = [(N[0][0] + 1.458, N[0][1], N[0][2])]
    C = [(CA[0][0] + 1.525 * 0.358, CA[0][1] - 1.525 * 0.934, CA[0][2])]
    O, CB, SIDE = [], [], []
    for i in range(3):
        O.append(_nerf(N[i], CA[i], C[i], 1.231, 120.5, psi + 180.0))
        CB.append(_nerf(C[i], N[i], CA[i], 1.530, 110.5, 122.0))
        if i < 2:
            n_next = _nerf(N[i], CA[i], C[i], 1.329, 116.2, psi)
            ca_next = _nerf(CA[i], C[i], n_next, 1.458, 121.7, omega)
            c_next = _nerf(C[i], n_next, ca_next, 1.525, 111.0, phi)
            N.append(n_next); CA.append(ca_next); C.append(c_next)
    og = _nerf(N[1], CA[1], CB[1], 1.417, 111.0, 62.0)       # residue 2 hydroxyl
    nz = _nerf(N[2], CA[2], CB[2], 1.489, 110.0, 65.0)       # residue 3 side amine
    oxt = _nerf(N[2], CA[2], C[2], 1.315, 118.0, psi + 180.0 + 122.0)
    if spurious_contact:
        c3, cb3 = C[2], CB[2]
        u = [cb3[k] - c3[k] for k in range(3)]
        Lu = math.sqrt(sum(t * t for t in u)) or 1.0
        u = [t / Lu for t in u]
        w = [oxt[k] - c3[k] for k in range(3)]
        dot = sum(w[k] * u[k] for k in range(3))
        perp = [w[k] - dot * u[k] for k in range(3)]
        Lp = math.sqrt(sum(t * t for t in perp)) or 1.0
        perp = [t / Lp for t in perp]
        # 1.315 A from C3 but only ~1.65 A from CB3
        ang = math.radians(37.0)
        d = [math.cos(ang) * u[k] + math.sin(ang) * perp[k] for k in range(3)]
        oxt = tuple(c3[k] + 1.315 * d[k] for k in range(3))

    if clash_amine:
        # jam the residue-3 side-chain amine onto an enzyme atom so per-pose
        # logic reads it as a covalent tether
        enz_ca = None
        for r in enz:
            if r.seqid.num == 91:
                enz_ca = next((a.pos for a in r if a.name == "CB"), None)
        if enz_ca is not None:
            nz = (enz_ca.x + 1.45, enz_ca.y, enz_ca.z)
    lig = []
    for i in range(truncate_to):
        # depsi replaces residue 3's alpha nitrogen with an ester oxygen
        if depsi and i == 2:
            lig.append(_mkatom("OE3", "O", *N[i]))
        else:
            lig.append(_mkatom(f"N{i+1}", "N", *N[i]))
        lig += [_mkatom(f"CA{i+1}", "C", *CA[i]),
                _mkatom(f"C{i+1}", "C", *C[i]), _mkatom(f"O{i+1}", "O", *O[i]),
                _mkatom(f"CB{i+1}", "C", *CB[i])]
    if not drop_hydroxyl and truncate_to >= 2:
        lig.append(_mkatom("OG2", "O", *og))
    if truncate_to >= 3:
        lig += [_mkatom("NZ3", "N", *nz)]
    lig.append(_mkatom("OXT", "O", *(oxt if truncate_to >= 3
                                     else _nerf(N[truncate_to-1], CA[truncate_to-1],
                                                C[truncate_to-1], 1.315, 118.0,
                                                psi + 302.0))))
    if extra_acid:
        cg = _nerf(N[0], CA[0], CB[0], 1.520, 113.0, 180.0)
        od1 = _nerf(CA[0], CB[0], cg, 1.250, 120.0, 0.0)
        od2 = _nerf(CA[0], CB[0], cg, 1.310, 118.0, 180.0)
        lig += [_mkatom("CG1", "C", *cg), _mkatom("OD1", "O", *od1),
                _mkatom("OD2", "O", *od2)]
    ch_l = gemmi.Chain("B")
    ch_l.add_residue(_mkres("PEP", 1, lig, het=True))
    model.add_chain(ch); model.add_chain(ch_l)
    st.add_model(model); st.setup_entities()
    return st


def _selftest_peptide(check) -> None:
    import tempfile
    opts = {"rank_atom": "CE1", "rank_by": "distance", "classes": "default",
            "bands": (3.5, 6.0), "include_rejected": True, "min_enzyme_len": 40,
            "min_ligand_atoms": 8, "his_window": None, "force_his": None,
            "enzyme_chain": None, "substrate_chain": None, "substrate_resname": None,
            "write_pml": False, "outdir": "/tmp", "no_topology_filter": False,
            "collect_bonds": True}
    with tempfile.TemporaryDirectory() as td:
        d = os.path.join(td, "pepjob", "seed-1_sample-0"); os.makedirs(d)
        path = os.path.join(d, "model.cif")
        build_synthetic_peptide().make_mmcif_document().write_file(path)
        res = scan_model(path, opts)
        byname = {c.atom_name: c for c in res.candidates}
        check("[peptide] backbone traced as 3 residues", res.backbone_units == 3,
              f"{res.backbone_units} ({res.backbone_note})")
        check("[peptide] not called cyclic", not res.backbone_cyclic, str(res.backbone_cyclic))
        check("[peptide] N1 is the N-terminal amine",
              "N1" in byname and byname["N1"].nuc_site == "n_terminal_amine"
              and byname["N1"].nuc_residue == 1,
              f'{byname.get("N1").nuc_site if "N1" in byname else "?"}'
              f'/{byname.get("N1").nuc_residue if "N1" in byname else "?"}')
        check("[peptide] OG2 is a residue-2 side chain hydroxyl",
              "OG2" in byname and byname["OG2"].nuc_class == "aliphatic_hydroxyl"
              and byname["OG2"].nuc_residue == 2 and byname["OG2"].nuc_site == "side_chain",
              f'{byname.get("OG2").nuc_class if "OG2" in byname else "?"}'
              f'/{byname.get("OG2").nuc_residue if "OG2" in byname else "?"}')
        check("[peptide] NZ3 is a residue-3 side chain amine",
              "NZ3" in byname and byname["NZ3"].nuc_class == "primary_amine"
              and byname["NZ3"].nuc_residue == 3,
              f'{byname.get("NZ3").nuc_class if "NZ3" in byname else "?"}')
        check("[peptide] backbone amide N2 rejected",
              "N2" in byname and byname["N2"].nuc_class == "amide_N",
              byname["N2"].nuc_class if "N2" in byname else "missing")
        check("[peptide] electrophile is the C-terminal acid",
              res.electrophile_atom == "C3" and "carboxyl" in res.electrophile
              and "terminus" in res.electrophile,
              f"{res.electrophile}/{res.electrophile_atom}")


def _selftest_depsi(check) -> None:
    """A depsipeptide (ester) linkage must not break the backbone trace, and a
    side-chain carboxylate must not steal the acyl assignment from the terminus."""
    import tempfile
    opts = {"rank_atom": "CE1", "rank_by": "distance", "classes": "default",
            "bands": (3.5, 6.0), "include_rejected": True, "min_enzyme_len": 40,
            "min_ligand_atoms": 8, "his_window": None, "force_his": None,
            "enzyme_chain": None, "substrate_chain": None, "substrate_resname": None,
            "write_pml": False, "outdir": "/tmp", "no_topology_filter": False,
            "collect_bonds": False}
    with tempfile.TemporaryDirectory() as td:
        for tag, kw in (("depsi", {"depsi": True}),
                        ("sidechain_acid", {"extra_acid": True})):
            d = os.path.join(td, tag, "seed-1_sample-0"); os.makedirs(d)
            path = os.path.join(d, "model.cif")
            build_synthetic_peptide(**kw).make_mmcif_document().write_file(path)
            res = scan_model(path, opts)
            if tag == "depsi":
                check("[depsi] ester linkage traced as one backbone",
                      res.backbone_units == 3 and res.backbone_fragments == 1
                      and res.backbone_branches == 0,
                      f"{res.backbone_units} residues, {res.backbone_fragments} chain(s), "
                      f"{res.backbone_branches} branch(es)")
                check("[depsi] ester link counted",
                      res.backbone_ester_links >= 1, str(res.backbone_ester_links))
                check("[depsi] coverage complete",
                      res.backbone_coverage >= 0.9, str(res.backbone_coverage))
                check("[depsi] ester oxygen not offered as a nucleophile",
                      all(not c.accepted for c in res.candidates if c.atom_name == "OE3"),
                      str([c.nuc_class for c in res.candidates if c.atom_name == "OE3"]))
            else:
                check("[sidechain acid] terminus keeps the acyl assignment",
                      res.electrophile_atom == "C3",
                      f"{res.electrophile}/{res.electrophile_atom}")


def _selftest_paths(check) -> None:
    """Path layouts seen in real AF3 runs, including ones that hid the seed."""
    cases = [
        ("<substrate>/<job>_seed<N>/<job>/<job>_model.cif",
         "/fs/out/prisanalogues/dap/pris1_seed0/pris1/pris1_model.cif",
         {"job": "pris1", "seed": 0, "substrate_dir": "dap", "collected_copy": False}),
        ("<job>/collected_models/<job>/<job>_seed000_model.cif",
         "/fs/out/ramo/pep1_on_chers/collected_models/pep1_on_chers/"
         "pep1_on_chers_seed000_model.cif",
         {"job": "pep1_on_chers", "seed": 0, "collected_copy": True}),
        ("<job>/seed-N_sample-M/model.cif",
         "/x/myjob/seed-3_sample-2/model.cif",
         {"job": "myjob", "seed": 3, "sample": 2, "collected_copy": False}),
    ]
    for label, path, want in cases:
        got = parse_identity(path)
        for k, v in want.items():
            check(f"[paths] {label} -> {k}", got.get(k) == v,
                  f"{k}={got.get(k)!r}")


def _selftest_anchor_direction(check) -> None:
    """A series truncated from the C-terminal end must be numbered from the N-end.

    The ramoshort data is shortened this way: the acyl carbon moves with every
    truncation while the N-terminal residues keep their positions. Assuming a
    C-terminal anchor scrambles the comparison, so the anchor is chosen from the
    data instead.
    """
    inv = []
    def add(sub, size, nucs):
        for atom, cls, rn, rc, site in nucs:
            inv.append({"series_key": "K", "job": f"{sub}_on_TE", "substrate": sub,
                        "variant": f"{sub}:{size}at", "n_substrate_atoms": size,
                        "atom": atom, "nuc_class": cls,
                        "class_group": CLASS_GROUP_OF.get(cls, cls), "accepted": True,
                        "residue": rn, "residue_from_c": rc, "site": site,
                        "min_d_ce1": 9.0, "class_agreement": 1.0})
    # six residues, nucleophiles at N-residues 2, 3 and 5; truncated to four residues
    add("long", 60, [("O1", "phenol", 2, 5, "side_chain"),
                     ("N1", "primary_amine", 3, 4, "side_chain"),
                     ("O2", "aliphatic_hydroxyl", 5, 2, "side_chain")])
    add("short", 40, [("O1", "phenol", 2, 3, "side_chain"),
                      ("N1", "primary_amine", 3, 2, "side_chain")])
    got = {r["substrate"]: r for r in series_rows(inv, "auto", "group")}
    check("[anchor] C-terminal truncation is detected and numbered from the N-end",
          all(r["anchor"] == "residue_from_n" for r in got.values()),
          str({k: v["anchor"] for k, v in got.items()}))
    check("[anchor] with the right anchor the positional diff agrees",
          got["short"]["positional_diff_agrees"],
          f"gained={got['short']['gained_positional']} "
          f"lost={got['short']['lost_positional']}")
    wrong = {r["substrate"]: r for r in series_rows(inv, "c", "group")}
    check("[anchor] forcing the wrong anchor visibly disagrees",
          not wrong["short"]["positional_diff_agrees"],
          str(wrong["short"]["gained_positional"]))
    check("[anchor] the class diff is right under either anchor",
          got["short"]["lost_classes"] == wrong["short"]["lost_classes"] == "hydroxyl:1",
          f"{got['short']['lost_classes']} / {wrong['short']['lost_classes']}")


def _selftest_series(check) -> None:
    """A substitution and truncation series on one enzyme.

    Three substrates: the full tripeptide, the same with the residue-2 hydroxyl
    replaced by nothing (substitution removes a nucleophile), and a dipeptide
    truncation. One pose of the full substrate has its side-chain amine jammed
    against the enzyme so that per-pose logic calls it a tether and drops it.
    The ensemble must keep it, and the series table must show what each variant
    gained or lost.
    """
    import tempfile
    opts = {"rank_atom": "CE1", "rank_by": "distance", "classes": "default",
            "bands": (3.5, 6.0), "include_rejected": False, "min_enzyme_len": 40,
            "min_ligand_atoms": 8, "his_window": None, "force_his": None,
            "enzyme_chain": None, "substrate_chain": None, "substrate_resname": None,
            "write_pml": False, "outdir": "/tmp", "no_topology_filter": False,
            "collect_bonds": False}
    with tempfile.TemporaryDirectory() as td:
        specs = [("full", {}), ("no_oh", {"drop_hydroxyl": True}),
                 ("trunc", {"truncate_to": 2})]
        for name, kw in specs:
            for k in range(6):
                d = os.path.join(td, f"{name}_on_TE", f"seed-{k}_sample-0")
                os.makedirs(d)
                st = build_synthetic_peptide(jitter=0.04 * k,
                                             clash_amine=(name == "full" and k == 3), **kw)
                st.make_mmcif_document().write_file(os.path.join(d, "model.cif"))
        models = discover_models([td], ["*model*.cif"])
        cons = run_survey(models, opts, 1, 60, 0.5, quiet=True)
        clash = os.path.join(td, "full_on_TE", "seed-3_sample-0", "model.cif")
        alone = scan_model(clash, opts)
        fixed = scan_model(clash, opts, consensus=cons)
        a_nz = next((c for c in alone.candidates if c.atom_name == "NZ3"), None)
        f_nz = next((c for c in fixed.candidates if c.atom_name == "NZ3"), None)
        check("[series] clashing pose alone drops the side-chain amine",
              a_nz is not None and not a_nz.accepted,
              f"{a_nz.nuc_class}/{a_nz.accepted}" if a_nz else "missing")
        check("[series] ensemble restores it",
              f_nz is not None and f_nz.accepted and f_nz.nuc_class == "primary_amine",
              f"{f_nz.nuc_class}/{f_nz.accepted}" if f_nz else "missing")
        check("[series] the disagreement is recorded, not hidden",
              f_nz is not None and f_nz.pose_class == "tethered"
              and f_nz.class_agreement is not None and f_nz.class_agreement < 1.0,
              f"{f_nz.pose_class}/{f_nz.class_agreement}" if f_nz else "missing")

        inv_pose, per_sub = [], {}
        for m in models:
            r = scan_model(m, opts, consensus=cons)
            inv_pose.extend(inventory_rows(r, opts))
        inv = aggregate_inventory(inv_pose)
        for r in inv:
            per_sub.setdefault(str(r["substrate"]), []).append(r)
        acc = {k: {r["atom"] for r in v if r["accepted"]} for k, v in per_sub.items()}
        check("[series] rejected heteroatoms still listed in the inventory",
              all(len(v) > len(acc[k]) for k, v in per_sub.items()),
              str({k: (len(v), len(acc[k])) for k, v in per_sub.items()}))
        check("[series] every substrate keeps all three nucleophiles it should",
              acc.get("PEP", set()) >= {"N1", "NZ3"},
              str(acc))
        srows = {r["variant"]: r for r in series_rows(inv)}
        check("[series] three substrates compared on one enzyme",
              len(srows) == 3 and len({r["series_key"] for r in srows.values()}) == 1,
              str({k: v["series_key"] for k, v in srows.items()}))
        full = max(srows.values(), key=lambda r: int(r["n_substrate_atoms"]))
        smaller = [r for r in srows.values()
                   if int(r["n_substrate_atoms"]) < int(full["n_substrate_atoms"])]
        check("[series] smaller variants report the nucleophiles they lost",
              all(r["n_lost"] >= 1 for r in smaller),
              str([(r["substrate"], r["n_lost"], r["lost_classes"]) for r in smaller]))
        check("[series] a subset never reports a gain",
              all(r["n_gained"] == 0 for r in smaller),
              str([(r["substrate"], r["n_gained"], r["gained_classes"])
                   for r in smaller]))
        check("[series] the anchor is chosen from the data, not assumed",
              all(r["anchor_selected"] == "auto" for r in srows.values())
              and len({r["anchor"] for r in srows.values()}) == 1,
              str({k: (v["anchor"], v["anchor_selected"]) for k, v in srows.items()}))
        fixed_c = series_rows(inv, "c", "group")
        fixed_n = series_rows(inv, "n", "group")
        auto = series_rows(inv, "auto", "group")
        def hits(rr):
            return sum(1 for r in rr if r["positional_diff_agrees"])
        check("[series] auto is at least as good as either fixed anchor",
              hits(auto) >= max(hits(fixed_c), hits(fixed_n)),
              f"auto={hits(auto)} c={hits(fixed_c)} n={hits(fixed_n)}")
        check("[series] the lost nucleophile is named by site and class, not atom name",
              all("res" in r["lost_positional"] and ":" in r["lost_positional"]
                  for r in smaller),
              str([r["lost_positional"] for r in smaller]))


def _selftest_aromaticity(check) -> None:
    """The same phenol ring must not read as an alcohol in a different variant.

    Two jobs hold the same ligand; in one of them the ring is slightly puckered so
    the geometric test fails. Per-job voting cannot see the disagreement because it
    is consistent within each job. The run-wide ring vote must.
    """
    import tempfile
    opts = {"rank_atom": "CE1", "rank_by": "distance", "classes": "default",
            "bands": (3.5, 6.0), "include_rejected": True, "min_enzyme_len": 40,
            "min_ligand_atoms": 8, "his_window": None, "force_his": None,
            "enzyme_chain": None, "substrate_chain": None, "substrate_resname": None,
            "write_pml": False, "outdir": "/tmp", "no_topology_filter": False,
            "collect_bonds": False}
    with tempfile.TemporaryDirectory() as td:
        for job, pucker in (("flat_job", 0.0), ("puckered_job", 0.16)):
            for k in range(6):
                d = os.path.join(td, job, f"seed-{k}_sample-0"); os.makedirs(d)
                st = build_synthetic(mutant=True)
                if pucker:
                    for r in st[0]["B"]:
                        for n, a in enumerate(r):
                            if a.name.startswith("CR") or a.name == "O5":
                                a.pos = gemmi.Position(
                                    a.pos.x, a.pos.y,
                                    a.pos.z + (pucker if n % 2 else -pucker))
                st.make_mmcif_document().write_file(os.path.join(d, "model.cif"))
        models = discover_models([td], ["*model*.cif"])
        flat = [m for m in models if "flat_job" in m][0]
        puck = [m for m in models if "puckered_job" in m][0]
        a = next(c for c in scan_model(flat, opts).candidates if c.atom_name == "O5")
        b = next(c for c in scan_model(puck, opts).candidates if c.atom_name == "O5")
        check("[aromaticity] pose geometry alone splits the two jobs",
              a.nuc_class != b.nuc_class, f"{a.nuc_class} vs {b.nuc_class}")
        cons = run_survey(models, opts, 1, 60, 0.5, quiet=True)
        check("[aromaticity] the ring disagreement is detected run-wide",
              bool(cons.get("ring_disagreements")), str(cons.get("ring_disagreements")))
        a2 = next(c for c in scan_model(flat, opts, consensus=cons).candidates
                  if c.atom_name == "O5")
        b2 = next(c for c in scan_model(puck, opts, consensus=cons).candidates
                  if c.atom_name == "O5")
        check("[aromaticity] both jobs agree once the ring vote is applied",
              a2.nuc_class == b2.nuc_class, f"{a2.nuc_class} vs {b2.nuc_class}")
        check("[aromaticity] the coarse group was stable throughout",
              a.class_group == b.class_group == "hydroxyl",
              f"{a.class_group}/{b.class_group}")


def _selftest_consensus(check) -> None:
    """Nine clean poses plus one with a squeezed terminal contact.

    Per-pose perception reads the odd pose as an ester and loses the acyl carbon.
    The ensemble vote restores it, which is the failure seen in the real data.
    """
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        root = os.path.join(td, "consensus_job")
        for k in range(10):
            d = os.path.join(root, f"seed-{k}_sample-0"); os.makedirs(d)
            st = build_synthetic_peptide(spurious_contact=(k == 7), jitter=0.03 * k)
            st.make_mmcif_document().write_file(os.path.join(d, "model.cif"))
        bad = os.path.join(root, "seed-7_sample-0", "model.cif")
        opts = {"rank_atom": "CE1", "rank_by": "distance", "classes": "default",
                "bands": (3.5, 6.0), "include_rejected": False, "min_enzyme_len": 40,
                "min_ligand_atoms": 8, "his_window": None, "force_his": None,
                "enzyme_chain": None, "substrate_chain": None,
                "substrate_resname": None, "write_pml": False, "outdir": td,
                "no_topology_filter": False, "collect_bonds": False}
        alone = scan_model(bad, opts)
        check("[consensus] squeezed pose alone loses the acid",
              alone.electrophile != "free_carboxyl_C_at_backbone_terminus",
              alone.electrophile)
        models = discover_models([root], ["*model*.cif", "*.cif"])
        cons = run_survey(models, opts, 1, 60, 0.5, quiet=True)
        gk = list(cons["bonds"])[0]
        check("[consensus] the spurious bond was outvoted",
              cons["bond_stats"][gk]["n_unstable"] >= 1
              and not any({"CB3", "OXT"} == {a, b} for a, b in cons["bonds"][gk]),
              str(cons["bond_stats"][gk]["unstable"][:3]))
        fixed = scan_model(bad, opts, consensus=cons)
        check("[consensus] same pose with the ensemble vote recovers the acid",
              fixed.electrophile.startswith("free_carboxyl_C_at_backbone_terminus")
              and fixed.electrophile_atom == "C3",
              f"{fixed.electrophile}/{fixed.electrophile_atom}")
        check("[consensus] acyl carbon was itself voted on",
              cons["electrophile"][gk]["atom"] == "C3"
              and cons["electrophile"][gk]["frac"] >= 0.9,
              str(cons["electrophile"].get(gk)))
        check("[consensus] triad agreed across the ensemble",
              cons["triad"].get("consensus_job", {}).get("frac") == 1.0
              and cons["triad"]["consensus_job"]["residues"] == [91, 190, 220],
              str(cons["triad"].get("consensus_job")))
        check("[consensus] forced triad reused in every pose",
              fixed.triad.forced and fixed.triad.method == "ensemble_consensus",
              fixed.triad.method)


def _selftest_topology(check) -> None:
    """A decoy Ser-His-Asp with the acid after the His must be rejected."""
    import tempfile
    st = build_synthetic(mutant=True)
    ch = gemmi.Chain("A")
    for r in st[0]["A"]:
        ch.add_residue(r)
    st2 = gemmi.Structure(); st2.spacegroup_hm = "P 1"
    m = gemmi.Model("1"); m.add_chain(ch)
    for c in st[0]:
        if c.name != "A":
            m.add_chain(c)
    st2.add_model(m); st2.setup_entities()
    t_ok = detect_triad(st2, {"A"}, None, require_topology=True)
    t_off = detect_triad(st2, {"A"}, None, require_topology=False)
    check("[topology] canonical triad passes the filter",
          t_ok.found and t_ok.topology_ok and t_ok.his_seqid == 220,
          f"{t_ok.method}/{t_ok.his_seqid}")
    check("[topology] residue key reported",
          t_ok.residue_key == "91|190|220", t_ok.residue_key)
    check("[topology] filter off still finds the same site",
          t_off.his_seqid == 220, str(t_off.his_seqid))

def run_selftest() -> int:
    import tempfile
    ok = True
    checks: List[Tuple[str, bool, str]] = []

    def check(label: str, cond: bool, detail: str = "") -> None:
        nonlocal ok
        ok = ok and bool(cond)
        checks.append((label, bool(cond), detail))

    for mutant in (False, True):
        st = build_synthetic(mutant=mutant)
        with tempfile.TemporaryDirectory() as td:
            job = os.path.join(td, "selftest_job", "seed-1_sample-0")
            os.makedirs(job)
            path = os.path.join(job, "model.cif")
            st.make_mmcif_document().write_file(path)
            with open(os.path.join(os.path.dirname(job), "ranking_scores.csv"), "w") as fh:
                fh.write("seed,sample,ranking_score\n1,0,0.77\n")
            opts = {"rank_atom": "CE1", "rank_by": "distance", "classes": "default",
                    "bands": (3.5, 6.0), "include_rejected": True,
                    "min_enzyme_len": 40, "min_ligand_atoms": 8, "his_window": None,
                    "force_his": None, "enzyme_chain": None, "substrate_chain": None,
                    "substrate_resname": None, "write_pml": False, "outdir": td}
            res = scan_model(path, opts)
            tag = "mutant" if mutant else "wildtype"
            check(f"[{tag}] model parsed", res.ok, res.error[:200])
            t = res.triad
            check(f"[{tag}] triad found geometrically", t.method == "geometric", t.method)
            check(f"[{tag}] catalytic His = 220", t.his_seqid == 220, str(t.his_seqid))
            check(f"[{tag}] acid = ASP190", t.acid_seqid == 190, str(t.acid_seqid))
            check(f"[{tag}] topology elbow<acid<His", t.topology_ok, str(t.topology_ok))
            check(f"[{tag}] elbow = residue 91", t.elbow_seqid == 91, str(t.elbow_seqid))
            check(f"[{tag}] elbow atom = {'CB' if mutant else 'OG'}",
                  t.elbow_atom == ("CB" if mutant else "OG"), t.elbow_atom)
            check(f"[{tag}] mutant flag", t.mutant == mutant, str(t.mutant))
            check(f"[{tag}] His-elbow distance ~2.79",
                  t.d_his_elbow is not None and abs(t.d_his_elbow - 2.79) < 0.05,
                  str(t.d_his_elbow))
            check(f"[{tag}] GxSxG motif recognised", t.motif_ok, t.motif)
            check(f"[{tag}] ring orientation canonical", t.orientation == "canonical",
                  t.orientation)

            byname = {c.atom_name: c for c in res.candidates}
            check(f"[{tag}] N1 = primary amine",
                  byname.get("N1") is not None and byname["N1"].nuc_class == "primary_amine",
                  byname["N1"].nuc_class if "N1" in byname else "missing")
            check(f"[{tag}] N1 distance to CE1 = 3.20",
                  "N1" in byname and abs(byname["N1"].d_ce1 - 3.20) < 0.02,
                  str(byname.get("N1").d_ce1 if "N1" in byname else None))
            check(f"[{tag}] N1 banded near_attack",
                  "N1" in byname and byname["N1"].band == "near_attack",
                  byname["N1"].band if "N1" in byname else "")
            check(f"[{tag}] N2 rejected as amide",
                  byname.get("N2") is not None and byname["N2"].nuc_class == "amide_N"
                  and not byname["N2"].accepted,
                  byname["N2"].nuc_class if "N2" in byname else "missing")
            check(f"[{tag}] O1 rejected as carbonyl",
                  "O1" in byname and not byname["O1"].accepted,
                  byname["O1"].nuc_class if "O1" in byname else "missing")
            check(f"[{tag}] O4 = aliphatic hydroxyl",
                  "O4" in byname and byname["O4"].nuc_class == "aliphatic_hydroxyl",
                  byname["O4"].nuc_class if "O4" in byname else "missing")
            check(f"[{tag}] O5 = phenol",
                  "O5" in byname and byname["O5"].nuc_class == "phenol",
                  byname["O5"].nuc_class if "O5" in byname else "missing")
            check(f"[{tag}] O3 = carboxyl OH",
                  "O3" in byname and byname["O3"].nuc_class == "carboxyl_OH",
                  byname["O3"].nuc_class if "O3" in byname else "missing")
            check(f"[{tag}] S1 = thiol",
                  "S1" in byname and byname["S1"].nuc_class == "thiol",
                  byname["S1"].nuc_class if "S1" in byname else "missing")
            check(f"[{tag}] electrophile = free carboxyl carbon",
                  res.electrophile.startswith("free_carboxyl_C")
                  and res.electrophile_atom == "C7",
                  f"{res.electrophile}/{res.electrophile_atom}")
            check(f"[{tag}] backbone traced", res.backbone_units >= 1,
                  f"{res.backbone_units} units, {res.backbone_note}")
            rows, summary = candidate_rows(res, opts)
            check(f"[{tag}] best nucleophile = N1",
                  summary.get("best_nuc_atom") == "N1", str(summary.get("best_nuc_atom")))
            check(f"[{tag}] ranking score read from sidecar",
                  str(summary.get("ranking_score")) == "0.77",
                  str(summary.get("ranking_score")))
            check(f"[{tag}] seed/sample parsed",
                  res.ident["seed"] == 1 and res.ident["sample"] == 0,
                  f"{res.ident['seed']}/{res.ident['sample']}")

    _selftest_peptide(check)
    _selftest_depsi(check)
    _selftest_paths(check)
    _selftest_aromaticity(check)
    _selftest_anchor_direction(check)
    _selftest_series(check)
    _selftest_consensus(check)
    _selftest_topology(check)

    width = max(len(c[0]) for c in checks) + 2
    for label, passed, detail in checks:
        mark = "PASS" if passed else "FAIL"
        extra = "" if passed else f"   got: {detail}"
        print(f"{label:<{width}} {mark}{extra}")
    print(f"\n{sum(1 for c in checks if c[1])}/{len(checks)} checks passed")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
