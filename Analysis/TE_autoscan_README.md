[README.md](https://github.com/user-attachments/files/31183909/README.md)
# te_autoscan

Hands-off nucleophile detection and catalytic-His geometry scanning for
AlphaFold3 models of NRPS thioesterase (TE) domain complexes, plus the
downstream R analysis and PyMOL figure scripts built on top of it.

```
.
├── te_autoscan.py     # the scan itself (this section of the README)
├── requirements.txt
├── R/
│   ├── te_analysis.R  # all Chapter 3 figures + statistics from te_autoscan output
│   ├── auroc_ci.R      # AUROC + 95% CI for Table 3.7
│   └── lcms_traces.R   # LC-MS chromatogram parsing/plotting (independent of te_autoscan)
└── pymol/
    └── chapter3_panels.py  # clean + labelled-guide structural panel renders
```

Point it at a folder of AlphaFold3 output and it will, with no per-job
input:

1. **Walk the AF3 output tree** and find every predicted model (`.cif`),
   pairing each with its seed/sample identity and confidence sidecars.
2. **Classify chains** into enzyme (polymer) and substrate (ligand / short
   peptide) without being told which is which.
3. **Locate the Ser-His-Asp catalytic triad** geometrically, tolerating a
   Ser→Ala active-site mutation (using the Ala CB as the Ser OG proxy) and
   histidine ring-flip ambiguity.
4. **Perceive every chemically plausible nucleophile** on the substrate
   from connectivity alone, so arbitrary CCD atom naming (`C1`/`N7`/`O12`…)
   doesn't matter. Amides, guanidines, phenol vs. alcohol, thioester
   sulfur, esters, and carboxylate oxygens are all distinguished and
   either scored or excluded.
5. **Measure nucleophile → His CE1 distance** (plus NE2, ND1, flip-invariant
   minimum, elbow atom, and attack distance / Bürgi–Dunitz angle to the
   electrophilic carbon), and assign near-attack / intermediate / too-far
   bands.
6. **Emit tidy CSVs**, a JSON run report, and optional PyMOL scripts.

Bond perception prefers the authoritative `_chem_comp_bond` block if the
CIF carries one, and otherwise falls back to covalent-radius distance
criteria. Nothing is keyed on residue or atom names for the substrate;
enzyme side chains do use standard PDB naming, which AF3 always writes
correctly. Every heuristic that fires is recorded in the output, so a
suspicious row can be traced back without rerunning.

## Requirements

- Python 3.9+
- [gemmi](https://gemmi.readthedocs.io/) (required)
- numpy (optional, but recommended)

```bash
pip install -r requirements.txt
```

## Usage

```bash
# scan an entire AF3 output tree, 16 parallel processes
python te_autoscan.py /path/to/af3_out -o scan_results -j 16

# one SLURM array task out of 20, then merge
python te_autoscan.py af3_out -o scan_results --shard $SLURM_ARRAY_TASK_ID/20
python te_autoscan.py --merge -o scan_results

# audit a single suspicious model, keeping rejected atoms in the output
python te_autoscan.py job/seed-1_sample-0/model.cif -o /tmp/audit \
    --include-rejected --classes all --write-pml

# verify the perception logic against a built-in synthetic structure
python te_autoscan.py --selftest
```

### Key options

| Flag | Default | Description |
|---|---|---|
| `inputs` | — | AF3 output directories, or single CIF/PDB files |
| `-o`, `--outdir` | `te_autoscan_out` | Output directory |
| `-j`, `--jobs` | `1` | Parallel processes |
| `--pattern` | *(auto)* | Glob for model files (repeatable). Default: `*model*.cif`, `*.cif`, `*.pdb` |
| `--rank-atom` | `CE1` | Distance used for ranking and banding (`CE1`, `NE2`, `ND1`, …) |
| `--rank-by` | `distance` | Rank candidates by `distance` or chemical `priority` |
| `--classes` | `default` | Nucleophile classes to keep: `default`, `all`, a group (`amine`, `hydroxyl`, `thiol`, `carboxyl`), or explicit class names |
| `--bands` | `3.5,6.0` | Near-attack and intermediate cutoffs (Å) |
| `--no-consensus` | off | Disable ensemble voting; perceive bonds/triad independently per pose |
| `--survey-per-job` | `60` | Poses sampled per job for the voting pass |
| `--bond-vote` | `0.5` | Fraction of surveyed poses a bond must appear in to be kept |
| `--series-class` | `group` | Compare substrates on coarse groups or `fine` classes (incl. phenol/alcohol) |
| `--series-anchor` | `auto` | Anchor residue numbering on the C-terminal acyl carbon, the N-terminus, or whichever lines the series up better |
| `--no-topology-filter` | off | Allow triads that violate the elbow < acid < His sequence order |
| `--include-rejected` | off | Also write non-nucleophilic heteroatoms with the rejection reason |
| `--min-enzyme-len` | `40` | Minimum residues for a chain to be classified as enzyme |
| `--min-ligand-atoms` | `8` | Minimum atoms for a chain to be classified as substrate |
| `--his-window` | *(none)* | Soft prior on the catalytic His residue number, e.g. `200,240` |
| `--force-his` | *(none)* | Force the catalytic His residue number |
| `--enzyme-chain` / `--substrate-chain` / `--substrate-resname` | *(none)* | Override auto-classification |
| `--write-pml` | off | Write a PyMOL script per model with the measured distances drawn |
| `--shard` | *(none)* | Process only shard `I/N` of the input set (SLURM array friendly) |
| `--merge` | off | Merge shard CSVs already written to `--outdir` |
| `--max-models` | *(none)* | Stop after N models |
| `--selftest` | off | Run the built-in self-test suite on a synthetic structure |
| `-q`, `--quiet` | off | Suppress progress output |
| `--version` | — | Print version and exit |

Run `python te_autoscan.py --help` for the full, current list.

## Output

Written to `--outdir`:

| File | Contents |
|---|---|
| `nucleophile_distances_long.csv` | One row per candidate nucleophile atom per model (long format) |
| `per_model_summary.csv` | One row per model, with the best-scoring nucleophile and triad |
| `per_job_summary.csv` | One row per job (aggregated across seeds/samples) |
| `nucleophile_inventory.csv` | Ensemble-voted inventory of nucleophile sites per substrate |
| `nucleophile_series.csv` | Substrate-to-substrate comparison of nucleophile classes gained/lost |
| `*.pml` (optional, `--write-pml`) | PyMOL script per model highlighting measured distances |
| a JSON run report | Run parameters, timing, and per-model status/errors |

Every accepted or rejected atom carries the heuristic(s) that decided it,
so any row can be traced back to a specific geometric or connectivity rule
without rerunning the scan.

## Design notes

- Bond perception prefers the `_chem_comp_bond` block in the CIF when
  present; otherwise it falls back to covalent-radius distance criteria
  with a valence cap, so a compressed or clashing pose can't invent
  chemically impossible bonds.
- The substrate side is perceived purely from connectivity — no residue
  or atom naming assumptions — so arbitrary CCD atom names don't matter.
  The enzyme side does rely on standard PDB atom naming, which AF3 writes
  correctly.
- Ensemble consensus (on by default) perceives bonds and the catalytic
  triad by voting across sampled poses within a job, so a single spurious
  contact in one seed/sample doesn't flip the call for the whole job.
  Disable with `--no-consensus` for the old, pose-independent behaviour.
- Ser→Ala active-site mutants are handled explicitly: the Ala Cβ is used
  as the Ser Oγ proxy for the elbow position.

## Testing

The script ships with a self-contained test suite that builds a synthetic
peptide/depsipeptide structure and checks perception, triad detection,
topology filtering, ensemble consensus, and series comparison end to end
— no external files needed:

```bash
python te_autoscan.py --selftest
```

## Downstream analysis

`te_autoscan.py` is the front of the pipeline. Everything in `R/` and
`pymol/` consumes its output (`*nucleophile_distances_long.csv` and the
CIF models) to make the figures, statistics, and panels for the thesis
chapter.

### `R/te_analysis.R`

Every figure and statistic for the AlphaFold3 results chapter, in one
script: reads the four `te_autoscan` long CSVs (one per substrate
dataset), restricts each to the ring-closing nucleophile, applies the
near-attack criterion (NE2 distance, attack distance, and Bürgi–Dunitz
angle all within threshold at once), and produces the confidence
comparisons, per-domain and per-pairing boxplots, heatmaps, the
position-2 variant comparison, and the AUROC-vs-confidence analysis —
each backed by the matching non-parametric test (Kruskal–Wallis,
Mann–Whitney with rank-biserial effect size, Dunn's test, Fisher's exact
test). Writes figures as PDF/PNG (and EMF if `devEMF` is installed) and
every statistic to `stats/all_statistics.csv`.

Edit the `BASE`/`F_*` paths near the top to point at your own CSVs, then
source the whole file in RStudio (do **not** run it with `Rscript` from
the terminal — see the header comment for why). Requires:

```r
install.packages(c("tidyverse", "rstatix", "patchwork", "pROC", "scales"))
install.packages("devEMF")  # optional, vector figures for Word
```

### `R/auroc_ci.R`

Standalone script for Table 3.7: AUROC with DeLong 95% confidence
intervals for each confidence metric (nucleophile pLDDT, pTM, ipTM,
ranking score) against the near-attack outcome, pooled and per-domain.
Point `df` at the same data frame `te_analysis.R` computes `roc()` on.
Requires `pROC`, `dplyr`, `readr`, `purrr`, `tibble`.

### `R/lcms_traces.R`

Parses a Shimadzu LabSolutions Excel export (extracted ion chromatograms)
and plots retention-time traces per sample, with three normalisation
modes (per-trace, per-sample, raw). Includes the target/deletion-ion
figures used in Chapter 3 and a simple peak-picking table. Set `XLSX`
and `OUT` near the top, then source in RStudio. Requires:

```r
install.packages(c("tidyverse", "readxl", "patchwork"))
```

### `pymol/chapter3_panels.py`

PyMOL script that renders each structural figure panel twice from an
identical camera: a clean geometry-only PNG and a labelled "guide" PNG
at the same crop, so labels can be typeset in Illustrator/Inkscape
without PyMOL's label collisions. Also writes `measurements.txt` with
the exact distance/angle strings to type onto each panel. Edit the
`ROOT`/`PANELS` list at the top (input CIF folder, per-panel atom
selections, zoom, expected distances) then run inside PyMOL:

```
run chapter3_panels.py
```

### A note on paths before you push

All four scripts hardcode local filesystem paths near the top (`BASE`,
`ROOT`, `XLSX`, `F_*`, or a bare `read_csv(...)`) — that's intentional,
since each one is meant to be edited to point at your own copy of the
data. The paths currently checked in are generic placeholders
(`/path/to/your/...`); swap them for your own before running, and if you
fork this for someone else, make sure no local machine paths or personal
account details creep back in.

## License

This project is licensed under the MIT License

## Author

Written by Anusara Kannangara for the Cryle Lab TE cyclisation pipeline.
