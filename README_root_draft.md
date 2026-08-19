# akbms5008mastersthesis

Code accompanying the Masters thesis *BMS5008 — Understanding Thioesterase Domain Involved in Peptide Antibiotic Cyclisation*, investigating
whether AlphaFold3 (AF3) predicts near-attack conformations for
thioesterase (TE) domain-mediated macrocyclisation of non-ribosomal
peptide substrates, across the ramoplanin and pristinamycin/virginiamycin
case study systems.

Full methodological detail is in the thesis, Sections 2.2–2.5 and
Appendices A–C. This README covers how the code fits together and how to
run it.

---

## Pipeline overview

```
substrate SMILES + TE domain FASTA
        │
        ▼
  af3_input_prep.py   ── AF3 input JSON + CCD ligand definition
        │                (Appendix C, Methods 2.2.3)
        ▼
  AlphaFold3 (external, not in this repo)
        │
        ▼
  te_autoscan.py       ── active-site geometry extraction
        │                (Appendix A, Methods 2.2.4–2.2.5)
        ▼
  R/te_analysis.R       ── statistics + Figures 3.1–3.6 (Appendix B)
  R/lcms_traces.R        ── LC-MS parsing + Figure 3.7, Appendix F
  pymol/*.pml             ── structural panel renders (e.g. Fig. 3.4.1)
```

Each stage is independent and can be re-run on its own output — you don't
need to regenerate AF3 predictions to re-run the statistics, for example.

---

## Repository structure

```
.
├── README.md
├── requirements.txt
├── te_autoscan.py            # geometry extraction (Appendix A)
├── R/
│   ├── te_analysis.R         # stats + Figures 3.1–3.6 (Appendix B)
│   ├── auroc_ci.R            # AUROC + 95% CI, Table 3.7
│   └── lcms_traces.R         # LC-MS parsing, Figure 3.7 & Appendix F
├── pymol/
│   └── chapter3_panels.py    # structural figure panels
├── batch_scans/              # te_autoscan batch runs, one per substrate/case study
│   ├── dab_batch.py
│   ├── dap_batch.py
│   ├── dthr_batch.py
│   ├── lthr_batch.py
│   ├── ser_batch.py
│   └── ramoplanin_casestudy_batch.py
└── af3_screening/            # pre-AF3 input prep and cross-screening
    ├── af3_input_prep.py     # AF3 JSON + CCD generation (Appendix C)
    ├── cross_screen.py       # pairwise substrate × TE domain screen
    ├── run_batch_modifications.py
    └── FASTA/
        ├── ramo_te199.fasta
        ├── endura_te_199a.fasta
        ├── che_nohisortag.fa
        └── pris/              # pristinamycin case study substrates
```

---

## What each stage does

### 1. Input preparation — `af3_screening/`

`af3_input_prep.py` takes a TE domain sequence, a substrate SMILES string,
and the covalent linkage definition (ligand acyl-carbon atom index +
catalytic serine residue number, modelled as alanine to block AF3 from
resolving the native covalent bond), and writes:

- the AF3 input JSON
- a Chemical Component Dictionary (CCD) entry for the ligand

`cross_screen.py` builds the full cognate/non-cognate matrix — each
substrate modelled in its own domain and cross-modelled into every other
domain in the same family. `run_batch_modifications.py` runs this across
a set of position-variant substrates (e.g. the position-2 series in
Section 3.5).

Substrate conformers themselves (RDKit distance-geometry embedding +
force-field minimisation, one conformer per random seed) are prepared
separately — see Methods 2.2.2.

### 2. Structure prediction — AlphaFold3 (external)

Not part of this repository. AF3 v3.0.0 (parameter set
AlphaFold-beta-20231127), 50 random seeds per domain/substrate pairing
(25 with 5 diffusion samples, 25 with 1), giving 150 structures per
pairing. See Methods 2.2.3 and Appendix G for the full structure
inventory.

### 3. Geometry extraction — `te_autoscan.py`

Parses AF3 output with no per-substrate configuration required:

- classifies enzyme vs. substrate chains (≥40-residue polymer chains =
  enzyme)
- locates the catalytic Ser(→Ala)-His-Asp triad geometrically, not by
  residue number
- perceives candidate nucleophiles on the substrate from molecular
  connectivity alone (free amines, aliphatic/phenolic hydroxyls, thiols
  accepted; amides, guanidinium, esters, thioesters excluded)
- measures nucleophile→His NE2 distance, nucleophile→carbonyl carbon
  distance, and the Bürgi–Dunitz attack angle
- applies the near-attack conformation criteria: ≤3.5 Å to NE2, ≤4.0 Å to
  the electrophilic carbon, and 95–115° attack angle, all three
  simultaneously (Methods 2.2.5)

Bond perception and active-site assignment are voted across sampled poses
per ensemble rather than re-derived per structure, since these calls are
sensitive to individual conformation.

The `batch_scans/` scripts run this over each case-study dataset
(ramoplanin family, truncation series, pristinamycin native, position-2
variants).

### 4. Statistics and figures — `R/`

`te_analysis.R` runs every statistical test in Chapter 3 and produces
Figures 3.1–3.6: Kruskal–Wallis with Dunn's test for 3+ group
comparisons, Mann–Whitney with rank-biserial effect size for two-group
contrasts, Spearman correlation for the length-vs-distance relationship,
and Fisher's exact test for near-attack counts (Methods 2.5). Four input
CSV paths and one output directory are set at the top of the script — no
other configuration needed. Run with `source("te_analysis.R")` in
RStudio (not `Rscript` from the terminal — see the script header).

`auroc_ci.R` computes AUROC with DeLong 95% CIs for Table 3.7.

`lcms_traces.R` parses the Shimadzu LabSolutions Excel export and
produces Figure 3.7 and the supplementary traces in Appendix F.

### 5. Structural figures — `pymol/`

`chapter3_panels.py` renders each structural panel twice from an
identical camera — a clean geometry-only PNG and a labelled "guide" PNG —
so labels can be typeset externally without PyMOL's label collisions.
Also writes the exact measured distance/angle strings used in each panel.

---

## Requirements

```
pip install -r requirements.txt          # te_autoscan.py: gemmi (required), numpy (recommended)
```

```r
install.packages(c("tidyverse", "rstatix", "patchwork", "pROC", "scales"))
install.packages("devEMF")   # optional, vector figures for Word
```

PyMOL is required to run `pymol/chapter3_panels.py`.

---

## A note on paths before you push

Several scripts hardcode local filesystem paths near the top (`BASE`,
`ROOT`, `XLSX`, `F_*`) — that's intentional, since each is meant to be
edited to point at your own copy of the data. Swap the checked-in
placeholders for your own paths before running, and strip any local
machine paths or personal account details before pushing.

---

## Code and data availability

Deposited at github.com/akan0038-hue/akbms5008mastersthesis, per Appendix
H of the thesis. Data are available from the author on request.

## License

MIT

## Author

Anusara Kannangara, Cryle Laboratory, Monash University.
