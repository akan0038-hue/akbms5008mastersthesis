# 5008 Masters Thesis by Anusara Kannangara

Code accompanying the Masters thesis *BMS5008 — Understanding Thioesterase
Domain Involved in Peptide Antibiotic Cyclisation*, investigating whether
AlphaFold3 (AF3) predicts near-attack conformations for thioesterase (TE)
domain-mediated macrocyclisation of non-ribosomal peptide substrates,
across the ramoplanin and pristinamycin/virginiamycin case study systems.

Full methodological detail is in the thesis, Sections 2.2–2.5 and
Appendices A–C. This README covers how the code fits together; each
folder has its own README with the specifics for that stage.

---

## Pipeline overview

```
FASTA/ (TE domain sequences) + substrate SMILES
        │
        ▼
  AF3 job scripts/        ── AF3 input JSON generation + SLURM submission
        │                    (via private af3_te_substrate.py /
        │                     run_af3_ppant_pipeline.py — see below)
        ▼
  AlphaFold3 (external, not in this repo)
        │
        ▼
  Analysis/te_autoscan.py ── active-site geometry extraction
        │                    (Appendix A, Methods 2.2.4–2.2.5)
        ▼
  Analysis/R/te_analysis.R  ── statistics + Figures 3.1–3.6 (Appendix B)
  Analysis/R/lcms_traces.R  ── LC-MS parsing + Figure 3.7, Appendix F
  Analysis/pymol/*.py       ── structural panel renders (e.g. Fig. 3.4.1)
```

Each stage can be re-run on its own output — you don't need to
regenerate AF3 predictions to re-run the statistics, for example.

---

## Repository structure

```
.
├── README.md
├── AF3 job scripts/              # AF3 input generation + job submission
│   ├── README.md
│   ├── ramoplanin_casestudy_batch.py
│   ├── dab_batch.py
│   ├── dap_batch.py
│   ├── dthr_batch.py
│   ├── lthr_batch.py
│   └── ser_batch.py
├── FASTA/                        # TE domain sequence inputs
│   ├── README.md
│   ├── ramo/
│   │   ├── ramo_te199.fasta
│   │   ├── endura_te_199a.fasta
│   │   ├── che_nohisortag.fa
│   │   └── virginiamycin.fa
│   └── pris/
│       ├── pris1.fa
│       ├── pris2.fa
│       ├── vir1.fa
│       ├── vir2.fa
│       ├── vir3.fa
│       └── vir4.fa
└── Analysis/                     # post-AF3 geometry scan, stats, figures
    ├── README.md
    ├── requirements.txt
    ├── te_autoscan.py
    ├── R/
    │   ├── te_analysis.R
    │   ├── auroc_ci.R
    │   └── lcms_traces.R
    └── pymol/
        └── chapter3_panels.py
```

---

## External dependencies (not in this repository)

The scripts in `AF3 job scripts/` call two lab-internal scripts that are
**not distributed in this public repo**:

* **`af3_te_substrate.py`** — builds the AF3 input JSON and CCD ligand
  definition for a given TE domain sequence, substrate SMILES, and
  covalent linkage.
* **`run_af3_ppant_pipeline.py`** — submits the generated JSONs as SLURM
  jobs and manages the AF3 run.

Both must be present on your `$PYTHONPATH` or in the working directory
before running anything in `AF3 job scripts/`. See that folder's README
for how each script uses them.

---

## Requirements

```bash
cd Analysis
pip install -r requirements.txt      # te_autoscan.py: gemmi (required), numpy (recommended)
```

```r
install.packages(c("tidyverse", "rstatix", "patchwork", "pROC", "scales"))
install.packages("devEMF")   # optional, vector figures for Word
```

PyMOL is required to run `Analysis/pymol/chapter3_panels.py`.

---

## A note on paths before you push

Several scripts hardcode local filesystem paths near the top (`BASE`,
`ROOT`, `XLSX`, `F_*`, `json_dir`, `out_dir`) — that's intentional, since
each is meant to be edited to point at your own copy of the data. Swap
the checked-in placeholders for your own paths before running, and strip
any local machine paths or personal account details before pushing.

---

## Code and data availability

Deposited at github.com/akan0038-hue/akbms5008mastersthesis, per Appendix
H of the thesis. Data are available from the author on request.

## License

MIT

## Author

Anusara Kannangara, Cryle Laboratory, Monash University.
