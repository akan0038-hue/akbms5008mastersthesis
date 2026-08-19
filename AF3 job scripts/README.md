# AF3 job scripts

Batch scripts that generate AlphaFold3 (AF3) input JSONs and, in most
cases, submit them as SLURM jobs — the pre-AF3-run half of the pipeline.
See the root README for how this fits into the overall workflow.

---

## Prerequisites

Every script here calls two lab-internal scripts that are **not
distributed in this repository**:

* **`af3_te_substrate.py`** — generates the AF3 input JSON and CCD ligand
  definition from a protein FASTA file, a substrate SMILES string, and a
  covalent modification site (chain index + residue number).
* **`run_af3_ppant_pipeline.py`** — submits the generated JSONs to SLURM
  (only called by `ser_batch.py` and `ramoplanin_casestudy_batch.py`; see
  below).

Both must be on your `$PYTHONPATH` or in the working directory. Run these
scripts from the repository root so the relative `FASTA/`, `json_dir/`,
and `out_dir/` paths resolve correctly.

---

## Scripts

### `ramoplanin_casestudy_batch.py`

Builds the full cognate/non-cognate cross-screen for the ramoplanin
family (Section 3.2 of the thesis, "AlphaFold3 does not generate
near-attack geometry for ramoplanin-family pairings" — 1,800 structures).
Each of four substrates — ramoplanin, enduracidin, chersinamycin, and a
synthetic construct — is modelled against every non-cognate TE domain in
the family. Runs both steps: JSON generation via `af3_te_substrate.py`,
then job submission via `run_af3_ppant_pipeline.py` (25 seeds, GPU
partition, 5 hour limit, 154G memory).

### `dab_batch.py`, `dap_batch.py`, `dthr_batch.py`, `lthr_batch.py`, `ser_batch.py`

The five position-2 substitution variants from Section 3.5 ("Substitution
at position 2 alters attack geometry according to nucleophile
chemistry" — 4,500 structures total): Dab (2,4-diaminobutyric acid), Dap
(2,3-diaminopropionic acid), D-Thr, L-Thr, and Ser. Each script applies
one substituent's SMILES fragment to the same six TE domain FASTA files
(`FASTA/pris/`: `pris1`, `pris2`, `vir1`–`vir4`) at the same modification
site.

`dab_batch.py`, `dap_batch.py`, `dthr_batch.py`, and `lthr_batch.py` only
run the JSON-generation step — submit the resulting `json_dir/` to
`run_af3_ppant_pipeline.py` yourself. `ser_batch.py` runs both steps.

---

## A note on paths before you push

`json_dir` and `out_dir` output paths, and the `FASTA` input directory,
are relative to wherever you run the script from (repo root, by
convention). Nothing else is hardcoded — but check for local machine
paths before committing if you customise these.
