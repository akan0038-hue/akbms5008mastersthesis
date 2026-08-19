# FASTA

TE domain sequence inputs used by the scripts in `AF3 job scripts/`.
Purification tags are removed from all sequences to match the expression
constructs described in Methods 2.4 (thesis).

---

## `ramo/` — ramoplanin-family domains

Used by `ramoplanin_casestudy_batch.py` for the cognate/non-cognate
cross-screen (Section 3.2).

| File | Domain |
|---|---|
| `ramoplanin.fasta` | Ramoplanin TE domain |
| `enduracidin.fasta` | Enduracidin TE domain |
| `chersinamycin.fa` | Chersinamycin TE domain |

## `pris/` — pristinamycin/virginiamycin domains

Used by `dab_batch.py`, `dap_batch.py`, `dthr_batch.py`, `lthr_batch.py`,
and `ser_batch.py` for the position-2 substitution series (Section 3.5)
and by the native pristinamycin dataset (Section 3.4).

| File | Domain |
|---|---|
| `pris1.fa`, `pris2.fa` | Pristinamycin TE domains |
| `vir1.fa`–`vir4.fa` | Virginiamycin TE domains |

Two pristinamycin and four virginiamycin domains — six domains total,
matching the six-domain streptogramin dataset referenced throughout
Chapter 3.
