# PDAC-RNA-Protein-Concordance

**Where do mRNA and protein agree in pancreatic cancer, where do they
disagree, and does it matter for subtypes, pathways and survival?**

An independent, fully reproducible re-analysis of the CPTAC-3 pancreatic
ductal adenocarcinoma (PDAC) discovery cohort (Cao et al., *Cell* 2021),
integrating RNA-seq, TMT proteomics, somatic mutations, copy number and
overall survival from the NCI Genomic Data Commons (GDC) and Proteomic Data
Commons (PDC). Every sample link between the two portals is made explicitly
and stored in an SQLite database.

`R` · `Python` · `Bioconductor` · `DESeq2` · `limma` · `fgsea` · `GSVA` ·
`ConsensusClusterPlus` · `survival` · `glmnet` · `SQLite`

---

## Relation to published work

The cohort was published by the CPTAC consortium in Cao et al., *Cell*
184:5031 (2021), which reported mRNA-protein correlations, subtypes and
phosphoproteomics. Later studies re-used it for pan-cancer proteogenomics
and PDAC subtyping. This repository is an **independent re-analysis** from the
public files: it does not use the consortium's harmonised tables, makes the
GDC-PDC sample matching explicit, adds a tissue-identity QC, quantifies
composition effects on concordance, and tests (with cross-validation)
whether the proteome adds prognostic value. Results agree with the original
where they overlap (median gene-wise rho 0.30 in tumours and 0.23 in NATs here, vs 0.36 and 0.26 reported). See
[`manuscript/manuscript.md`](manuscript/manuscript.md) for the full write-up
and verified references.

## Repository layout

```
config/config.yaml        cohort IDs, paths, thresholds, seeds
src/cancerpg/             Python package: GDC/PDC clients, matching, stats, plotting
R/                        R helpers (data access, figure theme, survival coding)
scripts/python/, scripts/R/   pipeline steps 00-14
run_pipeline.py           runs all steps in order, logs to logs/
db/schema.sql             SQLite schema (database built by steps 03 and 14)
resources/                pinned PDAC subtype signatures and PurIST model
results/tables/           all result tables (TSV)
results/figures/          all figures (SVG/PNG/JPG)
manuscript/               manuscript, references (DOIs verified)
docs/                     methods, figure guide, data provenance
tests/                    pytest + testthat
envs/, renv.lock          exact software environment
data/                     downloaded/derived data (not committed; see data/README.md)
```

## Reproduce

Requirements: R 4.5.x, conda, ~3 GB disk, internet access to GDC/PDC.

```bash
conda env create -p ./.cache/env -f envs/environment.yml
conda activate ./.cache/env
Rscript envs/setup_r_env.R          # or: Rscript -e "renv::restore()"
python run_pipeline.py              # all steps, ~1 h
python run_pipeline.py --list       # list steps
python -m pytest tests/python -q
```

Set `RSCRIPT` to the full path of `Rscript` if it is not on `PATH`. Software
versions are frozen (R 4.5.2, Bioconductor 3.22, CRAN snapshot 2025-12-15,
Python 3.12); random seeds are in `config/config.yaml`.

## Data and licence

Data: CPTAC-3 open-access files from the GDC and PDC, used under the NCI data
use policies; please cite Cao et al. (2021) when using the cohort. Code: MIT
licence (see `LICENSE`).

## Author

Muhammad Ahmed Mushtaq
