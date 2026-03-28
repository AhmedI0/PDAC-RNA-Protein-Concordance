# Data

Nothing in `raw/`, `interim/` or `processed/` is committed: every file is
downloaded or rebuilt by the pipeline (`python run_pipeline.py`). Total size
after a full run is about 1.7 GB, almost all of it the raw GDC downloads.

| Folder | Written by | Contents |
|---|---|---|
| `raw/gdc/` | step 01 | GDC case records (`cases.json`), file metadata (`files_*.json`), STAR gene counts (`rna/`), masked somatic MAFs (`maf/`), ASCAT gene-level copy number (`cn/`); every file md5-checked |
| `raw/pdc/` | step 02 | PDC000270 study record, biospecimens, TMT11 design, CDAP gene-level log2-ratio matrix |
| `raw/external/` | step 00 | `pdacR` gene lists and the PurIST model, pinned to commit `5b044b2` |
| `interim/` | step 03 | flattened clinical table, GDC file table, PDC aliquot table, RNA-protein sample matches |
| `processed/` | steps 04-05, 10 | analysis matrices (see `processed/README.md`) |

## Sources

* **GDC** (https://portal.gdc.cancer.gov), project `CPTAC-3`, primary site
  `pancreas`, open-access files only: *Gene Expression Quantification*
  (`STAR - Counts`, GENCODE v36), *Masked Somatic Mutation* (`Aliquot Ensemble
  Somatic Variant Merging and Masking`, WXS) and *Gene Level Copy Number*
  (`AscatNGS`, WGS). The GDC data release that was current at download time is
  recorded in `raw/gdc/gdc_status.json`.
* **PDC** (https://pdc.cancer.gov), study `PDC000270` (*CPTAC PDA Discovery
  Study - Proteome*, TMT11, 25 plexes). Gene-level quantification from the PDC
  Common Data Analysis Pipeline (CDAP): log2 ratio of each sample to the
  pooled reference channel of its plex.

No controlled-access data are used. Use of CPTAC data is subject to the NCI
data-use policies of the GDC and PDC; please cite Cao et al., *Cell* 2021
(doi:10.1016/j.cell.2021.08.023) when using this cohort.
