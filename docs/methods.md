# Methods

Each section names the script that implements it. Reference keys refer to
[`manuscript/references.bib`](../manuscript/references.bib); every DOI was
verified against Crossref (`scripts/python/verify_references.py`).

## 1. Cohort and data sources (steps 01-02)

The CPTAC-3 pancreatic ductal adenocarcinoma (PDAC) discovery cohort
[cao2021] was chosen because it combines RNA-seq and TMT proteomics from the
same cryopulverised tissue, has normal-adjacent tissue (NAT) for part of the
patients, and has mature overall-survival follow-up (129 of 170 GDC cases had
died at download). Only open-access data were used.

* **GDC** [grossman2016], project `CPTAC-3`, primary site `pancreas`: case
  records with demographics, diagnoses, pathology details, treatments and
  follow-ups; RNA-seq *STAR - Counts* (GENCODE v36) [dobin2013]; masked
  somatic mutations from the *Aliquot Ensemble Somatic Variant Merging and
  Masking* workflow (WXS); gene-level copy number from *AscatNGS* (WGS)
  [vanloo2010]. Files were downloaded through the GDC API and md5-checked.
* **PDC**, study `PDC000270` (TMT11, 25 plexes): biospecimen table, TMT
  experimental design, and the gene-level log2-ratio matrix of the Common Data
  Analysis Pipeline (CDAP) [rudnick2016], i.e. each sample relative to the
  pooled reference channel of its plex.

## 2. Sample matching and database (step 03)

GDC and PDC share CPTAC case (`C3L-…`/`C3N-…`) and sample (`C3L-…-01`)
identifiers but not aliquots. Matching was therefore done at the sample level
and recorded per case and tissue type (tumour / normal):

| level | definition |
|---|---|
| `same_sample` | the proteome aliquot comes from exactly the RNA-seq tissue piece |
| `pooled_sample` | the proteome aliquot is a pool of 2-3 tumour pieces that includes the RNA-seq piece |
| `same_case` | same patient and tissue type, different piece (none occurred) |

Tumour RNA was never paired with normal protein. Reference pools, QC cell
lines and the nine normal-duct donors of the study (PDC disease type "Other")
were excluded, as were cases not diagnosed as PDAC in GDC. All tables (cases,
follow-ups, treatments, GDC files, PDC aliquots, aliquot-to-sample links,
matches with inclusion flags) are stored in an SQLite database with a fixed
schema (`db/schema.sql`); analysis results are added in step 14.

## 3. Matrices (step 04)

* **RNA**: STAR counts of the reverse-strand column; strandedness was
  determined per library from the ratio of second- to first-strand counts
  (all 213 libraries reverse-stranded).
* **Protein**: CDAP log2 ratios, median-centred per sample. Proteins were
  analysed when quantified in >= 50% of the samples of each compared group.
* **Mutations**: protein-altering variant classes (missense, nonsense,
  frameshift, in-frame indel, splice site, start/stop loss) from one tumour
  MAF per case (the WXS sample matching the RNA-seq sample when available).
* **Copy number**: ASCAT absolute gene-level copy number; log2(CN + 0.5) /
  (median CN of the sample + 0.5) for correlation analyses; homozygous
  deletion = absolute copy number 0.

## 4. Quality control and purity (step 05)

* **RNA**: DESeq2 size factors and variance-stabilising transformation
  (blind) on protein-coding genes with >= 10 counts in >= 10% of samples
  [love2014]; PCA on the 2,000 most variable genes.
* **Tissue identity**: Whipple resections include duodenum. For each RNA
  sample the mean VST of 17 enterocyte/goblet/Paneth markers (e.g. *SI*,
  *APOA4*, *DEFA5*, *FABP2*, *CDX2*) and of 10 acinar markers (e.g. *PRSS1*,
  *CPA1*, *PNLIP*) was computed. Samples closer to the duodenal cluster than
  to the median of normal pancreas (cut-off half-way between the two centres)
  were flagged and excluded: 10 normal samples (all RNA-only) and one tumour.
* **Protein**: quantification rate, PCA on proteins quantified in all
  samples; the variance of each PC explained by tissue and by TMT plex
  (adjusted R^2), to decide whether plex must be modelled.
* **Purity**: ESTIMATE stromal and immune scores from RNA [yoshihara2013]
  (tidyestimate implementation); neoplastic cellularity as 2 x KRAS variant
  allele fraction (heterozygous, diploid assumption) in KRAS-mutant tumours.
