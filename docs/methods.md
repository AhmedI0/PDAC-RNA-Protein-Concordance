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

## 5. mRNA-protein concordance (step 07)

On matched, QC-pass samples (136 tumours, 21 NATs) and 9,254 genes with both
layers (protein quantified in >= 50% of tumours):

* gene-wise Spearman correlation across tumours (BH-FDR [benjamini1995]),
  across NATs, and across tumours + NATs pooled;
* null distribution from 20 permutations of tumour labels in the protein
  matrix;
* partial Spearman correlation removing ESTIMATE stromal and immune scores and
  the acinar marker score (rank-based residuals), to separate co-variation due
  to tissue composition from gene-intrinsic coupling;
* determinants: correlation of gene-wise rho with mRNA mean, mRNA SD, protein
  SD and protein missingness;
* sample-wise correlation of gene-centred (z-scored) profiles per tumour, and
  its association with KRAS-based cellularity;
* gene-set level: two-sided Mann-Whitney test of gene-wise rho inside vs
  outside each Hallmark [liberzon2015], Reactome [gillespie2022] and selected
  GO cellular-component complex set.

## 6. Genomic alterations (step 08)

For genes with a copy-number change (|log2 ratio| > 0.2) in >= 10 tumours,
Spearman correlation of copy number with the gene's own mRNA and protein;
attenuation = rho(CN, mRNA) - rho(CN, protein) [goncalves2017]; complex
subunits (members of GO CC "complex" sets with <= 400 genes) compared with
other genes. For *KRAS*, *TP53*, *CDKN2A* and *SMAD4*, mRNA and protein were
compared between alteration classes (missense/in-frame, truncating,
homozygous deletion) and wild type (Mann-Whitney, BH-FDR).

## 7. Differential expression (step 09)

* **mRNA**: DESeq2 Wald test, design `~ tissue`, all QC-pass PDAC RNA samples
  (160 tumours, 42 NATs); log2 fold changes shrunk with apeglm [zhu2019].
* **Protein**: limma [ritchie2015], design `~ tissue + plex`, patients as
  random block via `duplicateCorrelation` [smyth2005] (some tumours and NATs
  come from the same patient), robust empirical Bayes with intensity trend;
  136 tumours, 66 NATs.
* Significance: FDR < 5% and |log2FC| > 1 (mRNA) or > 0.5 (protein; TMT
  ratios are compressed by co-isolation interference). Genes were classified
  as concordant up/down, mRNA-only, protein-only or opposite. A sensitivity
  analysis repeated both models on the identical matched samples.

## 8. Pathways (step 10)

* Pre-ranked GSEA (fgsea multilevel [korotkevich2016]) on the DESeq2 Wald
  statistic and the limma moderated t, Hallmark and Reactome sets (15-500
  genes). As a sensitivity analysis, Reactome sets in which >= 25% of genes are
  ribosomal, proteasomal, spliceosomal or respiratory-chain subunits were
  excluded before comparing normalised enrichment scores between layers.
* GSVA [hanzelmann2013] per-tumour Hallmark scores from mRNA and from protein
  over the same 9,254 genes (missing protein values set to the gene minimum
  - 0.5 for this step only). Pathway-level concordance = Spearman correlation
  of the mRNA and protein score of a pathway across tumours, compared with the
  median gene-level rho of its members (paired Wilcoxon).

## 9. Subtypes (step 11)

* **PurIST** [rashid2020]: 8 weighted top-scoring gene pairs on log2(TPM+1);
  probability >= 0.5 = basal-like. Model taken unchanged from the authors'
  `pdacR` package (commit `5b044b2`).
* **Moffitt** tumour-specific signatures [moffitt2015] (25 basal-like, 25
  classical genes; renamed symbols mapped through org.Hs.eg.db aliases):
  continuous score = mean z(basal) - mean z(classical) per layer; categorical
  call per layer by consensus clustering (k = 2, hierarchical, Pearson
  distance, 1,000 resamplings of 80% of tumours) [wilkerson2010].
  Agreement: Cohen's kappa [cohen1960]; AUC of each score against PurIST.
* **Unsupervised**: consensus k-means on the 1,500 most variable features of
  each layer (features complete in all tumours), k = 2-6, k chosen by minimal
  proportion of ambiguous clustering (PAC) [senbabaoglu2014]; layer agreement
  by adjusted Rand index [hubert1985] at every k.

## 10. Survival (steps 12-13)

* Overall survival from diagnosis, administratively censored at 5 years;
  133 matched tumours, 105 deaths. Cox models adjusted for age (per 10 years),
  stage (I-II / III-IV / unknown), resection margin (R0 / R1-R2 / RX-unknown)
  and grade (G1-G2 / G3-G4); proportional hazards checked with Schoenfeld
  residuals [grambsch1994].
* Gene-level Cox per SD of each gene in each layer, adjusted for stage.
* Prognostic models: elastic-net Cox (alpha = 0.5) [friedman2010; simon2011]
  on the 2,000 most variable genes per layer (z-scored), with clinical
  covariates unpenalised where included; 50 x repeated 5-fold CV with the
  penalty tuned by an inner 5-fold CV in each training set; Harrell's C on
  held-out folds [harrell1996]. Signature genes = selected in >= 50% of the
  250 training-fold fits. Out-of-fold risk groups = median split of each
  tumour's mean out-of-fold risk rank.

## 11. Software

R 4.5.2 with Bioconductor 3.22 and CRAN packages frozen at the 15 Dec 2025
Posit Package Manager snapshot (`renv.lock`); Python 3.12 (`envs/environment.yml`).
Random seeds are fixed in `config/config.yaml`.
