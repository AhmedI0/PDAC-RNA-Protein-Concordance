# Agreement and disagreement between mRNA and protein in pancreatic ductal adenocarcinoma: an independent proteogenomic re-analysis of the CPTAC-3 cohort

**Muhammad Ahmed Mushtaq**

*This is an independent re-analysis of public data. The cohort was generated
and first described by the Clinical Proteomic Tumor Analysis Consortium (Cao
et al., 2021). Every number below is produced by the code in this repository;
methods are given in full in [`docs/methods.md`](../docs/methods.md).*

---

## Abstract

Pancreatic ductal adenocarcinoma (PDAC) remains one of the most lethal solid
tumours, and almost all of its molecular classification is based on mRNA.
Proteins carry out the functions that mRNA only encodes, so the extent to
which the transcriptome reports the proteome decides how far transcriptomic
subtypes, pathways and biomarkers can be trusted. RNA-seq, TMT proteomics,
somatic mutation, copy-number and survival data of the CPTAC-3 PDAC cohort
were retrieved from the NCI Genomic and Proteomic Data Commons and linked at
the level of the tissue sample. In total 136 tumours and 21 normal-adjacent
tissues (NATs) had both layers from the same tissue; 10 NATs and one tumour
were excluded because their transcriptome was that of duodenal mucosa. The
median gene-wise mRNA-protein correlation was 0.30 in tumours, and 68% of
9,254 genes were significantly coupled. Coupling was strongest for
extracellular matrix, interferon and immunoglobulin genes and weakest for the
subunits of ribosomes, respiratory chain, spliceosome and proteasome.
Copy-number effects on mRNA reached the protein for only 18% of genes, and
least for complex subunits. The
basal-like/classical axis was nearly identical in both layers (rho 0.96) and
equally prognostic when measured by protein (adjusted HR 1.70). Genes with
concordant mRNA and protein carried the same prognostic information in both
layers. However, genome-wide elastic-net models from mRNA, protein or both did
not outperform clinical pathology in 50 × 5-fold cross-validation, while a
single pre-specified subtype score did. Protein measurement therefore adds
mechanistic rather than prognostic information in resected PDAC, and
pre-specified, biologically grounded scores should be preferred over
data-driven signatures in cohorts of this size.

---

## 1. Introduction

Pancreatic cancer has the poorest outcome of all common cancers. Its 5-year
relative survival is 13% and it is the third leading cause of cancer death in
the United States (Siegel et al., 2024). It is projected to become the second
leading cause of cancer-related death by 2030 (Rahib et al., 2014). Only a
minority of patients present with resectable disease. Even after resection
and adjuvant modified FOLFIRINOX, median overall survival is 54.4 months
(Conroy et al., 2018).

Molecular subtyping has been the main route to stratify this disease. Bulk
transcriptomes first separated classical, quasi-mesenchymal and exocrine-like
tumours (Collisson et al., 2011). Virtual microdissection later separated
tumour-intrinsic classical and basal-like subtypes from stromal signals, with
a shorter survival for basal-like tumours (Moffitt et al., 2015). An
integrated genomic analysis described squamous, pancreatic progenitor,
immunogenic and ADEX classes (Bailey et al., 2016). Basal-like tumours respond
less often to first-line modified FOLFIRINOX (Aung et al., 2018), and a
single-sample classifier, PurIST, was built to make the call usable in the
clinic (Rashid et al., 2020). All of these classifications are defined by
mRNA.

mRNA is an imperfect proxy for protein. In cultured cells, protein levels
depend as much on translation and degradation rates as on mRNA abundance
(Schwanhäusser et al., 2011). Across samples, the dependency of protein on
mRNA is moderate and gene-specific (Liu et al., 2016). Copy-number gains of
genes encoding complex subunits are attenuated at the protein level
(Gonçalves et al., 2017). Pan-cancer CPTAC analyses show that mRNA-protein
correlation differs between tumour types, pathways and drug targets (Arad
and Geiger, 2023), and that patients with weaker post-transcriptional control
have shorter overall survival in several cancer types (Ghoshdastider and
Sendoel, 2023). Pan-cancer proteogenomics has also widened the set of
candidate therapeutic targets beyond what mRNA predicts (Savage et al., 2024).

The CPTAC-3 PDAC cohort provides the material to ask these questions in
pancreatic cancer. The tumours were cryopulverised and the same tissue powder
was used for sequencing and mass spectrometry (Cao et al., 2021). The data
are distributed through two portals, the Genomic Data Commons (Grossman et
al., 2016) and the Proteomic Data Commons, whose proteomes are processed by a
common pipeline (Rudnick et al., 2016). Harmonised tables exist (Li et al.,
2023), although the link between an RNA-seq file and a proteome aliquot, and
the tissue identity of normal samples, are not explicit in them.

If transcriptomic subtypes and signatures are applied to PDAC without asking
which of them survive translation, biomarkers will continue to be selected
on a layer that may not report the functional state of the tumour. To our
knowledge, the relationship between mRNA-protein agreement and the
prognostic value of each layer has not been tested with cross-validated
models in PDAC. Therefore, this study is designed to re-analyse the CPTAC-3
PDAC cohort from the public files, with explicit sample matching, and to
determine where mRNA and protein agree, why they disagree, and whether
protein adds information for subtypes, pathways and survival.

**Objectives**

* to link RNA-seq and proteome data at the tissue-sample level and check the
  tissue identity of every sample;
* to quantify gene-, sample- and pathway-level mRNA-protein concordance and
  its genomic and compositional determinants;
* to compare tumour-versus-normal changes, pathways and molecular subtypes
  between layers;
* to assess the prognostic value of each layer with pre-specified and
  cross-validated models.

---

## 2. Results

### 2.1 Cohort, sample matching and tissue identity

The GDC listed 170 CPTAC-3 pancreas cases, and the PDC study PDC000270
contained 250 patient aliquots in 25 TMT11 plexes (Figure 1). RNA-seq and
proteome were linked through the shared sample identifiers. For 80 tumours
the proteome aliquot was cut from the same tissue piece as the RNA-seq
aliquot. For 57 tumours it was a pool of two or three pieces that included
the RNA-seq piece. No match had to be made across different pieces of the
same patient. Twenty-one NATs were matched in the same way.

An RNA principal component analysis showed a separate group of normal
samples. These samples expressed enterocyte, goblet-cell and Paneth-cell
genes (*SI*, *APOA4*, *FABP2*, *DEFA5*, *CDX2*). A score of 17 intestinal
markers separated them from pancreas with a wide gap (centres 6.5 and 15.1;
Figure 2A-B). Ten NATs, all without proteome, and one tumour (C3L-04479) were
flagged as duodenal mucosa and excluded. The final matched set comprised 136
tumours and 21 NATs.

The protein data were dominated by tissue type. The first principal component
explained 39.3% of variance and 68% of it was explained by tumour versus
normal (Figure 2C-D). TMT plex explained up to 31% of the variance of minor
components, so plex was modelled as a covariate. The ESTIMATE immune and
stromal score (Yoshihara et al., 2013) correlated negatively with the
cellularity derived from *KRAS* allele fraction (rho = −0.49, n = 131,
p = 2 × 10⁻⁹).

Of the 136 patients, 51% were male and the median age was 65 years (Table 1).
Most tumours were stage II (52%), grade 2 (71%) and in the pancreatic head
(74%). Overall survival was available for 133 patients. Median overall
survival was 19.7 months (95% CI 13.6-23.6), with 105 deaths within 5 years.
*KRAS*, *TP53*, *CDKN2A* and *SMAD4* were mutated in 82%, 66%, 18% and 14% of
tumours.

### 2.2 mRNA and protein are moderately coupled, and coupling is gene-class specific

The gene-wise Spearman correlation between mRNA and protein across the 136
tumours had a median of 0.30 for 9,254 genes (Figure 3A). The same analysis
with randomly paired tumours gave a median of 0.00 and a 95th percentile of
0.15. In total 68% of genes were positively correlated at FDR < 5%, and 0.6%
negatively. The median was 0.23 in the 21 NATs (Figure 3B).

Correlation increased with the variability of a gene's mRNA across tumours
(rho = 0.49 between gene-wise correlation and mRNA standard deviation) and
decreased with the fraction of tumours in which the protein was missing
(rho = −0.16; Figure 3C). Removing the ESTIMATE stromal and immune scores and
the acinar marker score lowered the median only to 0.26.

Coupling differed sharply between gene classes (Figure 3D-E). The highest
Hallmark medians were for epithelial-mesenchymal transition (0.60), apical
surface (0.57) and the interferon-gamma (0.56) and interferon-alpha (0.55)
responses. Immunoglobulin (0.63), collagen trimer (0.57) and MHC class I
(0.53) complexes were also highly coupled. The lowest were MYC targets
(0.21), oxidative phosphorylation (0.22), the cytosolic ribosome (0.08), the
respiratory chain (0.14), the spliceosome (0.17) and the proteasome (0.21).
The interferon and EMT sets remained at the top after adjustment for
composition.

The sample-wise correlation of each tumour's centred mRNA and protein
profiles had a median of 0.34. It increased with *KRAS*-derived cellularity
(rho = 0.40, p = 1.7 × 10⁻⁵; Figure 3F). Tumours whose proteome came from a
pool of pieces did not differ from same-piece tumours (median 0.33 versus
0.35, p = 0.50).

### 2.3 Copy-number dosage is buffered at the protein level, and drivers act on different layers

Copy number was correlated with the gene's own mRNA and protein in 135
tumours for 6,529 genes with recurrent changes (Figure 7A). The median
correlation was 0.25 with mRNA and 0.10 with protein. Copy number changed the
mRNA of 4,304 genes (FDR < 5%), and the change was transmitted to protein for
793 of them (18%). Attenuation was larger for subunits of protein complexes
than for other genes (median 0.22 versus 0.17, p = 2.6 × 10⁻³¹; Figure 7B).

Driver alterations acted on different layers (Figure 7C). *TP53* missense
mutations raised p53 protein (log2 difference 0.64, FDR = 1.4 × 10⁻⁴) without
changing mRNA. *TP53* truncating mutations lowered mRNA (−0.28,
FDR = 0.004) without changing protein. Homozygous deletion of *SMAD4* lowered
both mRNA (−0.77, FDR = 2 × 10⁻⁴) and protein (−0.22, FDR = 0.005).
Homozygous deletion of *CDKN2A* lowered mRNA (−1.39, FDR = 2 × 10⁻⁴), while
missense and truncating *CDKN2A* mutations raised it (+1.79 and +1.41). *KRAS*
mutation raised *KRAS* mRNA (+0.39, FDR = 2 × 10⁻⁴) but not protein.

### 2.4 Tumour-versus-normal changes agree in direction but differ in magnitude

DESeq2 compared 160 tumours with 42 NATs and found 1,266 genes up and 1,259
down (FDR < 5%, |log2FC| > 1; Figure 4A). limma, with plex as covariate and
patient as a random block, compared 136 tumours with 66 NATs and found 1,229
proteins up and 1,269 down (FDR < 5%, |log2FC| > 0.5; Figure 4B). Fold changes
of the 9,213 genes in both layers correlated at rho = 0.40, and at 0.49 when
both models used the identical matched samples (Figure 4C).

Concordant increases (313 genes) were led by squamous and invasion genes
(*LGALS7*, *LY6D*, *KRT13*, *SPRR1A*, *A2ML1*). Concordant decreases (251)
were led by acinar enzymes (*CTRL*, *CEL*, *CTRB1*, *CELA3B*, *AMY2A*). Changes
in protein only (1,829) outnumbered changes in mRNA only (306). With the
identical matched samples these categories became 661 and 380. Only 16 genes
changed in opposite directions (Figure 4D).

### 2.5 Pathway activity is more concordant than individual genes

Pre-ranked GSEA (Korotkevich et al., 2016) of Hallmark sets (Liberzon et al.,
2015) gave the same direction in both layers for 17 of 50 sets (Figure 5A).
These included EMT (NES 2.2 in mRNA, 3.3 in protein), the interferon
responses, inflammatory response, TNF-alpha signalling, hypoxia and *KRAS*
signalling. Oxidative phosphorylation (NES 1.6 in mRNA, −3.1 in protein) and
MYC targets (2.5 and −2.1) changed in opposite directions. Agreement of
normalised enrichment scores was 0.30 for Hallmark sets but 0.04 for Reactome
sets (Gillespie et al., 2022). It rose to 0.12 when Reactome sets dominated by
ribosome, proteasome, spliceosome or respiratory-chain subunits were removed
(Figure 5B).

GSVA (Hänzelmann et al., 2013) scores of Hallmark pathways correlated between
layers with a median rho of 0.50 across tumours. The median gene-level
correlation of the same pathways was 0.43 (paired Wilcoxon p = 9.2 × 10⁻⁵;
Figure 5C). EMT reached 0.83 at pathway level against 0.60 for its genes.

### 2.6 The basal-like/classical axis is preserved in the proteome

PurIST called 17% of tumours basal-like. The Moffitt basal-minus-classical
score correlated between mRNA and protein at rho = 0.96 (Figure 6A). Both
scores identified the PurIST call with an AUC of 0.985 (mRNA) and 0.969
(protein). Consensus calls (Wilkerson and Hayes, 2010) agreed between layers
in 85% of tumours (Cohen's kappa 0.70; Figure 6B). The same genes separated the PurIST
groups in both layers (Figure 6D). They agreed less with PurIST (kappa
0.45 for mRNA, 0.38 for protein), because consensus clustering assigned
42-45% of tumours to the basal-like group.

Unsupervised consensus clustering of the 1,500 most variable features chose
two clusters for mRNA by the proportion of ambiguous clustering (PAC 0.07;
Şenbabaoğlu et al., 2014). For protein, PAC decreased up to k = 6 without a
clear optimum. Agreement between layers was absent at k = 2 (adjusted Rand
index −0.01) and modest at k = 3-6 (0.22-0.28; Figure 6C). These clusters
followed acinar content and cellularity rather than the tumour subtypes
(Figure 6E).

### 2.7 Subtypes and concordant genes carry the same prognostic information in both layers

Grade G3-G4 (HR 2.36, p = 2 × 10⁻⁴) and stage III-IV (HR 1.56, p = 0.04) were
the independent clinical predictors of overall survival. Basal-like tumours by
the Moffitt consensus call had shorter survival in both layers (log-rank
p = 2.7 × 10⁻⁴ for mRNA, 5.2 × 10⁻⁴ for protein; Figure 8A-B). Adjusted for
age, stage, margin and grade, the hazard ratio was 1.70 in both layers
(mRNA 95% CI 1.06-2.72; protein 1.07-2.70; Figure 8D). PurIST did not reach
significance (log-rank p = 0.065; adjusted HR 1.01). The per-tumour
mRNA-protein correlation was not associated with survival (HR 0.97 per SD,
p = 0.79).

In gene-level Cox models adjusted for stage, 770 genes in mRNA and 298 in
protein were associated with survival at p < 0.01. At FDR < 10% these numbers
were 650 and 15. Hazard ratios of the two layers correlated at rho = 0.44
across genes (Figure 8C). This agreement increased with the gene's
mRNA-protein correlation, from 0.07 for genes with rho < 0.2 to 0.39, 0.64
and 0.87 for genes with rho > 0.6 (Figure 8E).

### 2.8 Genome-wide models do not improve on clinical prediction

Elastic-net Cox models (Friedman et al., 2010; Simon et al., 2011) were
evaluated by Harrell's C on held-out folds (Harrell et al., 1996) in 50
repeats of 5-fold cross-validation, with the penalty tuned inside each
training set (Figure 8F). The clinical model reached a C-index of 0.622.
Adding 2,000 mRNA genes, 2,000 proteins, or 50 pathway scores of either layer
did not improve it (0.617-0.620). Models of molecular features alone reached
0.557-0.562. Adding only the pre-specified Moffitt mRNA score raised the
C-index to 0.632 and improved on the clinical model in 84% of repeats. The
protein version of the score improved it in 56%. Only two mRNA genes and three
proteins were selected in at least half of the training fits. Out-of-fold
risk groups of the protein model were not prognostic after adjustment (HR
1.14, p = 0.54; Figure 8G).

---
