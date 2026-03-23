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

## 3. Discussion

PDAC is still classified, staged for trials and searched for biomarkers almost
entirely at the mRNA level. Resected PDAC is a low-cellularity tumour whose
bulk profile mixes cancer cells, a dense desmoplastic stroma, immune cells and
residual acinar parenchyma (Moffitt et al., 2015). Every bulk measurement is
therefore an average over cell types whose mRNA and protein are regulated
differently. Proteomics of the same tissue is the direct way to ask which
transcriptomic conclusions survive translation, and CPTAC-3 has made such data
public for PDAC (Cao et al., 2021).

In this study, 136 tumours and 21 NATs with both layers from the same tissue
powder were analysed after explicit sample matching and the removal of 11
samples of duodenal origin. The median gene-wise mRNA-protein correlation was
0.30 in tumours and 0.23 in NATs, and 68% of genes were significantly coupled.
Copy-number effects reached the protein for 18% of genes. The basal-like axis
was preserved in protein (rho 0.96) and was equally prognostic (HR 1.70).
Genome-wide models did not exceed a clinical C-index of 0.622, whereas a
single subtype score reached 0.632.

The level of mRNA-protein coupling found here agrees with earlier CPTAC
studies. In a similar analysis of the same cohort, the median correlation
was 0.36 in tumours and 0.26 in NATs (Cao et al., 2021). The lower value here
is probably due to processing: we used the public common-pipeline ratios
(Rudnick et al., 2016) rather than the consortium's own quantification, and
the 136-tumour set includes low-purity tumours. Across tumour types, the
correlation varies widely and is lowest for housekeeping machinery (Arad and
Geiger, 2023). In cultured cells, translation and degradation explain as much
protein variation as mRNA does (Schwanhäusser et al., 2011). In accordance with
all these studies and the study under discussion we can infer that mRNA is a
useful but partial reporter of protein, and that its reliability is a
property of the gene class. Composition accounted for only a small part of the
coupling (median 0.30 to 0.26), so most of the agreement for stromal and
interferon genes reflects transcriptional programmes that are translated
faithfully, not only shifts in cell proportions.

The genes least coupled were the subunits of large complexes, and the same
genes buffered copy-number changes. Gonçalves et al. (2017) reported
post-transcriptional attenuation of copy-number variation across cancer types,
strongest for complex subunits. Liu et al. (2016) attributed this to the
degradation of subunits produced in excess of their partners. Here,
attenuation was larger for complex members (0.22 versus 0.17), and only 18% of
mRNA dosage effects reached the protein. It can be inferred that aneuploidy in
PDAC is largely buffered at the protein level, so copy-number-based target
nomination must be confirmed in protein. The driver analysis adds a direct
illustration. Missense *TP53* raised protein but not mRNA, which is consistent
with stabilisation of the mutant protein, while truncating *TP53* lowered mRNA,
which is consistent with nonsense-mediated decay. Neither effect would be
seen from one layer alone.

Tumour-versus-normal comparisons agreed in direction (rho 0.40-0.49), and
pathway activity agreed better than single genes (0.50 versus 0.43). The
opposite signs for oxidative phosphorylation and MYC targets are probably
due to tissue composition. Normal pancreas is dominated by acinar cells with
abundant mitochondria and ribosomes, so the tumour has less of this protein
mass, while proliferating cancer cells transcribe more of these genes. The
low Reactome agreement (0.04) is mostly explained by the redundancy of
machinery-dominated sets, because it tripled when those sets were removed.
For this reason, curated non-redundant collections such as Hallmark
(Liberzon et al., 2015) should be preferred for cross-layer comparison.

The tumour subtypes were the most robust cross-layer signal. Moffitt et al.
(2015) derived the basal-like and classical programmes from transcriptomes
after virtual microdissection, and Bailey et al. (2016) found a corresponding
squamous class. In this cohort the protein score correlated with the mRNA
score at 0.96 and carried the same hazard ratio. Basal-like tumours also
respond less to modified FOLFIRINOX (Aung et al., 2018). In accordance with
these studies and the study under discussion we can infer that the
basal-like state is a biological state of the tumour expressed through
abundant structural and secreted proteins, not a transcriptomic artefact.
PurIST was not significant here although it is prognostic elsewhere (Rashid
et al., 2020), probably because it calls only the clearly basal minority
(17%) and the confidence interval in about 23 patients was wide.

The prognostic analyses connect concordance and outcome. Genes with
concordant mRNA and protein carried the same prognostic information in both
layers (agreement 0.87 for genes with rho > 0.6), and discordant genes carried
unrelated information. Ghoshdastider and Sendoel (2023) linked weaker
post-transcriptional control to shorter survival across cancer types, but the
per-tumour correlation was not prognostic here (HR 0.97). This is probably
because it mainly follows cellularity in PDAC. Cross-validation showed that
neither layer improved prediction when used genome-wide, and that a single
pre-specified score did. Harrell et al. (1996) warned that models with many
candidate predictors and few events overfit. The instability of the selected
genes here supports that warning. It can be inferred that, in cohorts of this
size, the value of the proteome is mechanistic, and prognostic modelling
should use pre-specified scores.

Spearman correlation was used for concordance because TMT ratios are
compressed and a few extreme tumours can dominate a linear correlation. A
permutation null was preferred to the analytical one because shared
missingness and correlated genes inflate nominal significance. Cross-validated
C-indices with nested tuning were reported instead of resubstitution
estimates because only performance on unseen patients is informative (Harrell
et al., 1996). The concordance values found here are consistent with the
original analysis of the cohort (Cao et al., 2021), which supports the use of
the public common-pipeline data for re-analysis.

**Limitations.** This is a single cohort without external validation. Most
patients had stage II disease and survival is overall rather than
disease-specific. Normal-adjacent tissue is not normal pancreas, and 11
samples had to be removed for duodenal contamination; histology was not
available to confirm tissue identity. TMT ratio compression reduces protein
effect sizes and may explain part of the weaker protein prognostic signal.
Purity was estimated, not measured, and both estimates are noisy. Phospho- and
glycoproteomes, which the original study analysed, were not used.

Based on this study, it can be inferred that mRNA reliably reports protein
for stromal, immune and subtype-defining genes but not for complex machinery
or copy-number-driven changes, and that protein adds mechanistic rather than
prognostic information in resected PDAC; therefore, transcriptomic biomarkers
should be checked for protein concordance before clinical use, normal
samples from Whipple specimens should be screened for duodenal contamination,
and prognostic models should rest on pre-specified biological scores
validated in independent cohorts.

---

## 4. Methods (summary)

Full methods are in [`docs/methods.md`](../docs/methods.md).

**Data.** Open-access CPTAC-3 pancreas files were downloaded from the GDC
(Grossman et al., 2016): STAR gene counts (Dobin et al., 2013), masked somatic
mutations, ASCAT gene-level copy number (Van Loo et al., 2010) and clinical
records. The gene-level log2 ratios of PDC study PDC000270 were downloaded
from the PDC (Rudnick et al., 2016). All files were md5-verified.

**Matching and QC.** Samples were matched through shared sample identifiers
and recorded in an SQLite database with match level and exclusion reason.
Duodenal samples were flagged with a 17-gene intestinal marker score.
Cellularity was estimated from *KRAS* allele fraction and ESTIMATE (Yoshihara
et al., 2013).

**Statistics.** Gene-wise Spearman correlations were corrected by the
Benjamini-Hochberg method (Benjamini and Hochberg, 1995). RNA differential
expression used DESeq2 (Love et al., 2014) with apeglm shrinkage (Zhu et al.,
2019). Protein differential expression used limma (Ritchie et al., 2015) with
duplicateCorrelation (Smyth et al., 2005). Pathways used fgsea (Korotkevich et
al., 2016) and GSVA (Hänzelmann et al., 2013). Subtypes used PurIST (Rashid et
al., 2020), the Moffitt signatures (Moffitt et al., 2015) and
ConsensusClusterPlus (Wilkerson and Hayes, 2010) with PAC (Şenbabaoğlu et al.,
2014), Cohen's kappa (Cohen, 1960) and the adjusted Rand index (Hubert and
Arabie, 1985). Survival analysis used Cox models with Schoenfeld residual
checks (Grambsch and Therneau, 1994), elastic-net Cox (Friedman et al., 2010;
Simon et al., 2011) and Harrell's C (Harrell et al., 1996).

**Software.** R 4.5.2, Bioconductor 3.22, CRAN snapshot of 15 December 2025
(`renv.lock`) and Python 3.12 (`envs/environment.yml`). The full pipeline is
run with `python run_pipeline.py`.

## Data and code availability

All data are open access at the GDC (project CPTAC-3) and the PDC (study
PDC000270). Code, result tables, figures and the database schema are in this
repository.

## Figures

1. Study design, data availability, sample matching and overall survival.
2. Quality control: RNA and protein PCA, tissue identity, plex effects, purity.
3. mRNA-protein concordance at gene, gene-set and tumour level.
4. Tumour versus normal differential expression in mRNA and protein.
5. Pathway enrichment and pathway-level concordance.
6. PDAC subtypes in mRNA and protein; unsupervised clustering.
7. Copy-number dosage attenuation and driver cis effects.
8. Survival: subtypes, gene-level prognosis and cross-validated models.

S1. Survival by subtype classifier and unsupervised cluster.

## References

* Arad G, Geiger T (2023). Functional Impact of Protein–RNA Variation in Clinical Cancer Analyses. *Molecular & Cellular Proteomics* 22(7):100587. https://doi.org/10.1016/j.mcpro.2023.100587
* Aung KL, Fischer SE, Denroche RE, Jang GH, Dodd A, Creighton S, et al. (2018). Genomics-Driven Precision Medicine for Advanced Pancreatic Cancer: Early Results from the COMPASS Trial. *Clinical Cancer Research* 24(6):1344-1354. https://doi.org/10.1158/1078-0432.ccr-17-2994
* Bailey P, Chang DK, Nones K, Johns AL, Patch AM, Gingras MC, et al. (2016). Genomic analyses identify molecular subtypes of pancreatic cancer. *Nature* 531(7592):47-52. https://doi.org/10.1038/nature16965
* Benjamini Y, Hochberg Y (1995). Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. *Journal of the Royal Statistical Society Series B: Statistical Methodology* 57(1):289-300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
* Cao L, Huang C, Cui Zhou D, Hu Y, Lih TM, Savage SR, et al. (2021). Proteogenomic characterization of pancreatic ductal adenocarcinoma. *Cell* 184(19):5031-5052.e26. https://doi.org/10.1016/j.cell.2021.08.023
* Cohen J (1960). A Coefficient of Agreement for Nominal Scales. *Educational and Psychological Measurement* 20(1):37-46. https://doi.org/10.1177/001316446002000104
* Collisson EA, Sadanandam A, Olson P, Gibb WJ, Truitt M, Gu S, et al. (2011). Subtypes of pancreatic ductal adenocarcinoma and their differing responses to therapy. *Nature Medicine* 17(4):500-503. https://doi.org/10.1038/nm.2344
* Conroy T, Hammel P, Hebbar M, Ben Abdelghani M, Wei AC, Raoul JL, et al. (2018). FOLFIRINOX or Gemcitabine as Adjuvant Therapy for Pancreatic Cancer. *New England Journal of Medicine* 379(25):2395-2406. https://doi.org/10.1056/nejmoa1809775
* Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, et al. (2013). STAR: ultrafast universal RNA-seq aligner. *Bioinformatics* 29(1):15-21. https://doi.org/10.1093/bioinformatics/bts635
* Friedman J, Hastie T, Tibshirani R (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. *Journal of Statistical Software* 33(1). https://doi.org/10.18637/jss.v033.i01
* Ghoshdastider U, Sendoel A (2023). Exploring the pan-cancer landscape of posttranscriptional regulation. *Cell Reports* 42(10):113172. https://doi.org/10.1016/j.celrep.2023.113172
* Gillespie M, Jassal B, Stephan R, Milacic M, Rothfels K, Senff-Ribeiro A, et al. (2022). The reactome pathway knowledgebase 2022. *Nucleic Acids Research* 50(D1):D687-D692. https://doi.org/10.1093/nar/gkab1028
* Gonçalves E, Fragoulis A, Garcia-Alonso L, Cramer T, Saez-Rodriguez J, Beltrao P (2017). Widespread Post-transcriptional Attenuation of Genomic Copy-Number Variation in Cancer. *Cell Systems* 5(4):386-398.e4. https://doi.org/10.1016/j.cels.2017.08.013
* Grambsch PM, Therneau TM (1994). Proportional hazards tests and diagnostics based on weighted residuals. *Biometrika* 81(3):515-526. https://doi.org/10.1093/biomet/81.3.515
* Grossman RL, Heath AP, Ferretti V, Varmus HE, Lowy DR, Kibbe WA, et al. (2016). Toward a Shared Vision for Cancer Genomic Data. *New England Journal of Medicine* 375(12):1109-1112. https://doi.org/10.1056/nejmp1607591
* Hänzelmann S, Castelo R, Guinney J (2013). GSVA: gene set variation analysis for microarray and RNA-Seq data. *BMC Bioinformatics* 14(1). https://doi.org/10.1186/1471-2105-14-7
* Harrell FE, Lee KL, Mark DB (1996). Multivariable prognostic models: issues in developing models, evaluating assumptions and adequacy, and measuring and reducing errors. *Statistics in Medicine* 15(4):361-387. https://doi.org/10.1002/(sici)1097-0258(19960229)15:4<361::aid-sim168>3.0.co;2-4
* Hubert L, Arabie P (1985). Comparing partitions. *Journal of Classification* 2(1):193-218. https://doi.org/10.1007/bf01908075
* Korotkevich G, Sukhov V, Budin N, Shpak B, Artyomov MN, Sergushichev A (2016). Fast gene set enrichment analysis. *bioRxiv*. https://doi.org/10.1101/060012
* Li Y, Dou Y, Da Veiga Leprevost F, Geffen Y, Calinawan AP, Aguet F, et al. (2023). Proteogenomic data and resources for pan-cancer analysis. *Cancer Cell* 41(8):1397-1406. https://doi.org/10.1016/j.ccell.2023.06.009
* Liberzon A, Birger C, Thorvaldsdóttir H, Ghandi M, Mesirov JP, Tamayo P (2015). The Molecular Signatures Database Hallmark Gene Set Collection. *Cell Systems* 1(6):417-425. https://doi.org/10.1016/j.cels.2015.12.004
* Liu Y, Beyer A, Aebersold R (2016). On the Dependency of Cellular Protein Levels on mRNA Abundance. *Cell* 165(3):535-550. https://doi.org/10.1016/j.cell.2016.03.014
* Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15(12). https://doi.org/10.1186/s13059-014-0550-8
* Moffitt RA, Marayati R, Flate EL, Volmar KE, Loeza SGH, Hoadley KA, et al. (2015). Virtual microdissection identifies distinct tumor- and stroma-specific subtypes of pancreatic ductal adenocarcinoma. *Nature Genetics* 47(10):1168-1178. https://doi.org/10.1038/ng.3398
* Rahib L, Smith BD, Aizenberg R, Rosenzweig AB, Fleshman JM, Matrisian LM (2014). Projecting Cancer Incidence and Deaths to 2030: The Unexpected Burden of Thyroid, Liver, and Pancreas Cancers in the United States. *Cancer Research* 74(11):2913-2921. https://doi.org/10.1158/0008-5472.can-14-0155
* Rashid NU, Peng XL, Jin C, Moffitt RA, Volmar KE, Belt BA, et al. (2020). Purity Independent Subtyping of Tumors (PurIST), A Clinically Robust, Single-sample Classifier for Tumor Subtyping in Pancreatic Cancer. *Clinical Cancer Research* 26(1):82-92. https://doi.org/10.1158/1078-0432.ccr-19-1467
* Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, et al. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research* 43(7):e47-e47. https://doi.org/10.1093/nar/gkv007
* Rudnick PA, Markey SP, Roth J, Mirokhin Y, Yan X, Tchekhovskoi DV, et al. (2016). A Description of the Clinical Proteomic Tumor Analysis Consortium (CPTAC) Common Data Analysis Pipeline. *Journal of Proteome Research* 15(3):1023-1032. https://doi.org/10.1021/acs.jproteome.5b01091
* Savage SR, Yi X, Lei JT, Wen B, Zhao H, Liao Y, et al. (2024). Pan-cancer proteogenomics expands the landscape of therapeutic targets. *Cell* 187(16):4389-4407.e15. https://doi.org/10.1016/j.cell.2024.05.039
* Schwanhäusser B, Busse D, Li N, Dittmar G, Schuchhardt J, Wolf J, et al. (2011). Global quantification of mammalian gene expression control. *Nature* 473(7347):337-342. https://doi.org/10.1038/nature10098
* Șenbabaoğlu Y, Michailidis G, Li JZ (2014). Critical limitations of consensus clustering in class discovery. *Scientific Reports* 4(1). https://doi.org/10.1038/srep06207
* Siegel RL, Giaquinto AN, Jemal A (2024). Cancer statistics, 2024. *CA: A Cancer Journal for Clinicians* 74(1):12-49. https://doi.org/10.3322/caac.21820
* Simon N, Friedman J, Hastie T, Tibshirani R (2011). Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent. *Journal of Statistical Software* 39(5). https://doi.org/10.18637/jss.v039.i05
* Smyth GK, Michaud J, Scott HS (2005). Use of within-array replicate spots for assessing differential expression in microarray experiments. *Bioinformatics* 21(9):2067-2075. https://doi.org/10.1093/bioinformatics/bti270
* Van Loo P, Nordgard SH, Lingjærde OC, Russnes HG, Rye IH, Sun W, et al. (2010). Allele-specific copy number analysis of tumors. *Proceedings of the National Academy of Sciences* 107(39):16910-16915. https://doi.org/10.1073/pnas.1009843107
* Wilkerson MD, Hayes DN (2010). ConsensusClusterPlus: a class discovery tool with confidence assessments and item tracking. *Bioinformatics* 26(12):1572-1573. https://doi.org/10.1093/bioinformatics/btq170
* Yoshihara K, Shahmoradgoli M, Martínez E, Vegesna R, Kim H, Torres-Garcia W, et al. (2013). Inferring tumour purity and stromal and immune cell admixture from expression data. *Nature Communications* 4(1). https://doi.org/10.1038/ncomms3612
* Zhu A, Ibrahim JG, Love MI (2019). Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. *Bioinformatics* 35(12):2084-2092. https://doi.org/10.1093/bioinformatics/bty895
