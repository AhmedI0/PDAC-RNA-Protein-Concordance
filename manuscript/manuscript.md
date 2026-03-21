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
