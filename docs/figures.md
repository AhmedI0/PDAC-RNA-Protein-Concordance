# Figure guide

Every figure is in `results/figures/` as SVG (editable vector), PNG and JPG
(300 dpi). The script that draws it is given in brackets.

## Figure 1. Cohort and design (`06_cohort_overview.py`)

* **A** Data flow: GDC files (RNA-seq, WXS, WGS copy number, clinical) and the
  PDC proteome are linked in an SQLite database at the case > sample >
  aliquot level, followed by QC and the matched analysis cohort.
* **B** Number of PDAC cases with each data type after QC.
* **C** Matched RNA-protein samples by match level: dark = proteome aliquot
  from the same tissue piece as the RNA; light = proteome aliquot is a pool
  of tumour pieces that contains the RNA piece.
* **D** Kaplan-Meier overall survival of the matched patients, with 95%
  confidence band and number at risk.

## Figure 2. Quality control and purity (`05_qc_purity.R`)

* **A** RNA PCA (VST, 2,000 most variable genes); crosses = samples flagged
  as duodenal mucosa (the isolated group on PC2).
* **B** Intestinal versus acinar marker scores; dashed line = cut-off half-way
  between the pancreas and duodenal centres.
* **C** Protein PCA on proteins quantified in every sample; PC1 separates
  tumour from normal.
* **D** Adjusted R² of each protein PC explained by tissue and by TMT plex
  (given tissue).
* **E** Fraction of samples in which each protein is quantified; dashed line =
  50% filter.
* **F** ESTIMATE score against KRAS-VAF cellularity in tumours.

## Figure 3. mRNA-protein concordance (`07_rna_protein_concordance.py`)

* **A** Gene-wise Spearman rho across 136 tumours (red) and the permutation
  null (grey).
* **B** Density of gene-wise rho in tumours, normals, tumours + normals, and
  tumours after removing composition (stromal, immune, acinar scores).
* **C** rho by decile of mRNA standard deviation across tumours.
* **D** Hallmark sets with the highest and lowest median rho (red = FDR < 5%
  versus all other genes; dashed line = overall median).
* **E** Protein complexes and assemblies (GO cellular component).
* **F** Per-tumour (sample-wise) agreement against KRAS-VAF cellularity.

## Figure 4. Tumour versus normal (`09_differential_expression.R`)

* **A** mRNA volcano plot (DESeq2, apeglm-shrunk log2 fold change).
* **B** Protein volcano plot (limma with plex and patient block).
* **C** mRNA versus protein log2 fold change, coloured by category, with
  labelled examples of each.
* **D** Number of genes per category.

## Figure 5. Pathways (`10_pathways.R`)

* **A** Hallmark GSEA normalised enrichment scores, mRNA versus protein,
  coloured by where the set is significant (FDR < 5%).
* **B** The 24 Reactome pathways with the largest enrichment in either layer;
  asterisk = FDR < 5%.
* **C** Pathway-level agreement (Spearman rho of mRNA and protein GSVA scores
  across tumours) against the median gene-level rho of the pathway's members;
  points above the diagonal agree better as pathways than as genes.
* **D** Examples of high, median and low pathway-level agreement.

## Figure 6. Subtypes (`11_subtypes.R`)

* **A** Moffitt basal-minus-classical score in mRNA versus protein, coloured
  by the PurIST call.
* **B** Moffitt consensus calls, mRNA versus protein (Cohen's kappa).
* **C** Unsupervised consensus clustering: PAC by k for each layer (solid) and
  adjusted Rand index between layers (dashed).
* **D** Heatmap of the Moffitt basal-like and classical genes in both layers,
  tumours ordered by PurIST probability; top row = PurIST call.
* **E** Cellularity, stromal score, acinar score and Moffitt protein score of
  the six protein clusters.

## Figure 7. Genomic cis effects (`08_genomic_cis_effects.py`)

* **A** Copy-number-mRNA versus copy-number-protein correlation per gene
  (hexbin density); dashed line = identity.
* **B** Attenuation (CN-mRNA rho minus CN-protein rho) for complex subunits
  and other genes, among genes with an mRNA dosage effect.
* **C** mRNA (green) and protein (purple) of TP53, SMAD4, CDKN2A and KRAS by
  alteration class: WT, missense/in-frame (mis), truncating (trunc),
  homozygous deletion (HD); numbers = tumours.

## Figure 8. Survival and prognostic models (`12_survival.R`, `13_prognostic_models.R`)

* **A-B** Kaplan-Meier curves by Moffitt consensus subtype in mRNA and in
  protein, with risk tables.
* **C** Stage-adjusted log hazard ratio per SD for every gene, mRNA versus
  protein; red = p < 0.01 in both layers.
* **D** Hazard ratios from Cox models adjusted for age, stage, margin and
  grade.
* **E** Agreement of gene-level log hazard ratios between layers, by the
  gene's mRNA-protein correlation.
* **F** Held-out Harrell's C of each model over 50 repeats of 5-fold
  cross-validation; dotted line = 0.5, dashed line = clinical model.
* **G** Survival of out-of-fold risk groups from the protein elastic-net
  model.

## Figure S1. Survival by subtype classifier (`12_survival.R`)

* **A-C** PurIST and Moffitt consensus calls (mRNA, protein).
* **D-E** Unsupervised clusters (mRNA k = 2, protein k = 6).
* **F** Same as Figure 8E.
