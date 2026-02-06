# Step 05 - Quality control of both layers and tumour purity.
#
# RNA     : library size, detected genes, DESeq2 variance-stabilising transform
#           (VST) and PCA.
# Protein : missingness, distribution after median centring, PCA, and the share
#           of PC variance explained by TMT plex vs tissue.
# Purity  : ESTIMATE stromal/immune scores from RNA (tidyestimate) and a
#           genomic estimate of neoplastic cellularity, 2 x KRAS VAF.
#
#   Rscript scripts/R/05_qc_purity.R

source("R/utils.R")
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyestimate)
  library(patchwork)
  library(matrixStats)
})
set.seed(CFG$survival$seed)

samples <- read_samples()
genes   <- read_genes()

# ============================================================== RNA
log_step("RNA: VST")
counts <- read_matrix("rna_counts.tsv.gz")
pc_ids <- genes[gene_type == "protein_coding", gene_id]
counts_pc <- counts[rownames(counts) %in% pc_ids, ]
keep <- rowSums(counts_pc >= CFG$filters$rna_min_count) >= ceiling(0.1 * ncol(counts_pc))
coldata <- as.data.frame(samples[match(colnames(counts_pc), sample_key)])
rownames(coldata) <- coldata$sample_key
coldata$tissue <- factor(coldata$tissue, levels = c("Normal", "Tumor"))

dds <- DESeqDataSetFromMatrix(counts_pc[keep, ], coldata, design = ~tissue)
dds <- estimateSizeFactors(dds)
vsd <- assay(vst(dds, blind = TRUE))
vst_sym <- collapse_to_symbol(vsd, genes)

tpm <- read_matrix("rna_tpm.tsv.gz")
logtpm <- log2(collapse_to_symbol(tpm[rownames(tpm) %in% pc_ids, ], genes) + 1)

fwrite(data.table(gene = rownames(vst_sym), vst_sym), path_proc("rna_vst.tsv.gz"), sep = "\t")
fwrite(data.table(gene = rownames(logtpm), round(logtpm, 4)), path_proc("rna_log2tpm.tsv.gz"), sep = "\t")

qc_rna <- data.table(
  sample_key = colnames(counts_pc),
  tissue = coldata[colnames(counts_pc), "tissue"],
  assigned_reads_pc = colSums(counts_pc),
  detected_pc_genes = colSums(counts_pc > 0),
  size_factor = sizeFactors(dds)[colnames(counts_pc)]
)
write_tsv(qc_rna, "qc_rna_samples.tsv")

top <- order(rowVars(vst_sym), decreasing = TRUE)[1:2000]
pca_r <- prcomp(t(vst_sym[top, ]))
ve_r <- round(100 * pca_r$sdev^2 / sum(pca_r$sdev^2), 1)
pr <- data.table(sample_key = rownames(pca_r$x), PC1 = pca_r$x[, 1], PC2 = pca_r$x[, 2])
pr <- merge(pr, samples, by = "sample_key")

# ============================================================== tissue identity
# Whipple resections include the duodenum. Normal-adjacent samples dominated by
# small-intestinal mucosa are not pancreas and must not serve as the normal
# reference. Score = mean VST of enterocyte/Paneth/goblet markers vs acinar markers.
log_step("Tissue identity (duodenal mucosa check)")
INTESTINAL <- c("SI", "APOA4", "DEFA5", "DEFA6", "FABP2", "LCT", "MUC2", "CDX2", "VIL1",
                "FABP1", "APOB", "MTTP", "ALPI", "RBP2", "SLC26A3", "CDX1", "CDH17")
ACINAR <- c("PRSS1", "CPA1", "CELA3A", "CTRB1", "PNLIP", "CPB1", "CTRC", "CEL", "CELA2A", "PLA2G1B")
score <- function(m, g) colMeans(m[intersect(g, rownames(m)), , drop = FALSE], na.rm = TRUE)
ident <- data.table(sample_key = colnames(vst_sym),
                    intestinal_rna = score(vst_sym, INTESTINAL), acinar_rna = score(vst_sym, ACINAR))
ident <- merge(samples[, .(sample_key, tissue, has_rna, has_protein, included)], ident, by = "sample_key", all.x = TRUE)
cutoff <- 10   # provisional: gap seen in the score distribution
ident[, duodenal_flag := !is.na(intestinal_rna) & intestinal_rna > cutoff]
log_step(sprintf("intestinal score cut-off %.1f; flagged %s",
                 cutoff, paste(ident[duodenal_flag == TRUE, sample_key], collapse = ", ")))
write_tsv(ident[order(-intestinal_rna)], "qc_tissue_identity.tsv")

# analysis sample sheet used by every later step
an <- merge(samples, ident[, .(sample_key, intestinal_rna, acinar_rna, duodenal_flag)], by = "sample_key", all.x = TRUE)
an[is.na(duodenal_flag), duodenal_flag := FALSE]
an[, qc_pass := !duodenal_flag]
an[, matched := included == 1 & has_rna == 1 & has_protein == 1 & qc_pass]
fwrite(an, path_proc("analysis_samples.tsv"), sep = "\t")

# ============================================================== protein
log_step("Protein: QC")
prot <- read_matrix("protein_norm.tsv.gz")
pmeta <- samples[match(colnames(prot), sample_key)]
qc_prot <- data.table(
  sample_key = colnames(prot), tissue = pmeta$tissue, plex = pmeta$plex,
  n_quantified = colSums(!is.na(prot)),
  median_raw = NA_real_, iqr = colIQRs(prot, na.rm = TRUE)
)
raw <- read_matrix("protein_log2ratio.tsv.gz")
qc_prot$median_raw <- colMedians(raw[, colnames(prot)], na.rm = TRUE)
write_tsv(qc_prot, "qc_protein_samples.tsv")

# protein-level intestinal / acinar scores (same marker lists) for the identity table
pid <- data.table(sample_key = colnames(prot), intestinal_protein = score(prot, INTESTINAL),
                  acinar_protein = score(prot, ACINAR))
ident <- merge(ident, pid, by = "sample_key", all = TRUE)
ident[is.na(duodenal_flag), duodenal_flag := FALSE]
write_tsv(ident[order(-intestinal_rna)], "qc_tissue_identity.tsv")

frac_q <- rowMeans(!is.na(prot))
complete <- prot[frac_q == 1, ]
pca_p <- prcomp(t(complete[order(rowVars(complete), decreasing = TRUE)[1:2000], ]))
ve_p <- round(100 * pca_p$sdev^2 / sum(pca_p$sdev^2), 1)
pp <- data.table(sample_key = rownames(pca_p$x), PC1 = pca_p$x[, 1], PC2 = pca_p$x[, 2])
pp <- merge(pp, samples, by = "sample_key")

# variance of the first 10 PCs explained by tissue and by plex (adjusted R^2)
pc_r2 <- rbindlist(lapply(1:10, function(k) {
  d <- data.frame(y = pca_p$x[, k], tissue = pmeta$tissue, plex = pmeta$plex)
  data.table(PC = k, var_pct = ve_p[k],
             tissue_R2 = summary(lm(y ~ tissue, d))$adj.r.squared,
             plex_R2 = summary(lm(y ~ plex, d))$adj.r.squared,
             plex_given_tissue_R2 = summary(lm(y ~ tissue + plex, d))$adj.r.squared -
               summary(lm(y ~ tissue, d))$adj.r.squared)
}))
write_tsv(pc_r2, "qc_protein_pc_variance.tsv")

# ============================================================== purity
log_step("Purity: ESTIMATE + KRAS VAF")
est_in <- as.data.frame(vst_sym)   # tidy = FALSE: gene symbols are the row names
est <- est_in |>
  filter_common_genes(id = "hgnc_symbol", tidy = FALSE, tell_missing = FALSE, find_alias = TRUE) |>
  estimate_score(is_affymetrix = FALSE)
est <- as.data.table(est)
setnames(est, c("sample", "stromal", "immune", "estimate"), c("sample_key", "stromal_score", "immune_score", "estimate_score"))

mm <- fread(path_proc("mutation_matrix.tsv"))
pur <- merge(samples[, .(sample_key, case_submitter_id, tissue)], est, by = "sample_key", all.x = TRUE)
pur <- merge(pur, mm[, .(case_submitter_id, KRAS_allele, KRAS_vaf)], by = "case_submitter_id", all.x = TRUE)
pur[tissue == "Normal", c("KRAS_allele", "KRAS_vaf") := NA]
pur[, kras_cellularity := pmin(1, 2 * KRAS_vaf)]
write_tsv(pur[order(tissue, sample_key)], "purity_estimates.tsv")

ct <- pur[tissue == "Tumor" & !is.na(kras_cellularity) & !is.na(estimate_score) &
            !sample_key %in% ident[duodenal_flag == TRUE, sample_key]]
rho_pur <- cor.test(ct$estimate_score, ct$kras_cellularity, method = "spearman", exact = FALSE)
log_step(sprintf("ESTIMATE score vs KRAS cellularity (tumours, n=%d): rho=%.2f p=%.2g",
                 nrow(ct), rho_pur$estimate, rho_pur$p.value))

# ============================================================== figure 2
log_step("Figure 2")
pr <- merge(pr, ident[, .(sample_key, duodenal_flag)], by = "sample_key")
pr[, flag := factor(ifelse(duodenal_flag, "Duodenal mucosa", "Pancreas"), levels = c("Pancreas", "Duodenal mucosa"))]
pA <- ggplot(pr, aes(PC1, PC2, colour = tissue, shape = flag)) +
  geom_point(size = 1.2, alpha = 0.8) + scale_colour_manual(values = PAL$tissue) +
  scale_shape_manual(values = c(16, 4)) +
  labs(title = "RNA (VST, top 2,000 genes)", x = paste0("PC1 (", ve_r[1], "%)"),
       y = paste0("PC2 (", ve_r[2], "%)"), colour = NULL, shape = NULL) +
  theme_pub() + theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin())
idp <- ident[!is.na(intestinal_rna)]
pB <- ggplot(idp, aes(acinar_rna, intestinal_rna, colour = tissue, shape = duodenal_flag)) +
  geom_hline(yintercept = cutoff, linetype = 2, linewidth = 0.3) +
  geom_point(size = 1.1, alpha = 0.8) + scale_colour_manual(values = PAL$tissue, guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), guide = "none") +
  labs(title = "Tissue identity (RNA)", x = "Acinar marker score", y = "Intestinal marker score") + theme_pub()
pC <- ggplot(pp, aes(PC1, PC2, colour = tissue)) +
  geom_point(size = 1.2, alpha = 0.8) + scale_colour_manual(values = PAL$tissue, guide = "none") +
  labs(title = "Protein (top 2,000 complete)", x = paste0("PC1 (", ve_p[1], "%)"), y = paste0("PC2 (", ve_p[2], "%)")) +
  theme_pub()
r2l <- melt(pc_r2[PC <= 8], id.vars = c("PC", "var_pct"), measure.vars = c("tissue_R2", "plex_given_tissue_R2"))
r2l[, variable := factor(variable, labels = c("Tissue", "TMT plex | tissue"))]
pD <- ggplot(r2l, aes(factor(PC), pmax(value, 0), fill = variable)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("#B2182B", "grey55")) +
  labs(title = "Protein PC drivers", x = "Principal component", y = expression(Adjusted~R^2), fill = NULL) +
  theme_pub() + theme(legend.position = c(0.7, 0.85))
miss <- data.table(frac = frac_q)
pE <- ggplot(miss, aes(frac)) + geom_histogram(bins = 40, fill = PAL$layer[["Protein"]]) +
  geom_vline(xintercept = CFG$filters$protein_min_fraction_quantified, linetype = 2) +
  labs(title = "Protein quantification rate", x = "Fraction of samples quantified", y = "Proteins") + theme_pub()
pF <- ggplot(ct, aes(kras_cellularity, estimate_score)) +
  geom_point(size = 1, colour = "grey30") + geom_smooth(method = "lm", formula = y ~ x, colour = "#B2182B", linewidth = 0.6) +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.4, size = 2.8,
           label = sprintf("Spearman rho = %.2f\nn = %d", rho_pur$estimate, nrow(ct))) +
  labs(title = "Two purity proxies", x = "Neoplastic cellularity (2 x KRAS VAF)", y = "ESTIMATE score") + theme_pub()

fig2 <- (pA | pB | pC) / (pD | pE | pF) + plot_annotation(tag_levels = "A")
save_fig(fig2, "Fig2_qc_purity", width = 210, height = 145)

log_step(sprintf("RNA: %d protein-coding genes kept for VST, %d samples; protein: %d proteins, %d complete",
                 sum(keep), ncol(vsd), nrow(prot), nrow(complete)))
print(pc_r2)
