# Step 09 - Tumour vs normal-adjacent tissue at the mRNA and protein level.
#
# RNA     : DESeq2 negative-binomial GLM on raw STAR counts, design ~ tissue,
#           log2 fold changes shrunk with apeglm (Zhu et al. 2019).
# Protein : limma on median-centred TMT log2 ratios, design ~ tissue + plex,
#           patient as a random block (duplicateCorrelation) because some
#           tumours and normals come from the same patient; empirical Bayes
#           moderation with intensity trend.
# Primary analysis uses every QC-pass PDAC sample of each layer; a sensitivity
# analysis restricts both layers to the identical matched samples. A second
# sensitivity analysis asks whether ignoring tumour-NAT pairing in the RNA
# model matters: limma-voom with and without a patient block, and DESeq2 on the
# pairs alone with and without patient in the design.
#
#   Rscript scripts/R/09_differential_expression.R

source("R/utils.R")
suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
  library(fgsea)
  library(ggrepel)
  library(patchwork)
  library(BiocParallel)
})
register(SnowParam(workers = 6))
set.seed(CFG$survival$seed)
FDR <- CFG$filters$fdr
LFC_RNA <- CFG$filters$lfc_threshold      # |log2FC| > 1 for mRNA
LFC_PROT <- 0.5                            # TMT ratio compression: protein fold changes are smaller

an <- fread(path_proc("analysis_samples.tsv"))
genes <- read_genes()
an <- an[qc_pass == TRUE]

# ============================================================== RNA (DESeq2)
# Protein-coding genes with >= rna_min_count reads in at least
# rna_min_samples_fraction of the smaller group; shared by DESeq2 and voom.
rna_filtered <- function(keys) {
  counts <- read_matrix("rna_counts.tsv.gz")[, keys]
  pc <- genes[gene_type == "protein_coding", gene_id]
  counts <- counts[rownames(counts) %in% pc, ]
  cd <- as.data.frame(an[match(keys, sample_key)])
  rownames(cd) <- cd$sample_key
  cd$tissue <- factor(cd$tissue, levels = c("Normal", "Tumor"))
  cd$patient <- factor(make.names(cd$case_submitter_id))
  n_min <- ceiling(CFG$filters$rna_min_samples_fraction * min(table(cd$tissue)))
  keep <- rowSums(counts >= CFG$filters$rna_min_count) >= n_min
  list(counts = counts[keep, ], cd = cd)
}

run_deseq <- function(keys, label, design = ~tissue) {
  x <- rna_filtered(keys)
  cd <- x$cd
  dds <- DESeqDataSetFromMatrix(x$counts, cd, design = design)
  dds <- DESeq(dds, parallel = TRUE, quiet = TRUE)
  res <- results(dds, name = "tissue_Tumor_vs_Normal", alpha = FDR)
  shr <- lfcShrink(dds, coef = "tissue_Tumor_vs_Normal", type = "apeglm", res = res, quiet = TRUE)
  out <- data.table(gene_id = rownames(res), baseMean = res$baseMean, log2FC_mle = res$log2FoldChange,
                    log2FC = shr$log2FoldChange, lfcSE = shr$lfcSE, stat = res$stat,
                    pvalue = res$pvalue, padj = res$padj)
  out[, gene := genes$gene_name[match(gene_id, genes$gene_id)]]
  setorder(out, padj, -baseMean)
  out <- out[!duplicated(gene)]
  log_step(sprintf("DESeq2 %s: %d tumours vs %d normals, %d genes, %d up / %d down (FDR<%.2f, |LFC|>%g)",
                   label, sum(cd$tissue == "Tumor"), sum(cd$tissue == "Normal"), nrow(out),
                   out[padj < FDR & log2FC > LFC_RNA, .N], out[padj < FDR & log2FC < -LFC_RNA, .N], FDR, LFC_RNA))
  out[]
}

# ============================================================== protein (limma)
run_limma <- function(keys, label) {
  prot <- read_matrix("protein_norm.tsv.gz")[, keys]
  meta <- an[match(keys, sample_key)]
  tissue <- factor(meta$tissue, levels = c("Normal", "Tumor"))
  ok <- rowMeans(!is.na(prot[, tissue == "Tumor"])) >= CFG$filters$protein_min_fraction_quantified &
        rowMeans(!is.na(prot[, tissue == "Normal"])) >= CFG$filters$protein_min_fraction_quantified
  prot <- prot[ok, ]
  plex <- factor(make.names(meta$plex))
  design <- model.matrix(~ tissue + plex)
  block <- meta$case_submitter_id
  dc <- duplicateCorrelation(prot, design, block = block)
  fit <- lmFit(prot, design, block = block, correlation = dc$consensus.correlation)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  tt <- as.data.table(topTable(fit, coef = "tissueTumor", number = Inf, sort.by = "none"), keep.rownames = "gene")
  out <- tt[, .(gene, log2FC = logFC, AveExpr, t, pvalue = P.Value, padj = adj.P.Val,
                n_tumour = rowSums(!is.na(prot[gene, tissue == "Tumor", drop = FALSE])),
                n_normal = rowSums(!is.na(prot[gene, tissue == "Normal", drop = FALSE])))]
  log_step(sprintf("limma %s: %d tumours vs %d normals, %d proteins, block corr %.2f, %d up / %d down (FDR<%.2f, |LFC|>%g)",
                   label, sum(tissue == "Tumor"), sum(tissue == "Normal"), nrow(out), dc$consensus.correlation,
                   out[padj < FDR & log2FC > LFC_PROT, .N], out[padj < FDR & log2FC < -LFC_PROT, .N], FDR, LFC_PROT))
  out[]
}

classify <- function(m) {
  r_up <- m$padj_rna < FDR & m$log2FC_rna > LFC_RNA
  r_dn <- m$padj_rna < FDR & m$log2FC_rna < -LFC_RNA
  p_up <- m$padj_protein < FDR & m$log2FC_protein > LFC_PROT
  p_dn <- m$padj_protein < FDR & m$log2FC_protein < -LFC_PROT
  fifelse(r_up & p_up, "Concordant up",
  fifelse(r_dn & p_dn, "Concordant down",
  fifelse((r_up & p_dn) | (r_dn & p_up), "Opposite",
  fifelse((r_up | r_dn) & !(p_up | p_dn), "mRNA only",
  fifelse((p_up | p_dn) & !(r_up | r_dn), "Protein only", "Neither")))))
}

merge_layers <- function(r, p) {
  m <- merge(r[, .(gene, baseMean, log2FC_rna = log2FC, stat_rna = stat, padj_rna = padj)],
             p[, .(gene, AveExpr, log2FC_protein = log2FC, t_protein = t, padj_protein = padj)], by = "gene")
  m[, category := classify(m)]
  m[]
}

# ---------------------------------------------------------------- primary
rna_keys  <- an[has_rna == 1, sample_key]
prot_keys <- an[has_protein == 1, sample_key]
de_rna  <- run_deseq(rna_keys, "all QC-pass")
de_prot <- run_limma(prot_keys, "all QC-pass")
both <- merge_layers(de_rna, de_prot)
write_tsv(de_rna, "de_rna_tumour_vs_normal.tsv")
write_tsv(de_prot, "de_protein_tumour_vs_normal.tsv")
write_tsv(both[order(category, -abs(log2FC_rna))], "de_rna_protein_categories.tsv")

# ---------------------------------------------------------------- sensitivity: identical matched samples
mkeys <- an[matched == TRUE, sample_key]
de_rna_m  <- run_deseq(mkeys, "matched")
de_prot_m <- run_limma(mkeys, "matched")
both_m <- merge_layers(de_rna_m, de_prot_m)
write_tsv(both_m[order(category, -abs(log2FC_rna))], "de_rna_protein_categories_matched.tsv")

fc_cor  <- cor.test(both$log2FC_rna, both$log2FC_protein, method = "spearman", exact = FALSE)
fc_corm <- cor.test(both_m$log2FC_rna, both_m$log2FC_protein, method = "spearman", exact = FALSE)
cat_tab <- merge(both[, .(primary = .N), by = category], both_m[, .(matched_only = .N), by = category], by = "category", all = TRUE)
summ <- data.table(
  metric = c("RNA genes tested", "RNA up", "RNA down", "proteins tested", "protein up", "protein down",
             "genes in both layers", "Spearman rho of fold changes", "Spearman rho of fold changes (matched samples)"),
  value = c(nrow(de_rna), de_rna[padj < FDR & log2FC > LFC_RNA, .N], de_rna[padj < FDR & log2FC < -LFC_RNA, .N],
            nrow(de_prot), de_prot[padj < FDR & log2FC > LFC_PROT, .N], de_prot[padj < FDR & log2FC < -LFC_PROT, .N],
            nrow(both), round(fc_cor$estimate, 3), round(fc_corm$estimate, 3)))
write_tsv(summ, "de_summary.tsv")
write_tsv(cat_tab, "de_category_counts.tsv")
print(summ); print(cat_tab)

# ============================================================== figure 4
volcano <- function(d, lfc_thr, title, lab_n = 12) {
  d <- copy(d)[!is.na(padj)]
  d[, dir := fifelse(padj < FDR & log2FC > lfc_thr, "Up", fifelse(padj < FDR & log2FC < -lfc_thr, "Down", "NS"))]
  d[, nlp := pmin(-log10(padj), 60)]
  lab <- rbind(d[dir == "Up"][order(padj)][1:lab_n], d[dir == "Down"][order(padj)][1:lab_n])
  ggplot(d, aes(log2FC, nlp, colour = dir)) +
    geom_point(size = 0.35, alpha = 0.6) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = 3, linewidth = 0.3) +
    geom_hline(yintercept = -log10(FDR), linetype = 3, linewidth = 0.3) +
    geom_text_repel(data = lab, aes(label = gene), size = 2, max.overlaps = 30, segment.size = 0.2,
                    colour = "black", min.segment.length = 0) +
    scale_colour_manual(values = PAL$dir, guide = "none") +
    labs(title = title,
         subtitle = sprintf("%s up, %s down (FDR < %g, |log2FC| > %g); -log10 FDR capped at 60",
                            format(d[dir == "Up", .N], big.mark = ","), format(d[dir == "Down", .N], big.mark = ","),
                            FDR, lfc_thr),
         x = "log2 fold change (tumour / normal)", y = expression(-log[10]~FDR)) +
    theme_pub() + theme(plot.subtitle = element_text(size = 6.5))
}
pA <- volcano(de_rna, LFC_RNA, sprintf("mRNA (DESeq2, %d T vs %d N)", length(grep("_T$", rna_keys)), length(grep("_N$", rna_keys))))
pB <- volcano(de_prot, LFC_PROT, sprintf("Protein (limma, %d T vs %d N)", length(grep("_T$", prot_keys)), length(grep("_N$", prot_keys))))

cat_cols <- c(`Concordant up` = "#B2182B", `Concordant down` = "#2166AC", `mRNA only` = "#1B7837",
              `Protein only` = "#762A83", Opposite = "#E08214", Neither = "grey80")
both[, category := factor(category, levels = names(cat_cols))]
lab <- rbind(
  both[category == "Opposite"][order(-abs(log2FC_rna - 2 * log2FC_protein))][1:6],
  both[category == "Protein only"][order(-abs(log2FC_protein))][1:6],
  both[category == "mRNA only"][order(-abs(log2FC_rna))][1:6],
  both[category == "Concordant up"][order(-(log2FC_rna + 2 * log2FC_protein))][1:5],
  both[category == "Concordant down"][order(log2FC_rna + 2 * log2FC_protein)][1:5]
)[!is.na(gene)]
pC <- ggplot(both[order(category == "Neither", decreasing = TRUE)], aes(log2FC_rna, log2FC_protein, colour = category)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey50") + geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey50") +
  geom_point(size = 0.4, alpha = 0.7) +
  geom_text_repel(data = lab, aes(label = gene), size = 2, colour = "black", max.overlaps = 40,
                  segment.size = 0.2, min.segment.length = 0) +
  scale_colour_manual(values = cat_cols, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 1.8))) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.4, size = 2.4,
           label = sprintf("Spearman rho = %.2f\n%s genes", fc_cor$estimate, format(nrow(both), big.mark = ","))) +
  labs(title = "Tumour vs normal: mRNA against protein", x = "mRNA log2 fold change", y = "Protein log2 fold change", colour = NULL) +
  theme_pub() + theme(legend.position = "right")

cnt <- both[category != "Neither", .N, by = category]
pD <- ggplot(cnt, aes(reorder(category, N), N, fill = category)) + geom_col(width = 0.7) +
  geom_text(aes(label = N), hjust = -0.15, size = 2.4) + coord_flip() +
  scale_fill_manual(values = cat_cols, guide = "none") + scale_y_continuous(expand = expansion(c(0, 0.2))) +
  labs(title = "Genes per category", x = NULL, y = "Genes") + theme_pub()

fig4 <-(pA | pB) / ((pC | pD) + plot_layout(widths = c(1.7, 1))) + plot_annotation(tag_levels = "A")
save_fig(fig4, "Fig4_differential_expression", width = 200, height = 175)

# Placed after Figure 4 so that the random-number stream used by ggrepel for
# the figure labels is the same as without this section.
# ---------------------------------------------------------------- sensitivity: tumour-NAT pairing in RNA
# Every RNA NAT comes from a patient whose tumour is also in the RNA set, but the
# primary DESeq2 model treats all samples as independent. (a) limma-voom on the
# same samples and genes, without and with patient as a random block (two voom
# passes with duplicateCorrelation, as in the limma guide); comparing these two
# isolates the effect of pairing from the change of method. voom uses the
# DESeq2 median-of-ratios size factors as effective library sizes, so both RNA
# models share one normalisation. (b) DESeq2 on the tumour-NAT pairs alone,
# ~ tissue against ~ patient + tissue.
run_voom <- function(keys, blocked) {
  x <- rna_filtered(keys)
  sf <- estimateSizeFactorsForMatrix(x$counts)
  lib <- sf * exp(mean(log(colSums(x$counts))))
  design <- model.matrix(~tissue, x$cd)
  corr <- NA_real_
  v <- voom(x$counts, design, lib.size = lib)
  if (blocked) {
    b <- x$cd$patient
    corr <- duplicateCorrelation(v, design, block = b)$consensus.correlation
    v <- voom(x$counts, design, lib.size = lib, block = b, correlation = corr)
    corr <- duplicateCorrelation(v, design, block = b)$consensus.correlation
    fit <- lmFit(v, design, block = b, correlation = corr)
  } else {
    fit <- lmFit(v, design)
  }
  fit <- eBayes(fit)
  tt <- as.data.table(topTable(fit, coef = "tissueTumor", number = Inf, sort.by = "none"), keep.rownames = "gene_id")
  list(tab = tt[, .(gene_id, logFC, AveExpr, t, padj = adj.P.Val)], corr = corr)
}
calls <- function(lfc, padj, thr = LFC_RNA) fifelse(!is.na(padj) & padj < FDR & lfc > thr, "up",
                                                    fifelse(!is.na(padj) & padj < FDR & lfc < -thr, "down", "ns"))
jaccard <- function(a, b) { a <- a != "ns"; b <- b != "ns"; sum(a & b) / sum(a | b) }
rho <- function(a, b) round(cor(a, b, method = "spearman", use = "complete.obs"), 3)

vo  <- run_voom(rna_keys, blocked = FALSE)
vob <- run_voom(rna_keys, blocked = TRUE)
pg <- merge(de_rna[, .(gene_id, gene, log2FC_deseq2 = log2FC, log2FC_mle_deseq2 = log2FC_mle, stat_deseq2 = stat, padj_deseq2 = padj)],
            vo$tab[, .(gene_id, AveExpr, log2FC_voom = logFC, t_voom = t, padj_voom = padj)], by = "gene_id")
pg <- merge(pg, vob$tab[, .(gene_id, log2FC_voom_block = logFC, t_voom_block = t, padj_voom_block = padj)], by = "gene_id")
pg[, `:=`(call_deseq2 = calls(log2FC_deseq2, padj_deseq2), call_voom = calls(log2FC_voom, padj_voom),
          call_voom_block = calls(log2FC_voom_block, padj_voom_block))]
log_step(sprintf("voom: patient block correlation %.3f; DE genes %d unblocked, %d blocked",
                 vob$corr, pg[call_voom != "ns", .N], pg[call_voom_block != "ns", .N]))

# (b) the tumour-NAT pairs alone
pair_cases <- an[has_rna == 1, .(paired = all(c("Tumor", "Normal") %in% tissue)), by = case_submitter_id][paired == TRUE, case_submitter_id]
pair_keys <- an[has_rna == 1 & case_submitter_id %in% pair_cases, sample_key]
dp_u <- run_deseq(pair_keys, "pairs, ~ tissue")
dp_p <- run_deseq(pair_keys, "pairs, ~ patient + tissue", design = ~ patient + tissue)
pp <- merge(dp_u[, .(gene_id, log2FC_pairs = log2FC, stat_pairs = stat, padj_pairs = padj)],
            dp_p[, .(gene_id, log2FC_pairs_paired = log2FC, stat_pairs_paired = stat, padj_pairs_paired = padj)], by = "gene_id")
pp[, `:=`(call_pairs = calls(log2FC_pairs, padj_pairs), call_pairs_paired = calls(log2FC_pairs_paired, padj_pairs_paired))]
pg <- merge(pg, pp, by = "gene_id", all.x = TRUE)
setorder(pg, padj_deseq2, gene)
write_tsv(pg, "de_rna_pairing_genes.tsv")

# Downstream: RNA-protein categories and Hallmark GSEA with each RNA model
cats <- function(lfc, stat, padj) {
  r <- pg[, .(gene, baseMean = AveExpr, log2FC = lfc, stat = stat, padj = padj)]
  merge_layers(r, de_prot)
}
both_v  <- cats(pg$log2FC_voom, pg$t_voom, pg$padj_voom)
both_vb <- cats(pg$log2FC_voom_block, pg$t_voom_block, pg$padj_voom_block)
gs <- fread(path_proc("genesets.tsv.gz"))[collection == "HALLMARK"]
hall <- split(gs$gene_symbol, gs$gs_name)
hallmark <- function(stat, label) {
  ok <- !is.na(stat)
  set.seed(CFG$survival$seed)
  r <- as.data.table(fgsea(hall, sort(setNames(stat[ok], pg$gene[ok])), minSize = 15, maxSize = 500,
                           eps = 0, nPermSimple = 10000))
  setnames(r[, .(pathway, NES, padj)], c("pathway", paste0("NES_", label), paste0("padj_", label)))
}
hg <- Reduce(function(a, b) merge(a, b, by = "pathway"),
             list(hallmark(pg$stat_deseq2, "deseq2"), hallmark(pg$t_voom, "voom"), hallmark(pg$t_voom_block, "voom_block")))
write_tsv(hg[order(padj_voom_block)], "de_rna_pairing_hallmark.tsv")
sig_dir <- function(n, p) fifelse(p < FDR, sign(n), 0)
key_sets <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_MYC_TARGETS_V1")

n_up <- function(cl) sum(cl == "up"); n_dn <- function(cl) sum(cl == "down")
psumm <- data.table(metric = c(
    "RNA tumours", "RNA NATs", "tumour-NAT pairs (RNA)", "genes compared",
    "voom patient block consensus correlation",
    "DESeq2 up", "DESeq2 down", "voom up", "voom down", "voom blocked up", "voom blocked down",
    "DE calls kept after blocking (% of voom calls)", "DE calls gained by blocking",
    "Jaccard of DE calls, voom vs voom blocked", "Jaccard of DE calls, DESeq2 vs voom blocked",
    "Spearman rho of t, voom vs voom blocked", "Spearman rho of log2FC, voom vs voom blocked",
    "Spearman rho of log2FC, DESeq2 (MLE) vs voom blocked",
    "pairs only: DESeq2 ~ tissue DE genes", "pairs only: DESeq2 ~ patient + tissue DE genes",
    "pairs only: Jaccard of DE calls, unpaired vs paired", "pairs only: Spearman rho of Wald statistics, unpaired vs paired",
    "RNA-protein fold-change rho, voom", "RNA-protein fold-change rho, voom blocked",
    "Opposite-direction genes, voom", "Opposite-direction genes, voom blocked",
    "Hallmark NES rho, voom vs voom blocked", "Hallmark NES rho, DESeq2 vs voom blocked",
    "Hallmark sets significant, DESeq2", "Hallmark sets significant, voom blocked",
    "Hallmark sets changing significance or direction, voom vs voom blocked",
    "OXPHOS NES (voom blocked)", "MYC targets V1 NES (voom blocked)"),
  value = c(
    sum(an[has_rna == 1, tissue] == "Tumor"), sum(an[has_rna == 1, tissue] == "Normal"), length(pair_cases), nrow(pg),
    round(vob$corr, 3),
    n_up(pg$call_deseq2), n_dn(pg$call_deseq2), n_up(pg$call_voom), n_dn(pg$call_voom),
    n_up(pg$call_voom_block), n_dn(pg$call_voom_block),
    round(100 * pg[call_voom != "ns", mean(call_voom_block == call_voom)], 1),
    pg[call_voom == "ns" & call_voom_block != "ns", .N],
    round(jaccard(pg$call_voom, pg$call_voom_block), 3), round(jaccard(pg$call_deseq2, pg$call_voom_block), 3),
    rho(pg$t_voom, pg$t_voom_block), rho(pg$log2FC_voom, pg$log2FC_voom_block),
    rho(pg$log2FC_mle_deseq2, pg$log2FC_voom_block),
    pp[call_pairs != "ns", .N], pp[call_pairs_paired != "ns", .N],
    round(jaccard(pp$call_pairs, pp$call_pairs_paired), 3), rho(pp$stat_pairs, pp$stat_pairs_paired),
    rho(both_v$log2FC_rna, both_v$log2FC_protein), rho(both_vb$log2FC_rna, both_vb$log2FC_protein),
    both_v[category == "Opposite", .N], both_vb[category == "Opposite", .N],
    rho(hg$NES_voom, hg$NES_voom_block), rho(hg$NES_deseq2, hg$NES_voom_block),
    hg[padj_deseq2 < FDR, .N], hg[padj_voom_block < FDR, .N],
    hg[sig_dir(NES_voom, padj_voom) != sig_dir(NES_voom_block, padj_voom_block), .N],
    round(hg[pathway == key_sets[1], NES_voom_block], 2), round(hg[pathway == key_sets[2], NES_voom_block], 2)))
pcat <- rbind(both_v[, .(metric = paste0("category, voom: ", category), value = .N), by = category],
              both_vb[, .(metric = paste0("category, voom blocked: ", category), value = .N), by = category])
psumm <- rbind(psumm, pcat[order(category, metric), .(metric, value = as.numeric(value))])
write_tsv(psumm, "de_rna_pairing_summary.tsv")
print(psumm)
