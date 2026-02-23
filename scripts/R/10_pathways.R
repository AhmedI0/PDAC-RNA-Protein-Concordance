# Step 10 - Pathway-level view of tumour vs normal and of RNA-protein agreement.
#
# 1. Pre-ranked GSEA (fgsea, multilevel) on the DESeq2 Wald statistic and the
#    limma moderated t of step 09: Hallmark (50) and Reactome (15-500 genes).
# 2. GSVA (Hanzelmann et al. 2013) per-sample Hallmark scores in the matched
#    tumours, computed separately from mRNA (VST) and protein (log2 ratios)
#    restricted to the genes measured in both layers. The correlation of a
#    pathway's mRNA score with its protein score across tumours is compared with
#    the median gene-level correlation of its members: aggregation into pathways
#    averages out gene-specific noise and post-transcriptional effects.
#
#   Rscript scripts/R/10_pathways.R

source("R/utils.R")
suppressPackageStartupMessages({
  library(fgsea)
  library(GSVA)
  library(patchwork)
  library(ggrepel)
  library(BiocParallel)
})
set.seed(CFG$survival$seed)

gs <- fread(path_proc("genesets.tsv.gz"))
hall <- split(gs[collection == "HALLMARK", gene_symbol], gs[collection == "HALLMARK", gs_name])
reac <- split(gs[collection == "REACTOME", gene_symbol], gs[collection == "REACTOME", gs_name])
pretty <- function(x) {
  x <- sub("^HALLMARK_|^REACTOME_", "", x)
  x <- gsub("_", " ", tolower(x))
  paste0(toupper(substr(x, 1, 1)), substring(x, 2))
}

# ============================================================== fgsea
de_rna  <- fread(path_tab("de_rna_tumour_vs_normal.tsv"))
de_prot <- fread(path_tab("de_protein_tumour_vs_normal.tsv"))
rank_rna  <- sort(setNames(de_rna[!is.na(stat), stat], de_rna[!is.na(stat), gene]))
rank_prot <- sort(setNames(de_prot$t, de_prot$gene))

run_fgsea <- function(sets, ranks, label, layer) {
  set.seed(CFG$survival$seed)
  r <- fgsea(pathways = sets, stats = ranks, minSize = 15, maxSize = 500, eps = 0, nPermSimple = 10000)
  r <- as.data.table(r)[, .(pathway, size, NES, pval, padj, leadingEdge = vapply(leadingEdge, paste, "", collapse = ";"))]
  r[, `:=`(collection = label, layer = layer)]
  r[]
}
gsea <- rbindlist(list(
  run_fgsea(hall, rank_rna, "HALLMARK", "RNA"), run_fgsea(hall, rank_prot, "HALLMARK", "Protein"),
  run_fgsea(reac, rank_rna, "REACTOME", "RNA"), run_fgsea(reac, rank_prot, "REACTOME", "Protein")
))
write_tsv(gsea[order(collection, layer, padj)], "gsea_tumour_vs_normal.tsv")

gw <- dcast(gsea, collection + pathway ~ layer, value.var = c("NES", "padj"))
gw[, agreement := fifelse(padj_RNA < 0.05 & padj_Protein < 0.05 & sign(NES_RNA) == sign(NES_Protein), "both, same direction",
                  fifelse(padj_RNA < 0.05 & padj_Protein < 0.05, "both, opposite",
                  fifelse(padj_RNA < 0.05, "RNA only", fifelse(padj_Protein < 0.05, "protein only", "neither"))))]
write_tsv(gw, "gsea_rna_vs_protein.tsv")
# Many Reactome sets re-use the same ribosomal, proteasomal and respiratory-chain
# subunits; their protein-level behaviour (diluted in collagen-rich tumours when
# TMT data are normalised per total protein) then dominates the set-level comparison.
# Sensitivity: drop sets in which >= 25% of members belong to these complexes.
cc <- gs[collection == "GO_CC"]
machinery <- unique(cc[gs_name %in% c("GOCC_CYTOSOLIC_RIBOSOME", "GOCC_MITOCHONDRIAL_RIBOSOME", "GOCC_PROTEASOME_COMPLEX",
                                      "GOCC_RESPIRATORY_CHAIN_COMPLEX", "GOCC_SPLICEOSOMAL_COMPLEX"), gene_symbol])
frac_mach <- sapply(reac, function(g) mean(g %in% machinery))
gw[collection == "REACTOME", machinery_fraction := frac_mach[pathway]]
nes_sum <- data.table(
  set = c("Hallmark", "Reactome, all", "Reactome, <25% complex machinery"),
  n_sets = c(gw[collection == "HALLMARK", .N], gw[collection == "REACTOME", .N],
             gw[collection == "REACTOME" & machinery_fraction < 0.25, .N]),
  spearman_NES = c(gw[collection == "HALLMARK", cor(NES_RNA, NES_Protein, method = "spearman", use = "complete.obs")],
                   gw[collection == "REACTOME", cor(NES_RNA, NES_Protein, method = "spearman", use = "complete.obs")],
                   gw[collection == "REACTOME" & machinery_fraction < 0.25,
                      cor(NES_RNA, NES_Protein, method = "spearman", use = "complete.obs")]))
write_tsv(gw, "gsea_rna_vs_protein.tsv")
write_tsv(nes_sum, "gsea_nes_agreement.tsv")
print(nes_sum)
log_step(sprintf("Hallmark agreement: %s",
                 paste(names(table(gw[collection == "HALLMARK", agreement])), table(gw[collection == "HALLMARK", agreement]), collapse = ", ")))

# ============================================================== GSVA in matched tumours
an <- fread(path_proc("analysis_samples.tsv"))
T <- an[tissue == "Tumor" & matched == TRUE, sample_key]
rna  <- read_matrix("rna_vst.tsv.gz")[, T]
prot <- read_matrix("protein_norm.tsv.gz")[, T]
conc <- fread(path_tab("concordance_gene.tsv"))
common <- conc$gene                                   # genes in both layers, protein in >= 50% of tumours
prot_c <- prot[common, ]
# GSVA needs complete data: impute the few missing protein values per gene with the gene's minimum
# observed value minus 0.5 (missing-not-at-random: low-abundance), only for this score.
for (i in which(rowSums(is.na(prot_c)) > 0)) {
  v <- prot_c[i, ]; v[is.na(v)] <- min(v, na.rm = TRUE) - 0.5; prot_c[i, ] <- v
}
hall_c <- lapply(hall, intersect, common)
hall_c <- hall_c[lengths(hall_c) >= 10]
gsva_score <- function(m) gsva(gsvaParam(m, hall_c, minSize = 10, maxSize = 500, kcdf = "Gaussian"),
                               verbose = FALSE, BPPARAM = SerialParam())
g_rna <- gsva_score(rna[common, ])
g_pro <- gsva_score(prot_c)
fwrite(data.table(pathway = rownames(g_rna), g_rna), path_proc("gsva_hallmark_rna.tsv.gz"), sep = "\t")
fwrite(data.table(pathway = rownames(g_pro), g_pro), path_proc("gsva_hallmark_protein.tsv.gz"), sep = "\t")

pw <- rbindlist(lapply(rownames(g_rna), function(p) {
  ct <- cor.test(g_rna[p, ], g_pro[p, ], method = "spearman", exact = FALSE)
  data.table(pathway = p, n_genes = length(hall_c[[p]]), rho_pathway = unname(ct$estimate), p = ct$p.value,
             median_rho_genes = median(conc[gene %in% hall_c[[p]], rho_tumour], na.rm = TRUE))
}))
pw[, fdr := p.adjust(p, "BH")]
pw[, gain := rho_pathway - median_rho_genes]
write_tsv(pw[order(-rho_pathway)], "pathway_level_concordance.tsv")
log_step(sprintf("GSVA pathway-level rho: median %.2f vs median gene-level %.2f (paired Wilcoxon p=%.2g)",
                 median(pw$rho_pathway), median(pw$median_rho_genes),
                 wilcox.test(pw$rho_pathway, pw$median_rho_genes, paired = TRUE)$p.value))

# ============================================================== figure 5
hw <- gw[collection == "HALLMARK"]
hw[, label := pretty(pathway)]
hw[, sig := factor(agreement, levels = c("both, same direction", "RNA only", "protein only", "both, opposite", "neither"))]
pA <- ggplot(hw, aes(NES_RNA, NES_Protein, colour = sig)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey50") + geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey50") +
  geom_abline(linetype = 3, linewidth = 0.3) +
  geom_point(size = 1.4) +
  # label a readable subset: every opposite set and the strongest sets of each other class
  geom_text_repel(data = rbind(hw[sig == "both, opposite"],
                               hw[sig == "both, same direction"][order(-abs(NES_RNA + NES_Protein))][1:5],
                               hw[sig == "protein only"][order(-abs(NES_Protein))][1:3],
                               hw[sig == "RNA only"][order(-abs(NES_RNA))][1:3]),
                  aes(label = label), size = 1.9, colour = "black", max.overlaps = Inf, box.padding = 0.35,
                  min.segment.length = 0, segment.size = 0.2, seed = 1) +
  scale_colour_manual(values = c(`both, same direction` = "#B2182B", `RNA only` = PAL$layer[["RNA"]],
                                 `protein only` = PAL$layer[["Protein"]], `both, opposite` = "#E08214", neither = "grey70"),
                      drop = FALSE) +
  annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.6, size = 2.3,
           label = sprintf("Spearman rho = %.2f", nes_sum[set == "Hallmark", spearman_NES])) +
  labs(title = "Hallmark GSEA, tumour vs normal", x = "NES (mRNA)", y = "NES (protein)", colour = "FDR < 5%") +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 1.6))) +
  theme_pub() + theme(legend.position = "bottom", legend.text = element_text(size = 5.5),
                      legend.title = element_text(size = 6), legend.margin = margin(0, 0, 0, 0))

topR <- gw[collection == "REACTOME" & (padj_RNA < 0.05 | padj_Protein < 0.05)]
topR <- topR[order(-pmax(abs(NES_RNA), abs(NES_Protein)))][1:24]
tl <- melt(topR[, .(pathway, NES_RNA, NES_Protein, padj_RNA, padj_Protein)], id.vars = "pathway",
           measure.vars = list(NES = c("NES_RNA", "NES_Protein"), padj = c("padj_RNA", "padj_Protein")),
           variable.name = "layer")
tl[, layer := factor(c("mRNA", "Protein")[layer], levels = c("mRNA", "Protein"))]
tl[, label := substr(pretty(pathway), 1, 48)]
tl[, label := factor(label, levels = unique(tl[layer == "mRNA"][order(NES), label]))]
pB <- ggplot(tl, aes(layer, label, fill = NES)) + geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(padj < 0.05, "*", "")), size = 2.6, vjust = 0.75) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-3.5, 3.5), oob = scales::squish) +
  labs(title = "Reactome (top 24)", x = NULL, y = NULL) + theme_pub() +
  theme(axis.text.y = element_text(size = 5.2), axis.line = element_blank(), axis.ticks = element_blank())

pw[, label := pretty(pathway)]
pC <- ggplot(pw, aes(median_rho_genes, rho_pathway)) +
  geom_abline(linetype = 3, linewidth = 0.3) +
  geom_point(size = 1.3, colour = PAL$layer[["Protein"]]) +
  geom_text_repel(data = pw[order(-gain)][c(1:6, (.N - 2):.N)], aes(label = label), size = 1.9, max.overlaps = 20, segment.size = 0.2) +
  labs(title = "Pathway vs gene-level agreement", x = "Median gene-level rho (members)",
       y = "Pathway-score rho (GSVA mRNA vs protein)") +
  coord_cartesian(xlim = c(0, 1), ylim = c(-0.1, 1)) + theme_pub()

ex <- pw[order(-rho_pathway)][c(1, round(.N / 2), .N), pathway]
exd <- rbindlist(lapply(ex, function(p) data.table(pathway = pretty(p), rna = g_rna[p, ], protein = g_pro[p, ])))
exd[, pathway := factor(pathway, levels = pretty(ex))]
pD <- ggplot(exd, aes(rna, protein)) + geom_point(size = 0.6, colour = "grey25") +
  geom_smooth(method = "lm", formula = y ~ x, linewidth = 0.5, colour = "#B2182B", se = FALSE) +
  facet_wrap(~pathway, scales = "free", ncol = 1) +
  labs(title = "Examples (high / median / low)", x = "GSVA score (mRNA)", y = "GSVA score (protein)") +
  theme_pub() + theme(strip.text = element_text(size = 6))

fig5 <- (pA | pB) / (pC | pD) + plot_layout(heights = c(1, 1)) + plot_annotation(tag_levels = "A")
save_fig(fig5, "Fig5_pathways", width = 200, height = 190)
