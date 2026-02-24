# Step 11 - Molecular subtypes from mRNA and from protein.
#
# 1. PurIST (Rashid et al. 2020): single-sample classical / basal-like call from
#    8 weighted top-scoring gene pairs on log2(TPM+1). Reference subtype.
# 2. Moffitt 2015 tumour-specific signatures (25 basal-like + 25 classical genes):
#    a) continuous score per layer = mean z(basal) - mean z(classical);
#    b) consensus clustering (k = 2) on the same genes in each layer, as in the
#       original paper, giving one categorical call per layer.
#    Agreement: Cohen's kappa (RNA vs protein calls), AUC of each layer's score
#    against PurIST.
# 3. Unsupervised consensus clustering (top 1,500 most variable features per
#    layer, k = 2..6, k chosen by PAC); agreement between layers by adjusted
#    Rand index; clusters described by purity and Hallmark GSVA scores.
#
#   Rscript scripts/R/11_subtypes.R

source("R/utils.R")
suppressPackageStartupMessages({
  library(ConsensusClusterPlus)
  library(aricode)
  library(patchwork)
  library(matrixStats)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})
set.seed(CFG$survival$seed)

an <- fread(path_proc("analysis_samples.tsv"))
T <- an[tissue == "Tumor" & matched == TRUE, sample_key]
rna  <- read_matrix("rna_vst.tsv.gz")[, T]
prot <- read_matrix("protein_norm.tsv.gz")[, T]
tpm  <- read_matrix("rna_log2tpm.tsv.gz")[, T]
pur  <- fread(path_tab("purity_estimates.tsv"))

# map published symbols to GENCODE v36 symbols when they were renamed
to_current <- function(sym, universe) {
  out <- sym
  miss <- !sym %in% universe
  if (any(miss)) {
    al <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = sym[miss], keytype = "ALIAS", columns = "SYMBOL"))
    al <- al[al$SYMBOL %in% universe, ]
    out[miss] <- al$SYMBOL[match(sym[miss], al$ALIAS)]
  }
  out
}

# ============================================================== PurIST
pm <- fread(file.path(ROOT, "resources", "purist_tsp_model.tsv"))
b0 <- pm[pair == 0, beta]
pm <- pm[pair > 0 & beta != 0]
pm[, `:=`(a = to_current(gene_a, rownames(tpm)), b = to_current(gene_b, rownames(tpm)))]
stopifnot(!anyNA(pm$a), !anyNA(pm$b))
ind <- sapply(seq_len(nrow(pm)), function(i) as.numeric(tpm[pm$a[i], ] > tpm[pm$b[i], ]))
purist <- data.table(sample_key = T, purist_prob = plogis(b0 + drop(ind %*% pm$beta)))
purist[, purist_call := fifelse(purist_prob >= 0.5, "Basal-like", "Classical")]
log_step(sprintf("PurIST: %d pairs used; basal-like %d / %d (%.0f%%)", nrow(pm),
                 purist[purist_call == "Basal-like", .N], nrow(purist), 100 * purist[, mean(purist_call == "Basal-like")]))

# ============================================================== Moffitt signatures
sig <- fread(file.path(ROOT, "resources", "pdac_subtype_signatures.tsv"))
basal <- sig[signature == "Moffitt.Basal.25", gene_symbol]
class <- sig[signature == "Moffitt.Classical.25", gene_symbol]
cov_tab <- data.table(
  signature = c("basal-like", "classical"),
  n_genes = c(length(basal), length(class)),
  in_rna = c(sum(to_current(basal, rownames(rna)) %in% rownames(rna)), sum(to_current(class, rownames(rna)) %in% rownames(rna))),
  in_protein = c(sum(to_current(basal, rownames(prot)) %in% rownames(prot)), sum(to_current(class, rownames(prot)) %in% rownames(prot))))
write_tsv(cov_tab, "subtype_signature_coverage.tsv")
print(cov_tab)

sig_score <- function(m, up, dn, min_frac = 0.5) {
  up <- na.omit(to_current(up, rownames(m))); dn <- na.omit(to_current(dn, rownames(m)))
  up <- up[up %in% rownames(m)]; dn <- dn[dn %in% rownames(m)]
  keep <- function(g) g[rowMeans(!is.na(m[g, , drop = FALSE])) >= min_frac]
  up <- keep(up); dn <- keep(dn)
  z <- (m - rowMeans(m, na.rm = TRUE)) / rowSds(m, na.rm = TRUE)
  list(score = colMeans(z[up, , drop = FALSE], na.rm = TRUE) - colMeans(z[dn, , drop = FALSE], na.rm = TRUE),
       genes = c(up, dn), z = z[c(up, dn), , drop = FALSE], up = up, dn = dn)
}
s_rna <- sig_score(rna, basal, class)
s_pro <- sig_score(prot, basal, class)

cc_call <- function(s, label) {
  z <- s$z
  z[is.na(z)] <- 0
  dir.create(file.path(ROOT, ".cache", "tmp", "cc"), showWarnings = FALSE, recursive = TRUE)
  # maxK = 3 because ConsensusClusterPlus' tracking plot fails for maxK = 2; k = 2 is used
  cc <- ConsensusClusterPlus(z, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, clusterAlg = "hc",
                             innerLinkage = "average", finalLinkage = "average", distance = "pearson",
                             seed = CFG$survival$seed, plot = NULL, verbose = FALSE)
  cl <- cc[[2]]$consensusClass
  basal_cl <- names(which.max(tapply(s$score, cl, mean)))
  fifelse(as.character(cl) == basal_cl, "Basal-like", "Classical")
}
sub <- data.table(sample_key = T, score_rna = s_rna$score, score_protein = s_pro$score,
                  moffitt_rna = cc_call(s_rna, "RNA"), moffitt_protein = cc_call(s_pro, "Protein"))
sub <- merge(purist, sub, by = "sample_key")

kappa2 <- function(a, b) {
  t <- table(factor(a, levels = c("Classical", "Basal-like")), factor(b, levels = c("Classical", "Basal-like")))
  po <- sum(diag(t)) / sum(t); pe <- sum(rowSums(t) * colSums(t)) / sum(t)^2
  (po - pe) / (1 - pe)
}
auc <- function(score, truth) {
  r <- rank(score); n1 <- sum(truth); n0 <- sum(!truth)
  (sum(r[truth]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
agree <- data.table(
  comparison = c("Moffitt RNA vs Moffitt protein", "PurIST vs Moffitt RNA", "PurIST vs Moffitt protein"),
  kappa = c(kappa2(sub$moffitt_rna, sub$moffitt_protein), kappa2(sub$purist_call, sub$moffitt_rna),
            kappa2(sub$purist_call, sub$moffitt_protein)),
  percent_agreement = 100 * c(mean(sub$moffitt_rna == sub$moffitt_protein), mean(sub$purist_call == sub$moffitt_rna),
                              mean(sub$purist_call == sub$moffitt_protein)))
score_cor <- cor.test(sub$score_rna, sub$score_protein, method = "spearman", exact = FALSE)
aucs <- data.table(layer = c("RNA", "Protein"),
                   auc_vs_purist = c(auc(sub$score_rna, sub$purist_call == "Basal-like"),
                                     auc(sub$score_protein, sub$purist_call == "Basal-like")))
write_tsv(agree, "subtype_agreement.tsv")
write_tsv(aucs, "subtype_score_auc.tsv")
print(agree); print(aucs)
log_step(sprintf("Moffitt score RNA vs protein: Spearman %.2f", score_cor$estimate))

# ============================================================== unsupervised clustering
pac <- function(cm, lower = 0.1, upper = 0.9) {
  v <- cm[lower.tri(cm)]
  mean(v > lower & v < upper)
}
consensus_layer <- function(m, label, n_feat = 1500) {
  m <- m[rowMeans(!is.na(m)) == 1, ]
  m <- m[order(rowMads(m), decreasing = TRUE)[1:n_feat], ]
  m <- m - rowMedians(m)
  cc <- ConsensusClusterPlus(m, maxK = 6, reps = 1000, pItem = 0.8, pFeature = 1, clusterAlg = "km",
                             distance = "euclidean", seed = CFG$survival$seed, plot = NULL, verbose = FALSE)
  pacs <- sapply(2:6, function(k) pac(cc[[k]]$consensusMatrix))
  k <- (2:6)[which.min(pacs)]
  list(k = k, pac = data.table(layer = label, k = 2:6, PAC = pacs), cl = cc[[k]]$consensusClass,
       all = lapply(2:6, function(j) cc[[j]]$consensusClass))
}
uc_r <- consensus_layer(rna, "RNA")
uc_p <- consensus_layer(prot, "Protein")
write_tsv(rbind(uc_r$pac, uc_p$pac), "consensus_clustering_pac.tsv")
sub[, cluster_rna := paste0("R", uc_r$cl[sample_key])]
sub[, cluster_protein := paste0("P", uc_p$cl[sample_key])]
ari <- data.table(k = 2:6, ARI_rna_vs_protein = sapply(1:5, function(i) ARI(uc_r$all[[i]][T], uc_p$all[[i]][T])))
ari_best <- ARI(sub$cluster_rna, sub$cluster_protein)
write_tsv(ari, "consensus_clustering_ari.tsv")
log_step(sprintf("Consensus clustering: RNA k=%d, protein k=%d (by PAC); ARI between layers at chosen k = %.2f",
                 uc_r$k, uc_p$k, ari_best))
print(ari)

sub <- merge(sub, pur[, .(sample_key, kras_cellularity, estimate_score, stromal_score, immune_score)], by = "sample_key", all.x = TRUE)
sub <- merge(sub, an[, .(sample_key, case_submitter_id, acinar_rna)], by = "sample_key")
write_tsv(sub, "subtype_calls.tsv")

# cluster description with GSVA Hallmark scores (protein) and composition
g_pro <- read_matrix("gsva_hallmark_protein.tsv.gz")[, T]
g_rna <- read_matrix("gsva_hallmark_rna.tsv.gz")[, T]
describe <- function(g, cl, label) {
  rbindlist(lapply(rownames(g), function(p) {
    kw <- kruskal.test(g[p, ] ~ factor(cl))
    data.table(layer = label, pathway = p, p = kw$p.value,
               t(tapply(g[p, ], cl, mean)))
  }), fill = TRUE)[, fdr := p.adjust(p, "BH")][order(p)]
}
write_tsv(describe(g_pro, sub[match(T, sample_key), cluster_protein], "Protein"), "cluster_hallmark_protein.tsv")
write_tsv(describe(g_rna, sub[match(T, sample_key), cluster_rna], "RNA"), "cluster_hallmark_rna.tsv")

# ============================================================== figure 6
pA <- ggplot(sub, aes(score_rna, score_protein, colour = purist_call)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey60") + geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey60") +
  geom_point(size = 1.1, alpha = 0.85) + scale_colour_manual(values = PAL$subtype) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.3, size = 2.4,
           label = sprintf("Spearman rho = %.2f\nAUC vs PurIST: RNA %.2f, protein %.2f", score_cor$estimate,
                           aucs[layer == "RNA", auc_vs_purist], aucs[layer == "Protein", auc_vs_purist])) +
  labs(title = "Moffitt basal-like minus classical score", x = "mRNA score", y = "Protein score", colour = "PurIST (RNA)") +
  theme_pub() + theme(legend.position = "bottom")

ct <- sub[, .N, by = .(moffitt_rna, moffitt_protein)]
pB <- ggplot(ct, aes(moffitt_rna, moffitt_protein, fill = N)) + geom_tile(colour = "white") +
  geom_text(aes(label = N), size = 3) + scale_fill_gradient(low = "grey95", high = "#762A83", guide = "none") +
  labs(title = sprintf("Moffitt consensus calls\nkappa = %.2f", agree[1, kappa]), x = "mRNA", y = "Protein") + theme_pub()

ord <- sub[order(purist_prob), sample_key]
hm <- function(s, label) {
  d <- as.data.table(s$z[, ord], keep.rownames = "gene")
  d <- melt(d, id.vars = "gene", variable.name = "sample_key", value.name = "z")
  d[, set := fifelse(gene %in% s$up, "Basal-like", "Classical")]
  d[, layer := label]
  d
}
hd <- rbind(hm(s_rna, "mRNA"), hm(s_pro, "Protein"))
# PurIST call as a top row of each facet, so it lines up with the tumours of both layers
pur_row <- sub[, .(gene = "PurIST call", sample_key, z = fifelse(purist_call == "Basal-like", 2.5, -2.5), set = "PurIST")]
hd <- rbind(hd, pur_row[, layer := "mRNA"], copy(pur_row)[, layer := "Protein"], use.names = TRUE)
hd[, sample_key := factor(sample_key, levels = ord)]
hd[, set := factor(set, levels = c("PurIST", "Basal-like", "Classical"))]
hd[, gene := factor(gene, levels = unique(c("PurIST call", s_rna$up, s_rna$dn, s_pro$up, s_pro$dn)))]
pC <- ggplot(hd, aes(sample_key, gene, fill = pmax(pmin(z, 2.5), -2.5))) + geom_raster() +
  facet_grid(set ~ layer, scales = "free_y", space = "free_y",
             labeller = labeller(set = c(PurIST = "", `Basal-like` = "Basal-like", Classical = "Classical"))) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", na.value = "grey85", name = "z") +
  labs(title = "Moffitt genes, tumours ordered by PurIST probability",
       subtitle = "top row: PurIST call (red basal-like, blue classical)", x = NULL, y = NULL) + theme_pub() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 3.6),
        axis.line = element_blank(), panel.spacing = unit(1, "mm"), plot.subtitle = element_text(size = 6.5),
        strip.text.y = element_text(size = 6))

pacd <- rbind(uc_r$pac, uc_p$pac)
pacd <- merge(pacd, ari, by = "k")
pD <- ggplot(pacd, aes(k, PAC, colour = layer)) + geom_line(linewidth = 0.5) + geom_point(size = 1.2) +
  geom_line(data = ari, aes(k, ARI_rna_vs_protein), inherit.aes = FALSE, linetype = 2, linewidth = 0.5) +
  geom_point(data = ari, aes(k, ARI_rna_vs_protein), inherit.aes = FALSE, shape = 1, size = 1.3) +
  scale_colour_manual(values = c(RNA = PAL$layer[["RNA"]], Protein = PAL$layer[["Protein"]])) +
  annotate("text", x = 2.6, y = 0.02, label = "ARI, RNA vs protein", hjust = 0, size = 2.2) +
  labs(title = "Unsupervised clustering", subtitle = sprintf("chosen k: RNA %d, protein %d", uc_r$k, uc_p$k),
       x = "k", y = "PAC (solid) / ARI (dashed)", colour = NULL) +
  theme_pub() + theme(legend.position = c(0.75, 0.9), legend.background = element_blank(),
                      plot.subtitle = element_text(size = 6.5))

cd <- melt(sub[, .(cluster_protein, `Cellularity (KRAS VAF)` = kras_cellularity, `Stromal score` = stromal_score,
                   `Acinar score` = acinar_rna, `Moffitt score (protein)` = score_protein)],
           id.vars = "cluster_protein")
pE <- ggplot(cd, aes(cluster_protein, value, fill = cluster_protein)) +
  geom_boxplot(outlier.size = 0.3, width = 0.65, linewidth = 0.3) + facet_wrap(~variable, scales = "free_y", nrow = 1) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(title = "What separates the protein clusters", x = NULL, y = NULL) + theme_pub() +
  theme(strip.text = element_text(size = 6))

top <- (pA | pB | pD) + plot_layout(widths = c(1.4, 0.8, 0.8))
fig6 <- (top / free(pC) / pE) + plot_layout(heights = c(1, 1.25, 0.65)) + plot_annotation(tag_levels = "A")
save_fig(fig6, "Fig6_subtypes", width = 200, height = 230)
