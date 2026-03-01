# Step 12 - Clinical outcome: overall survival by subtype, by gene and pathway
# (mRNA vs protein) and by per-tumour RNA-protein agreement.
#
# Cohort: matched tumours with follow-up; OS censored at 5 years.
# Cox proportional hazards models adjusted for age, stage, resection margin and
# grade; proportional hazards checked with Schoenfeld residuals (cox.zph).
#
#   Rscript scripts/R/12_survival.R

source("R/utils.R")
source("R/survival_utils.R")
suppressPackageStartupMessages({
  library(ggsurvfit)
  library(patchwork)
  library(matrixStats)
})
set.seed(CFG$survival$seed)

an    <- fread(path_proc("analysis_samples.tsv"))
cases <- read_cases()
sv    <- survival_table(an, cases)
sub   <- fread(path_tab("subtype_calls.tsv"))
sw    <- fread(path_tab("concordance_sample.tsv"))
pur   <- fread(path_tab("purity_estimates.tsv"))
sv <- merge(sv, sub[, .(sample_key, purist_call, purist_prob, moffitt_rna, moffitt_protein, score_rna, score_protein,
                        cluster_rna, cluster_protein)], by = "sample_key")
sv <- merge(sv, sw[, .(sample_key, rho_samplewise)], by = "sample_key")
sv <- merge(sv, pur[, .(sample_key, estimate_score, kras_cellularity)], by = "sample_key")
for (v in c("purist_call", "moffitt_rna", "moffitt_protein")) sv[, (v) := factor(get(v), levels = c("Classical", "Basal-like"))]
log_step(sprintf("survival cohort: %d tumours, %d deaths within 5 years", nrow(sv), sum(sv$event)))

# ============================================================== clinical and subtype models
fits <- list(
  clinical = coxph(as.formula(paste("Surv(time, event) ~", CLIN_FORMULA)), sv),
  purist = coxph(as.formula(paste("Surv(time, event) ~ purist_call +", CLIN_FORMULA)), sv),
  moffitt_rna = coxph(as.formula(paste("Surv(time, event) ~ moffitt_rna +", CLIN_FORMULA)), sv),
  moffitt_protein = coxph(as.formula(paste("Surv(time, event) ~ moffitt_protein +", CLIN_FORMULA)), sv),
  score_rna = coxph(as.formula(paste("Surv(time, event) ~ scale(score_rna) +", CLIN_FORMULA)), sv),
  score_protein = coxph(as.formula(paste("Surv(time, event) ~ scale(score_protein) +", CLIN_FORMULA)), sv),
  concordance = coxph(as.formula(paste("Surv(time, event) ~ scale(rho_samplewise) + scale(estimate_score) +", CLIN_FORMULA)), sv),
  concordance_unadj = coxph(Surv(time, event) ~ scale(rho_samplewise), sv)
)
cox_tab <- rbindlist(lapply(names(fits), function(n) tidy_cox(fits[[n]], n)))
write_tsv(cox_tab, "survival_cox_models.tsv")
zph <- rbindlist(lapply(names(fits), function(n) {
  z <- cox.zph(fits[[n]])$table
  data.table(model = n, term = rownames(z), chisq = z[, "chisq"], p = z[, "p"])
}))
write_tsv(zph, "survival_ph_checks.tsv")
lr <- rbindlist(lapply(c("purist_call", "moffitt_rna", "moffitt_protein", "cluster_rna", "cluster_protein"), function(v) {
  s <- survdiff(as.formula(paste("Surv(time, event) ~", v)), sv)
  data.table(grouping = v, groups = length(s$n), chisq = s$chisq, p = pchisq(s$chisq, length(s$n) - 1, lower.tail = FALSE))
}))
write_tsv(lr, "survival_logrank.tsv")
print(cox_tab[grepl("purist|moffitt|score|rho", term)]); print(lr)

# ============================================================== gene-level Cox, both layers
T <- sv$sample_key
conc <- fread(path_tab("concordance_gene.tsv"))
rna  <- read_matrix("rna_vst.tsv.gz")[conc$gene, T]
prot <- read_matrix("protein_norm.tsv.gz")[conc$gene, T]
zs <- function(m) (m - rowMeans(m, na.rm = TRUE)) / rowSds(m, na.rm = TRUE)
rz <- zs(rna); pz <- zs(prot)
gene_cox <- function(z, label) {
  rbindlist(lapply(rownames(z), function(g) {
    x <- z[g, ]
    f <- tryCatch(coxph(Surv(sv$time, sv$event) ~ x + sv$stage), error = function(e) NULL)
    if (is.null(f)) return(NULL)
    s <- summary(f)$coefficients
    data.table(gene = g, layer = label, logHR = s[1, "coef"], se = s[1, "se(coef)"], p = s[1, "Pr(>|z|)"],
               n = sum(!is.na(x)))
  }))
}
gc_r <- gene_cox(rz, "RNA"); gc_p <- gene_cox(pz, "Protein")
gc <- rbind(gc_r, gc_p)
gc[, fdr := p.adjust(p, "BH"), by = layer]
gw <- dcast(gc, gene ~ layer, value.var = c("logHR", "p", "fdr"))
gw <- merge(gw, conc[, .(gene, rho_tumour)], by = "gene")
write_tsv(gc, "survival_gene_cox_long.tsv")
write_tsv(gw[order(pmin(p_RNA, p_Protein))], "survival_gene_cox.tsv")
hr_cor <- cor.test(gw$logHR_RNA, gw$logHR_Protein, method = "spearman", exact = FALSE)
# does prognostic agreement depend on RNA-protein concordance?
gw[, conc_bin := cut(rho_tumour, c(-1, 0.2, 0.4, 0.6, 1), labels = c("<0.2", "0.2-0.4", "0.4-0.6", ">0.6"))]
hr_by_bin <- gw[, .(genes = .N, hr_spearman = cor(logHR_RNA, logHR_Protein, method = "spearman")), by = conc_bin][order(conc_bin)]
write_tsv(hr_by_bin, "survival_hr_agreement_by_concordance.tsv")
gsum <- data.table(
  metric = c("genes tested", "RNA p<0.01", "protein p<0.01", "both p<0.01, same direction", "RNA FDR<0.1", "protein FDR<0.1",
             "Spearman logHR RNA vs protein"),
  value = c(nrow(gw), gw[p_RNA < 0.01, .N], gw[p_Protein < 0.01, .N],
            gw[p_RNA < 0.01 & p_Protein < 0.01 & sign(logHR_RNA) == sign(logHR_Protein), .N],
            gw[fdr_RNA < 0.1, .N], gw[fdr_Protein < 0.1, .N], round(hr_cor$estimate, 3)))
write_tsv(gsum, "survival_gene_cox_summary.tsv")
print(gsum); print(hr_by_bin)

# ============================================================== pathway-level Cox (GSVA)
g_r <- read_matrix("gsva_hallmark_rna.tsv.gz")[, T]
g_p <- read_matrix("gsva_hallmark_protein.tsv.gz")[, T]
pw_cox <- rbindlist(lapply(list(RNA = g_r, Protein = g_p), function(m) {
  rbindlist(lapply(rownames(m), function(p) {
    f <- coxph(as.formula(paste("Surv(time, event) ~ scale(x) +", CLIN_FORMULA)), cbind(sv, x = m[p, ]))
    s <- summary(f)
    data.table(pathway = p, HR = s$conf.int[1, 1], lower = s$conf.int[1, 3], upper = s$conf.int[1, 4],
               p = s$coefficients[1, "Pr(>|z|)"])
  }))
}), idcol = "layer")
pw_cox[, fdr := p.adjust(p, "BH"), by = layer]
write_tsv(pw_cox[order(p)], "survival_pathway_cox.tsv")
print(pw_cox[order(p)][1:12])

# ============================================================== figure 8 (part 1: subtypes and genes)
km <- function(var, title, pal) {
  f <- survfit2(as.formula(paste("Surv(time, event) ~", var)), data = sv)
  p <- lr[grouping == var, p]
  ggsurvfit(f, linewidth = 0.6) + add_censor_mark(size = 1.2) + add_risktable(risktable_stats = "n.risk", size = 2.2,
                                                                          theme = theme_risktable_default(axis.text.y.size = 6, plot.title.size = 6)) +
    scale_colour_manual(values = pal) + scale_x_continuous(breaks = seq(0, 60, 12)) +
    annotate("text", x = 60, y = 0.95, hjust = 1, size = 2.4, label = sprintf("log-rank p = %.2g", p)) +
    labs(title = title, x = "Months", y = "Overall survival") + theme_pub() + theme(legend.position = c(0.72, 0.75))
}
clpal <- function(n) setNames(RColorBrewer::brewer.pal(max(3, n), "Set2")[seq_len(n)], NULL)
kms <- list(
  km("purist_call", "PurIST (mRNA)", PAL$subtype),
  km("moffitt_rna", "Moffitt consensus (mRNA)", PAL$subtype),
  km("moffitt_protein", "Moffitt consensus (protein)", PAL$subtype),
  km("cluster_rna", "Unsupervised clusters (mRNA)", clpal(uniqueN(sv$cluster_rna))),
  km("cluster_protein", "Unsupervised clusters (protein)", clpal(uniqueN(sv$cluster_protein)))
)
pHR <- ggplot(hr_by_bin, aes(conc_bin, hr_spearman)) + geom_col(width = 0.65, fill = PAL$layer[["Protein"]]) +
  geom_text(aes(label = sprintf("%.2f\n(n=%d)", hr_spearman, genes)), vjust = -0.2, size = 2.2) +
  scale_y_continuous(limits = c(0, 1.05)) +
  labs(title = "Prognostic agreement vs\nRNA-protein concordance", x = "Gene-wise rho (mRNA vs protein)",
       y = "Spearman of log HR\n(mRNA vs protein)") + theme_pub()
kms_s <- kms
kms_s[[4]] <- kms_s[[4]] + labs(title = "Clusters (mRNA)")
kms_s[[5]] <- kms_s[[5]] + labs(title = "Clusters (protein)") +
  theme(legend.position = "right", legend.key.size = unit(2.5, "mm"), legend.text = element_text(size = 5.5))
figS <- wrap_plots(c(lapply(kms_s, function(p) free(wrap_elements(ggsurvfit_build(p)))), list(free(pHR))), ncol = 3) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 11))
save_fig(figS, "FigS1_survival_by_subtype", width = 210, height = 165)
pA <- kms[[2]]
pB <- kms[[3]]

fd <- cox_tab[model %in% c("purist", "moffitt_rna", "moffitt_protein", "score_rna", "score_protein", "concordance") &
                grepl("purist|moffitt|score|rho_samplewise", term)]
fd[, label := c(purist = "PurIST basal-like (mRNA)", moffitt_rna = "Moffitt basal-like (mRNA)",
                moffitt_protein = "Moffitt basal-like (protein)", score_rna = "Moffitt score, per SD (mRNA)",
                score_protein = "Moffitt score, per SD (protein)", concordance = "Per-tumour RNA-protein rho, per SD")[model]]
fd[, label := factor(label, levels = rev(label))]
pC <- ggplot(fd, aes(HR, label)) + geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, linewidth = 0.4) + geom_point(size = 1.6, colour = "#B2182B") +
  geom_text(aes(x = max(upper) * 1.05, label = sprintf("%.2f (%.2f-%.2f)  p=%.2g", HR, lower, upper, p)), hjust = 0, size = 2.1) +
  scale_x_log10() + coord_cartesian(xlim = c(min(fd$lower) * 0.9, max(fd$upper) * 4.5)) +
  labs(title = "Adjusted Cox models", subtitle = "age, stage, margin, grade", x = "Hazard ratio (log scale)", y = NULL) +
  theme_pub() + theme(plot.subtitle = element_text(size = 6.5))

pD <- ggplot(gw, aes(logHR_RNA, logHR_Protein)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey60") + geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey60") +
  geom_point(size = 0.25, alpha = 0.35, colour = "grey40") +
  geom_point(data = gw[p_RNA < 0.01 & p_Protein < 0.01], size = 0.9, colour = "#B2182B") +
  ggrepel::geom_text_repel(data = gw[p_RNA < 0.01 & p_Protein < 0.01][order(p_RNA * p_Protein)][1:12], aes(label = gene),
                           size = 1.9, max.overlaps = 20, segment.size = 0.2) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2, size = 2.3,
           label = sprintf("Spearman rho = %.2f\nred: p < 0.01 in both (n = %d)", hr_cor$estimate,
                           gw[p_RNA < 0.01 & p_Protein < 0.01, .N])) +
  labs(title = "Gene-level prognosis, stage-adjusted", x = "log HR per SD (mRNA)", y = "log HR per SD (protein)") + theme_pub()

# panels reused by Figure 8 (assembled in step 13 together with the prognostic models)
saveRDS(list(km_rna = ggsurvfit_build(pA), km_protein = ggsurvfit_build(pB), forest = pC, gene_hr = pD, hr_bin = pHR),
        file.path(ROOT, CFG$paths$interim, "fig8_panels.rds"))
