# Step 13 - Does the proteome add prognostic information to the transcriptome?
#
# Penalised Cox models (glmnet, elastic net alpha = 0.5) for overall survival,
# evaluated with 50 x repeated 5-fold cross-validation. Inside every training
# fold the penalty is tuned by an inner 5-fold CV (cv.glmnet, lambda.min), so
# the test fold never informs model fitting. Clinical covariates enter
# unpenalised (penalty.factor = 0) where included.
#
# Features: the genes measured in both layers (step 07); per layer the 2,000
# most variable across tumours (outcome-blind filter), z-scored. Protein values
# missing in a tumour are imputed with the gene's minimum - 0.5 (low abundance).
# Pathway models use the 50 Hallmark GSVA scores of each layer.
# Two low-dimensional comparators add the pre-specified Moffitt basal-minus-classical
# score (step 11, computed without outcome) to the clinical model, unpenalised.
#
# Performance: Harrell's C on each held-out fold, averaged per repeat.
# Signature: genes selected in >= 50% of the 250 training-fold fits (stability)
# and the coefficients of a final fit on all tumours.
# Risk groups: each tumour's mean out-of-fold risk rank across repeats, split
# at the median (no tumour is scored by a model that saw it).
#
#   Rscript scripts/R/13_prognostic_models.R

source("R/utils.R")
source("R/survival_utils.R")
suppressPackageStartupMessages({
  library(glmnet)
  library(parallel)
  library(patchwork)
  library(ggsurvfit)
  library(matrixStats)
})
SEED <- CFG$survival$seed
set.seed(SEED)
N_REP <- CFG$survival$cv_repeats
K <- CFG$survival$cv_folds
ALPHA <- 0.5
N_FEAT <- 2000

an <- fread(path_proc("analysis_samples.tsv"))
sv <- survival_table(an, read_cases())
T <- sv$sample_key
conc <- fread(path_tab("concordance_gene.tsv"))

zs <- function(m) (m - rowMeans(m)) / rowSds(m)
rna <- read_matrix("rna_vst.tsv.gz")[conc$gene, T]
prot <- read_matrix("protein_norm.tsv.gz")[conc$gene, T]
prot <- prot[rowMeans(!is.na(prot)) >= 0.8, ]
for (i in which(rowSums(is.na(prot)) > 0)) { v <- prot[i, ]; v[is.na(v)] <- min(v, na.rm = TRUE) - 0.5; prot[i, ] <- v }
top <- function(m) m[order(rowVars(m), decreasing = TRUE)[1:min(N_FEAT, nrow(m))], ]
X_rna <- t(zs(top(rna)));  colnames(X_rna) <- paste0("rna_", colnames(X_rna))
X_pro <- t(zs(top(prot))); colnames(X_pro) <- paste0("prot_", colnames(X_pro))
g_r <- t(zs(read_matrix("gsva_hallmark_rna.tsv.gz")[, T]));     colnames(g_r) <- paste0("gsvaRNA_", colnames(g_r))
g_p <- t(zs(read_matrix("gsva_hallmark_protein.tsv.gz")[, T])); colnames(g_p) <- paste0("gsvaProt_", colnames(g_p))
X_cli <- model.matrix(as.formula(paste("~", CLIN_FORMULA)), sv)[, -1]
y <- cbind(time = sv$time, status = sv$event)
# pre-specified, outcome-blind Moffitt basal-minus-classical scores (step 11), one covariate each
sub <- fread(path_tab("subtype_calls.tsv"))
i <- match(T, sub$sample_key)
X_mof_r <- cbind(X_cli, moffitt_score = as.numeric(scale(sub$score_rna[i])))
X_mof_p <- cbind(X_cli, moffitt_score = as.numeric(scale(sub$score_protein[i])))

MODELS <- list(
  `Clinical`               = list(X = X_cli, pen = FALSE),
  `Clinical + Moffitt score (mRNA)`    = list(X = X_mof_r, pen = FALSE),
  `Clinical + Moffitt score (protein)` = list(X = X_mof_p, pen = FALSE),
  `mRNA genes`             = list(X = X_rna, pen = TRUE),
  `Protein genes`          = list(X = X_pro, pen = TRUE),
  `mRNA + protein genes`   = list(X = cbind(X_rna, X_pro), pen = TRUE),
  `Clinical + mRNA`        = list(X = cbind(X_cli, X_rna), pen = TRUE, n_cli = ncol(X_cli)),
  `Clinical + protein`     = list(X = cbind(X_cli, X_pro), pen = TRUE, n_cli = ncol(X_cli)),
  `Clinical + mRNA pathways`    = list(X = cbind(X_cli, g_r), pen = TRUE, n_cli = ncol(X_cli)),
  `Clinical + protein pathways` = list(X = cbind(X_cli, g_p), pen = TRUE, n_cli = ncol(X_cli))
)

# fold assignment stratified by event status, fixed per repeat for all models (paired comparison)
fold_ids <- lapply(seq_len(N_REP), function(r) {
  set.seed(SEED + r)
  f <- integer(nrow(sv))
  for (e in 0:1) { idx <- which(sv$event == e); f[idx] <- sample(rep_len(seq_len(K), length(idx))) }
  f
})

fit_predict <- function(model, train, test, inner_seed) {
  m <- MODELS[[model]]
  if (!m$pen) {
    d <- data.frame(time = y[train, 1], status = y[train, 2], m$X[train, , drop = FALSE])
    f <- survival::coxph(survival::Surv(time, status) ~ ., data = d)
    return(list(lp = drop(m$X[test, , drop = FALSE] %*% coef(f)), sel = character(0)))
  }
  pf <- rep(1, ncol(m$X)); if (!is.null(m$n_cli)) pf[seq_len(m$n_cli)] <- 0
  set.seed(inner_seed)
  cv <- glmnet::cv.glmnet(m$X[train, ], survival::Surv(y[train, 1], y[train, 2]), family = "cox", alpha = ALPHA,
                          nfolds = 5, penalty.factor = pf)
  b <- as.matrix(coef(cv, s = "lambda.min"))[, 1]
  list(lp = drop(m$X[test, , drop = FALSE] %*% b), sel = names(b)[b != 0])
}

jobs <- expand.grid(model = names(MODELS), rep = seq_len(N_REP), stringsAsFactors = FALSE)
cl <- makeCluster(min(10, detectCores() - 2))
clusterExport(cl, c("MODELS", "y", "fold_ids", "fit_predict", "jobs", "K", "ALPHA", "SEED"))
log_step(sprintf("CV: %d models x %d repeats x %d folds on %d workers", length(MODELS), N_REP, K, length(cl)))
res <- parLapply(cl, seq_len(nrow(jobs)), function(j) {
  model <- jobs$model[j]; r <- jobs$rep[j]; f <- fold_ids[[r]]
  oof <- rep(NA_real_, nrow(y)); cs <- numeric(K); sel <- list()
  for (k in seq_len(K)) {
    tr <- which(f != k); te <- which(f == k)
    out <- fit_predict(model, tr, te, SEED + 1000 * r + k)
    oof[te] <- out$lp
    cs[k] <- survival::concordance(survival::Surv(y[te, 1], y[te, 2]) ~ out$lp, reverse = TRUE)$concordance
    sel[[k]] <- out$sel
  }
  list(model = model, rep = r, c_fold_mean = mean(cs), oof = oof, sel = unlist(sel))
})
stopCluster(cl)
log_step("CV done")

perf <- rbindlist(lapply(res, function(x) data.table(model = x$model, rep = x$rep, cindex = x$c_fold_mean)))
perf[, model := factor(model, levels = names(MODELS))]
clin <- perf[model == "Clinical", .(rep, c_clin = cindex)]
perf <- merge(perf, clin, by = "rep")
summ <- perf[, .(mean_cindex = mean(cindex), q025 = quantile(cindex, 0.025), q975 = quantile(cindex, 0.975),
                 mean_gain_vs_clinical = mean(cindex - c_clin),
                 frac_repeats_better_than_clinical = mean(cindex > c_clin)), by = model][order(model)]
write_tsv(perf[, .(model, rep, cindex)], "prognostic_cv_cindex_by_repeat.tsv")
write_tsv(summ, "prognostic_cv_summary.tsv")
print(summ)

# ---------------------------------------------------------------- signature stability and final fits
stab <- rbindlist(lapply(c("mRNA genes", "Protein genes", "mRNA + protein genes"), function(m) {
  s <- unlist(lapply(res[sapply(res, `[[`, "model") == m], `[[`, "sel"))
  if (!length(s)) return(NULL)
  data.table(model = m, feature = names(table(s)), frequency = as.integer(table(s)) / (N_REP * K))
}))
final <- rbindlist(lapply(c("mRNA genes", "Protein genes"), function(m) {
  X <- MODELS[[m]]$X
  set.seed(SEED)
  cv <- cv.glmnet(X, Surv(y[, 1], y[, 2]), family = "cox", alpha = ALPHA, nfolds = 5)
  b <- as.matrix(coef(cv, s = "lambda.min"))[, 1]
  data.table(model = m, feature = names(b)[b != 0], coef_full_data = b[b != 0])
}))
sig <- merge(stab, final, by = c("model", "feature"), all = TRUE)
sig[is.na(frequency), frequency := 0]
sig[, layer := fifelse(grepl("^rna_", feature), "mRNA", "Protein")]
sig[, gene := sub("^(rna|prot)_", "", feature)]
sig <- sig[frequency >= 0.5 | !is.na(coef_full_data)][order(model, -frequency)]
write_tsv(sig, "prognostic_signature_genes.tsv")
log_step(sprintf("stable signature genes (selected in >=50%% of fits): mRNA %d, protein %d",
                 sig[model == "mRNA genes" & frequency >= 0.5, .N], sig[model == "Protein genes" & frequency >= 0.5, .N]))

# ---------------------------------------------------------------- cross-validated risk groups
risk_group <- function(m) {
  r <- sapply(res[sapply(res, `[[`, "model") == m], function(x) rank(x$oof) / length(x$oof))
  g <- factor(ifelse(rowMeans(r) > median(rowMeans(r)), "High risk", "Low risk"), levels = c("Low risk", "High risk"))
  g
}
sv[, risk_protein := risk_group("Protein genes")]
sv[, risk_rna := risk_group("mRNA genes")]
write_tsv(sv[, .(sample_key, case_submitter_id, time, event, risk_rna, risk_protein)], "prognostic_cv_risk_groups.tsv")
lr_p <- function(v) { s <- survdiff(as.formula(paste("Surv(time, event) ~", v)), sv); pchisq(s$chisq, 1, lower.tail = FALSE) }
cox_risk <- rbind(tidy_cox(coxph(as.formula(paste("Surv(time, event) ~ risk_protein +", CLIN_FORMULA)), sv), "cv risk, protein"),
                  tidy_cox(coxph(as.formula(paste("Surv(time, event) ~ risk_rna +", CLIN_FORMULA)), sv), "cv risk, mRNA"))
write_tsv(cox_risk, "prognostic_cv_risk_cox.tsv")
print(cox_risk[grepl("risk", term)])

# ================================================================ figure 8
P <- readRDS(file.path(ROOT, CFG$paths$interim, "fig8_panels.rds"))
mcol <- c(Clinical = "grey55", `Clinical + Moffitt score (mRNA)` = "#F4A582", `Clinical + Moffitt score (protein)` = "#D6604D",
          `mRNA genes` = PAL$layer[["RNA"]], `Protein genes` = PAL$layer[["Protein"]],
          `mRNA + protein genes` = "#5E3C99", `Clinical + mRNA` = "#A6DBA0", `Clinical + protein` = "#C2A5CF",
          `Clinical + mRNA pathways` = "#5AAE61", `Clinical + protein pathways` = "#9970AB")
pE <- ggplot(perf, aes(cindex, reorder(model, cindex, median), fill = model)) +
  geom_vline(xintercept = 0.5, linetype = 3, linewidth = 0.3) +
  geom_vline(xintercept = summ[model == "Clinical", mean_cindex], linetype = 2, linewidth = 0.3, colour = "grey40") +
  geom_boxplot(outlier.size = 0.3, width = 0.6, linewidth = 0.3) + scale_fill_manual(values = mcol, guide = "none") +
  labs(title = sprintf("Cross-validated C-index (%d x %d-fold)", N_REP, K), x = "Harrell's C (held-out folds)", y = NULL) +
  theme_pub()
f_r <- survfit2(Surv(time, event) ~ risk_protein, data = sv)
pF <- ggsurvfit(f_r, linewidth = 0.6) + add_censor_mark(size = 1.2) +
  add_risktable(risktable_stats = "n.risk", size = 2.2, theme = theme_risktable_default(axis.text.y.size = 6, plot.title.size = 6)) +
  scale_colour_manual(values = c(`Low risk` = "#4393C3", `High risk` = "#B2182B")) + scale_x_continuous(breaks = seq(0, 60, 12)) +
  annotate("text", x = 60, y = 0.95, hjust = 1, size = 2.4, label = sprintf("log-rank p = %.2g", lr_p("risk_protein"))) +
  labs(title = "Protein signature, out-of-fold risk", x = "Months", y = "Overall survival") + theme_pub() +
  theme(legend.position = c(0.72, 0.75))

fig8 <- wrap_plots(list(P$km_rna, P$km_protein, P$forest, P$gene_hr, P$hr_bin, pE, wrap_elements(ggsurvfit_build(pF))),
                   design = "AABBCCC\nDDDEEFF\nGGG####") + plot_annotation(tag_levels = "A")
save_fig(fig8, "Fig8_survival", width = 210, height = 215)
