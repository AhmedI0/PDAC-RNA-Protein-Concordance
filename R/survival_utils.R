# Survival helpers shared by steps 12 and 13.

suppressPackageStartupMessages(library(survival))

# Analysis-ready survival table for the matched tumours: overall survival with
# administrative censoring at `max_days`, and the clinical covariates coded for Cox models.
survival_table <- function(an, cases, max_days = CFG$survival$max_followup_days) {
  T <- an[tissue == "Tumor" & matched == TRUE, .(sample_key, case_submitter_id)]
  d <- merge(T, cases, by = "case_submitter_id")
  d <- d[!is.na(os_days) & !is.na(os_event)]
  d[, time := pmin(os_days, max_days) / 30.4375]
  d[, event := as.integer(os_event == 1 & os_days <= max_days)]
  d[, age10 := age_years / 10]
  d[, stage := factor(fcase(stage_group %in% c("I", "II"), "I-II", stage_group %in% c("III", "IV"), "III-IV",
                            default = "unknown"), levels = c("I-II", "III-IV", "unknown"))]
  d[, margin := factor(fcase(residual_disease == "R0", "R0", residual_disease %in% c("R1", "R2"), "R1-R2",
                             default = "RX/unknown"), levels = c("R0", "R1-R2", "RX/unknown"))]
  d[, grade := factor(fcase(tumor_grade %in% c("G1", "G2"), "G1-G2", tumor_grade %in% c("G3", "G4"), "G3-G4",
                            default = "unknown"), levels = c("G1-G2", "G3-G4", "unknown"))]
  d[, sex := factor(sex)]
  for (v in c("stage", "margin", "grade")) d[, (v) := droplevels(get(v))]   # no empty levels in Cox models
  d[]
}

CLIN_FORMULA <- "age10 + stage + margin + grade"

# Tidy one Cox fit
tidy_cox <- function(fit, model) {
  s <- summary(fit)
  ci <- s$conf.int
  data.table(model = model, term = rownames(ci), HR = ci[, 1], lower = ci[, 3], upper = ci[, 4],
             p = s$coefficients[, "Pr(>|z|)"], n = s$n, events = s$nevent,
             concordance = unname(s$concordance[1]))
}

# Harrell's C of a risk score (higher = worse)
cindex <- function(time, event, risk) {
  unname(survival::concordance(Surv(time, event) ~ risk, reverse = TRUE)$concordance)
}
