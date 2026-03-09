# testthat checks for the shared R helpers.
#   Rscript -e "testthat::test_dir('tests/R')"

library(testthat)
source(file.path("..", "..", "R", "utils.R"), chdir = TRUE)
source(file.path("..", "..", "R", "survival_utils.R"), chdir = TRUE)

test_that("collapse_to_symbol keeps the most expressed row per symbol", {
  m <- matrix(c(1, 1, 5, 5, 2, 2), nrow = 3, byrow = TRUE, dimnames = list(c("E1", "E2", "E3"), c("s1", "s2")))
  genes <- data.table(gene_id = c("E1", "E2", "E3"), gene_name = c("A", "A", "B"))
  out <- collapse_to_symbol(m, genes)
  expect_equal(sort(rownames(out)), c("A", "B"))
  expect_equal(unname(out["A", ]), c(5, 5))
})

test_that("survival_table censors at the administrative limit and codes covariates", {
  an <- data.table(sample_key = c("C1_T", "C2_T", "C3_T"), case_submitter_id = c("C1", "C2", "C3"),
                   tissue = "Tumor", matched = TRUE)
  cases <- data.table(case_submitter_id = c("C1", "C2", "C3"), os_days = c(100, 4000, 500), os_event = c(1, 1, 0),
                      age_years = c(60, 70, 50), stage_group = c("II", "III", NA), residual_disease = c("R0", "R1", "RX"),
                      tumor_grade = c("G2", "G3", "G1"), sex = c("male", "female", "male"))
  d <- survival_table(an, cases, max_days = 1825)
  expect_equal(d[case_submitter_id == "C2", event], 0L)                # death after 5 years is censored
  expect_equal(d[case_submitter_id == "C2", round(time, 2)], round(1825 / 30.4375, 2))
  expect_equal(as.character(d[case_submitter_id == "C3", stage]), "unknown")
  expect_false("unknown" %in% levels(d$grade))                          # empty levels dropped
})
