# Build the project R library with renv.
#
# R 4.5.2, Bioconductor 3.22, and CRAN packages frozen at the Posit Package
# Manager snapshot of 2025-12-15 (Windows binaries). Run once from the project
# root:
#
#   Rscript envs/setup_r_env.R
#
# Afterwards `renv::restore()` reproduces the same library from renv.lock.

cran_snapshot <- "https://packagemanager.posit.co/cran/2025-12-15"
bioc_version  <- "3.22"

stopifnot(getRversion() >= "4.5.0", getRversion() < "4.6.0")
options(repos = c(CRAN = cran_snapshot), timeout = 600)

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", lib = .libPaths()[1])
}
if (!file.exists("renv/activate.R")) {
  renv::init(bare = TRUE, restart = FALSE, settings = list(
    bioconductor.version = bioc_version,
    snapshot.type = "implicit"
  ))
}
source("renv/activate.R")
options(repos = c(CRAN = cran_snapshot))
renv::settings$bioconductor.version(bioc_version)

cran_pkgs <- c(
  "BiocManager", "data.table", "yaml", "jsonlite", "DBI", "RSQLite",
  "ggplot2", "ggrepel", "patchwork", "scales", "svglite", "ragg",
  "survival", "ggsurvfit", "glmnet", "matrixStats", "tidyestimate",
  "aricode", "msigdbr", "testthat", "R.utils", "RColorBrewer"
)
bioc_pkgs <- c(
  "DESeq2", "limma", "apeglm", "fgsea", "GSVA",
  "ConsensusClusterPlus", "AnnotationDbi",
  "org.Hs.eg.db", "SummarizedExperiment", "BiocParallel"
)

renv::install(c(cran_pkgs, paste0("bioc::", bioc_pkgs)), prompt = FALSE)
renv::snapshot(packages = c("renv", cran_pkgs, bioc_pkgs), prompt = FALSE)

cat("\nInstalled:\n")
print(sapply(c(cran_pkgs, bioc_pkgs), function(p) as.character(packageVersion(p))))
