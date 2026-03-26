# Packages the pipeline needs but never attaches with library(), declared here so
# that renv (implicit snapshot) records them. This file is not sourced.
if (FALSE) {
  library(apeglm)    # DESeq2::lfcShrink(type = "apeglm") in step 09
  library(R.utils)   # data.table::fread() on .gz files
}
