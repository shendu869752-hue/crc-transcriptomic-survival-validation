#!/usr/bin/env Rscript

cran_packages <- c(
  "data.table", "digest", "future", "future.apply", "ggplot2", "glmnet",
  "metafor", "survival", "timeROC"
)

missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioconductor_packages <- c(
  "AnnotationDbi", "Biobase", "BiocGenerics", "Biostrings",
  "hgu133plus2.db", "IRanges", "org.Hs.eg.db", "S4Vectors"
)
missing_bioc <- bioconductor_packages[
  !vapply(bioconductor_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

cat("Dependencies are installed. Compare package versions with sessionInfo.txt.\n")

