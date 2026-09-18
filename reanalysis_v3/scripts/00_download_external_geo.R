#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

parse_cli <- function(values) {
  result <- list()
  for (value in values) {
    if (!grepl("^--[^=]+=", value)) stop("Arguments must have the form --name=value")
    fields <- strsplit(sub("^--", "", value), "=", fixed = TRUE)[[1L]]
    result[[fields[[1L]]]] <- paste(fields[-1L], collapse = "=")
  }
  result
}

as_flag <- function(value) tolower(value) %in% c("1", "true", "yes", "y")
cli <- parse_cli(commandArgs(trailingOnly = TRUE))
output_dir <- if (!is.null(cli[["output-dir"]])) cli[["output-dir"]] else "data/external_GEO"
force <- as_flag(if (!is.null(cli[["force"]])) cli[["force"]] else "false")

sources <- c(
  GSE14333 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14333/matrix/GSE14333_series_matrix.txt.gz",
  GSE17536 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE17nnn/GSE17536/matrix/GSE17536_series_matrix.txt.gz",
  GSE17537 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE17nnn/GSE17537/matrix/GSE17537_series_matrix.txt.gz"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
for (accession in names(sources)) {
  destination <- file.path(output_dir, paste0(accession, "_series_matrix.txt.gz"))
  if (file.exists(destination) && !force) {
    message("Keeping existing file: ", destination)
    next
  }
  temporary <- paste0(destination, ".download")
  on.exit(unlink(temporary), add = TRUE)
  status <- utils::download.file(sources[[accession]], temporary, mode = "wb", method = "libcurl")
  if (!identical(status, 0L) || !file.exists(temporary) || file.info(temporary)$size <= 0) {
    stop("Download failed for ", accession)
  }
  if (file.exists(destination)) unlink(destination)
  if (!file.rename(temporary, destination)) stop("Could not finalize ", destination)
  message("Downloaded: ", destination)
}

