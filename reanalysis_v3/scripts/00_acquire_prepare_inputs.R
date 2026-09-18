#!/usr/bin/env Rscript

# Reconstruct or audit the two gene-level expression matrices consumed by v3.
#
# Default behavior is read-only with respect to the source project: it rebuilds
# both matrices in memory from already-downloaded public raw files, compares
# them with the current derived inputs, and writes only small audit tables under
# reanalysis_v3/provenance/. Large downloads occur only with --download=true.
#
# Historical GSE39582 reconstruction is intentionally reproduced as implemented
# in the legacy code, not relabelled as a better method. It used the first
# AnnotationDbi SYMBOL edge per probe and then the maximum probe value within
# each gene separately for every sample. The A-level alternatives are exposed
# separately as unique_highest_mean and unique_mean. The fixed, annotation-only
# unique_mean matrix is the primary A-level input. A cohort-wide
# unique_highest_mean matrix is a transductive reference only; its planned
# sensitivity analysis must select the probe inside each training partition.

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(hgu133plus2.db)
})

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Package 'digest' is required for SHA-256 manifests", call. = FALSE)
}

SCRIPT_VERSION <- "2026-09-17.2"

PUBLIC_SOURCES <- data.frame(
  source_id = c(
    "GSE39582_series_matrix",
    "TCGA_COAD_HiSeqV2",
    "TCGA_PANCAN_survival",
    "TCGA_COAD_clinical"
  ),
  dataset_id = c(
    "GSE39582; GPL570",
    "TCGA.COAD.sampleMap/HiSeqV2",
    "Survival_SupplementalTable_S1_20171025_xena_sp",
    "TCGA.COAD.sampleMap/COAD_clinicalMatrix"
  ),
  version = c(
    "GEO accession; series metadata last update 2021-06-18; file is mutable",
    "2017-10-13",
    "2018-09-13",
    "2019-12-06"
  ),
  data_url = c(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE39nnn/GSE39582/matrix/GSE39582_series_matrix.txt.gz",
    "https://tcga.xenahubs.net/download/TCGA.COAD.sampleMap/HiSeqV2.gz",
    "https://pancanatlas.xenahubs.net/download/Survival_SupplementalTable_S1_20171025_xena_sp",
    "https://tcga.xenahubs.net/download/TCGA.COAD.sampleMap/COAD_clinicalMatrix"
  ),
  metadata_url = c(
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39582",
    "https://tcga.xenahubs.net/download/TCGA.COAD.sampleMap/HiSeqV2.json",
    "https://pancanatlas.xenahubs.net/download/Survival_SupplementalTable_S1_20171025_xena_sp.json",
    "https://tcga.xenahubs.net/download/TCGA.COAD.sampleMap/COAD_clinicalMatrix.json"
  ),
  raw_filename = c(
    "GSE39582_series_matrix.txt.gz",
    "COAD_HiSeqV2.gz",
    "TCGA_PANCAN_survival.tsv",
    "COAD_clinicalMatrix.tsv"
  ),
  stringsAsFactors = FALSE
)

parse_cli <- function(x) {
  out <- list()
  for (item in x) {
    if (!grepl("^--[^=]+=", item)) {
      stop("Arguments must have the form --name=value: ", item, call. = FALSE)
    }
    kv <- strsplit(sub("^--", "", item), "=", fixed = TRUE)[[1]]
    out[[kv[1]]] <- paste(kv[-1], collapse = "=")
  }
  out
}

as_flag <- function(x, name) {
  z <- tolower(trimws(x))
  if (z %in% c("true", "1", "yes", "y")) return(TRUE)
  if (z %in% c("false", "0", "no", "n")) return(FALSE)
  stop("Invalid logical value for --", name, ": ", x, call. = FALSE)
}

cli <- parse_cli(commandArgs(trailingOnly = TRUE))
arg <- function(name, default) if (!is.null(cli[[name]])) cli[[name]] else default

cmd_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_all, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "00_acquire_prepare_inputs.R"
script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
v3_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)

action <- arg("action", "verify-current")
valid_actions <- c(
  "verify-current", "write-current", "audit-a-level", "write-a-level",
  "metadata-only"
)
if (!action %in% valid_actions) {
  stop("--action must be one of: ", paste(valid_actions, collapse = ", "), call. = FALSE)
}

source_root <- arg(
  "source-root",
  Sys.getenv("CRC_SOURCE_ROOT", normalizePath(file.path(v3_root, ".."), winslash = "/", mustWork = FALSE))
)
derived_dir <- arg("derived-dir", file.path(source_root, "data"))
current_results_dir <- arg("current-results-dir", file.path(source_root, "results"))
raw_dir <- arg("raw-dir", Sys.getenv("CRC_RAW_DIR", file.path(v3_root, "raw_inputs")))
output_dir <- arg("output-dir", file.path(v3_root, "prepared_inputs"))
provenance_dir <- arg("provenance-dir", file.path(v3_root, "provenance"))
download_requested <- as_flag(arg("download", "false"), "download")
force <- as_flag(arg("force", "false"), "force")

dir.create(provenance_dir, recursive = TRUE, showWarnings = FALSE)

normalize_existing <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path, call. = FALSE)
  # Keep an ASCII junction/alias spelling when supplied. normalizePath() can
  # resolve it to a non-ASCII target that older Windows locale handling in R
  # cannot pass back to file.info()/digest reliably.
  gsub("\\\\", "/", path)
}

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

file_record <- function(role, path) {
  path <- normalize_existing(path)
  data.frame(
    role = role,
    path = path,
    bytes = unname(file.info(path)$size),
    md5 = unname(tools::md5sum(path)),
    sha256 = sha256_file(path),
    stringsAsFactors = FALSE
  )
}

download_one <- function(url, dest) {
  if (file.exists(dest) && !force) {
    message("Keeping existing raw file: ", dest)
    return(invisible(dest))
  }
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(dest, ".part")
  if (file.exists(tmp)) unlink(tmp)
  message("Downloading: ", url)
  status <- utils::download.file(url, tmp, mode = "wb", method = "libcurl", quiet = FALSE)
  if (!identical(status, 0L) || !file.exists(tmp) || file.info(tmp)$size <= 0) {
    stop("Download failed: ", url, call. = FALSE)
  }
  if (file.exists(dest)) unlink(dest)
  if (!file.rename(tmp, dest)) stop("Could not finalize download: ", dest, call. = FALSE)
  invisible(dest)
}

download_public_sources <- function() {
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_len(nrow(PUBLIC_SOURCES))) {
    download_one(PUBLIC_SOURCES$data_url[i], file.path(raw_dir, PUBLIC_SOURCES$raw_filename[i]))
    if (nzchar(PUBLIC_SOURCES$metadata_url[i]) && grepl("\\.json$", PUBLIC_SOURCES$metadata_url[i])) {
      download_one(
        PUBLIC_SOURCES$metadata_url[i],
        file.path(raw_dir, paste0(PUBLIC_SOURCES$source_id[i], ".metadata.json"))
      )
    }
  }
}

if (download_requested) download_public_sources()

fwrite(PUBLIC_SOURCES, file.path(provenance_dir, "source_manifest_runtime.csv"))

find_input <- function(dir, candidates, label) {
  paths <- file.path(dir, candidates)
  hit <- paths[file.exists(paths)]
  if (!length(hit)) {
    stop(
      "Missing ", label, " in ", dir, ". Tried: ",
      paste(candidates, collapse = ", "),
      ". Use --download=true or provide --raw-dir.",
      call. = FALSE
    )
  }
  hit[1]
}

open_text <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
}

scan_geo_header <- function(path) {
  con <- open_text(path)
  on.exit(close(con), add = TRUE)
  line_no <- 0L
  begin_line <- NA_integer_
  row_count <- NA_integer_
  series <- platform <- processing <- series_last_update <- NA_character_
  sample_accessions <- character()
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) break
    line_no <- line_no + 1L
    if (identical(line, "!series_matrix_table_begin")) {
      begin_line <- line_no
      break
    }
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (!length(fields) || is.na(fields[1]) || !nzchar(fields[1])) next
    key <- fields[1]
    values <- if (length(fields) > 1L) gsub('^"|"$', "", fields[-1]) else character()
    if (key == "!Series_geo_accession") series <- values[1]
    if (key == "!Series_platform_id") platform <- values[1]
    if (key == "!Series_last_update_date") series_last_update <- values[1]
    if (key == "!Sample_data_processing") processing <- unique(values)[1]
    if (key == "!Sample_data_row_count") row_count <- as.integer(unique(values)[1])
    if (key == "!Sample_geo_accession") sample_accessions <- values
  }
  if (!is.finite(begin_line) || !is.finite(row_count)) {
    stop("Could not locate a complete GEO series-matrix header in ", path, call. = FALSE)
  }
  list(
    begin_line = begin_line,
    row_count = row_count,
    series = series,
    platform = platform,
    processing = processing,
    series_last_update = series_last_update,
    sample_accessions = sample_accessions
  )
}

read_geo_probe_matrix <- function(path) {
  meta <- scan_geo_header(path)
  message("Reading ", meta$series, " probe matrix: ", meta$row_count, " rows")
  if (!grepl("\\.gz$", path, ignore.case = TRUE)) {
    dat <- fread(
      path,
      sep = "\t",
      skip = meta$begin_line,
      nrows = meta$row_count,
      header = TRUE,
      check.names = FALSE,
      quote = '"',
      data.table = FALSE
    )
  } else {
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    dat <- read.delim(
      con,
      sep = "\t",
      skip = meta$begin_line,
      nrows = meta$row_count,
      header = TRUE,
      check.names = FALSE,
      quote = '"',
      comment.char = "",
      stringsAsFactors = FALSE
    )
  }
  if (!identical(names(dat)[1], "ID_REF")) stop("Expected ID_REF as the first GEO column")
  probe_ids <- as.character(dat[[1]])
  if (anyDuplicated(probe_ids)) stop("Duplicated PROBEID values in GSE39582 series matrix")
  expr <- as.matrix(dat[, -1, drop = FALSE])
  storage.mode(expr) <- "double"
  rownames(expr) <- probe_ids
  if (length(meta$sample_accessions) && !identical(colnames(expr), meta$sample_accessions)) {
    stop("GEO matrix columns do not equal !Sample_geo_accession order")
  }
  if (!identical(meta$series, "GSE39582") || !identical(meta$platform, "GPL570")) {
    stop("Unexpected GEO identity: ", meta$series, " / ", meta$platform)
  }
  list(expr = expr, meta = meta)
}

annotation_maps <- function(probe_ids) {
  edges <- AnnotationDbi::select(
    hgu133plus2.db,
    keys = probe_ids,
    columns = "SYMBOL",
    keytype = "PROBEID"
  )
  edges <- edges[!is.na(edges$SYMBOL) & nzchar(edges$SYMBOL), c("PROBEID", "SYMBOL")]
  edges <- edges[edges$PROBEID %in% probe_ids, , drop = FALSE]

  # Exact historical rule: first AnnotationDbi row for each observed PROBEID.
  legacy <- edges[!duplicated(edges$PROBEID), , drop = FALSE]

  # A-level annotation rule: retain only probes with one distinct nonempty
  # SYMBOL. This annotation-only eligibility step is fixed before modelling;
  # expression-dependent choice among eligible probes is handled separately.
  edge_dt <- unique(as.data.table(edges))
  edge_split <- split(edge_dt$SYMBOL, edge_dt$PROBEID)
  edge_split <- lapply(edge_split, function(x) sort(unique(x)))
  symbols <- vapply(probe_ids, function(id) {
    z <- edge_split[[id]]
    if (is.null(z)) "" else paste(z, collapse = ";")
  }, character(1))
  n_symbols <- lengths(edge_split)[probe_ids]
  n_symbols[is.na(n_symbols)] <- 0L
  audit <- data.frame(
    PROBEID = probe_ids,
    symbols = unname(symbols),
    n_symbols = as.integer(n_symbols),
    eligible = as.integer(n_symbols) == 1L,
    stringsAsFactors = FALSE
  )
  audit$SYMBOL <- ifelse(audit$eligible, audit$symbols, NA_character_)
  one <- audit[audit$eligible, c("PROBEID", "n_symbols", "SYMBOL"), drop = FALSE]

  list(edges = edges, legacy = legacy, one = one, audit = audit)
}

aggregate_geo <- function(expr, method, maps = NULL) {
  stopifnot(method %in% c("legacy_per_sample_max", "unique_highest_mean", "unique_mean"))
  if (is.null(maps)) maps <- annotation_maps(rownames(expr))
  samples <- colnames(expr)

  if (method == "legacy_per_sample_max") {
    symbols <- maps$legacy$SYMBOL[match(rownames(expr), maps$legacy$PROBEID)]
    keep <- !is.na(symbols) & nzchar(symbols)
    dt <- as.data.table(expr[keep, , drop = FALSE])
    dt[, SYMBOL := symbols[keep]]
    out <- dt[, lapply(.SD, max, na.rm = TRUE), by = SYMBOL, .SDcols = samples]
    setorder(out, SYMBOL)
    mat <- as.matrix(out[, ..samples])
    storage.mode(mat) <- "double"
    rownames(mat) <- out$SYMBOL
    mapping <- data.frame(
      method = method,
      note = paste(
        "Historical reproduction only: first SYMBOL edge per PROBEID,",
        "then per-sample maximum across probes; probe identity may vary by sample"
      ),
      stringsAsFactors = FALSE
    )
    return(list(expr = mat, mapping = mapping, maps = maps))
  }

  eligible <- maps$one
  eligible$idx <- match(eligible$PROBEID, rownames(expr))
  eligible$probe_mean <- rowMeans(expr[eligible$idx, , drop = FALSE], na.rm = TRUE)
  gene_probe_counts <- table(eligible$SYMBOL)
  eligible$n_probes_gene <- as.integer(gene_probe_counts[eligible$SYMBOL])
  eligible <- eligible[order(eligible$SYMBOL, -eligible$probe_mean, eligible$PROBEID), , drop = FALSE]

  if (method == "unique_highest_mean") {
    chosen <- eligible[!duplicated(eligible$SYMBOL), , drop = FALSE]
    mat <- expr[chosen$idx, , drop = FALSE]
    rownames(mat) <- chosen$SYMBOL
    mapping <- chosen[, c(
      "SYMBOL", "PROBEID", "probe_mean", "n_symbols", "n_probes_gene"
    )]
    mapping$method <- method
    return(list(expr = mat, mapping = mapping, maps = maps))
  }

  dt <- as.data.table(expr[eligible$idx, , drop = FALSE])
  dt[, SYMBOL := eligible$SYMBOL]
  out <- dt[, lapply(.SD, mean, na.rm = TRUE), by = SYMBOL, .SDcols = samples]
  setorder(out, SYMBOL)
  mat <- as.matrix(out[, ..samples])
  storage.mode(mat) <- "double"
  rownames(mat) <- out$SYMBOL
  mapping <- eligible[, c(
    "SYMBOL", "PROBEID", "probe_mean", "n_symbols", "n_probes_gene"
  )]
  mapping$method <- method
  list(expr = mat, mapping = mapping, maps = maps)
}

read_geo_clinical <- function(path, expected_samples = NULL) {
  con <- open_text(path)
  on.exit(close(con), add = TRUE)
  sample_accessions <- character()
  fields <- list()
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line) || identical(line, "!series_matrix_table_begin")) break
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (!length(parts)) next
    if (identical(parts[1], "!Sample_geo_accession")) {
      sample_accessions <- gsub('^"|"$', "", parts[-1])
    }
    if (identical(parts[1], "!Sample_characteristics_ch1")) {
      values <- gsub('^"|"$', "", parts[-1])
      if (!length(values)) next
      key <- trimws(sub(":.*$", "", values[1]))
      if (key %in% names(fields)) stop("Duplicated GEO characteristic: ", key)
      fields[[key]] <- trimws(sub("^[^:]*:", "", values))
    }
  }
  required <- c("os.event", "os.delay (months)", "tnm.stage", "age.at.diagnosis (year)")
  if (!length(sample_accessions) || !all(required %in% names(fields))) {
    stop(
      "GSE39582 header lacks sample IDs or required clinical fields: ",
      paste(setdiff(required, names(fields)), collapse = ", ")
    )
  }
  bad_lengths <- required[lengths(fields[required]) != length(sample_accessions)]
  if (length(bad_lengths)) {
    stop("GEO clinical field length mismatch: ", paste(bad_lengths, collapse = ", "))
  }
  if (!is.null(expected_samples) && !identical(sample_accessions, expected_samples)) {
    stop("GEO clinical sample order does not equal expression sample order")
  }
  data.frame(
    sample = sample_accessions,
    os_event = fields[["os.event"]],
    os_delay_months = fields[["os.delay (months)"]],
    tnm_stage = fields[["tnm.stage"]],
    age = fields[["age.at.diagnosis (year)"]],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

read_xena_matrix <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  dat <- read.delim(
    con,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  )
  genes <- as.character(dat[[1]])
  expr <- as.matrix(dat[, -1, drop = FALSE])
  storage.mode(expr) <- "double"
  rownames(expr) <- genes
  if (anyDuplicated(rownames(expr))) {
    dt <- as.data.table(expr)
    dt[, gene := genes]
    sample_cols <- colnames(expr)
    agg <- dt[, lapply(.SD, mean), by = gene, .SDcols = sample_cols]
    setorder(agg, gene)
    expr <- as.matrix(agg[, ..sample_cols])
    storage.mode(expr) <- "double"
    rownames(expr) <- agg$gene
  }
  expr
}

read_survival <- function(path) {
  dat <- fread(path, sep = "\t", header = TRUE, data.table = FALSE, na.strings = c("NA"))
  required <- c("sample", "OS", "OS.time")
  if (!all(required %in% names(dat))) stop("Survival file lacks: ", paste(setdiff(required, names(dat)), collapse = ", "))
  dat$sample <- as.character(dat$sample)
  dat$OS <- suppressWarnings(as.numeric(dat$OS))
  dat$OS.time <- suppressWarnings(as.numeric(dat$OS.time))
  dat[is.finite(dat$OS.time) & dat$OS.time > 0 & dat$OS %in% 0:1, , drop = FALSE]
}

prepare_tcga_current <- function(expr_path, survival_path) {
  expr <- read_xena_matrix(expr_path)
  survival <- read_survival(survival_path)
  is_tumor <- vapply(strsplit(colnames(expr), "-", fixed = TRUE), function(parts) {
    length(parts) >= 4L && identical(parts[4], "01")
  }, logical(1))
  tumor_samples <- colnames(expr)[is_tumor]
  common <- intersect(tumor_samples, survival$sample)
  if (!length(common)) stop("No TCGA type-01 samples overlap eligible survival records")
  list(
    expr = expr[, common, drop = FALSE],
    survival = survival[match(common, survival$sample), , drop = FALSE],
    common = common
  )
}

prepare_tcga_clinical <- function(path, common_samples) {
  dat <- read.delim(
    path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE,
    quote = "", comment.char = ""
  )
  keep <- c(
    "sampleID", "pathologic_stage", "pathologic_T", "pathologic_N",
    "pathologic_M", "gender", "age_at_initial_pathologic_diagnosis"
  )
  if (!all(keep %in% names(dat))) {
    stop("TCGA clinical matrix lacks: ", paste(setdiff(keep, names(dat)), collapse = ", "))
  }
  dat <- dat[, keep, drop = FALSE]
  if (anyDuplicated(dat$sampleID)) stop("Duplicated sampleID in TCGA clinical matrix")
  idx <- match(common_samples, dat$sampleID)
  if (anyNA(idx)) {
    stop("TCGA clinical matrix is missing matched expression samples: ",
         paste(common_samples[is.na(idx)], collapse = ", "))
  }
  dat <- dat[idx, , drop = FALSE]
  rownames(dat) <- NULL
  if (!identical(as.character(dat$sampleID), common_samples)) {
    stop("TCGA clinical projection order mismatch")
  }
  dat
}

read_derived_expr <- function(path) {
  dat <- fread(path, check.names = FALSE, data.table = FALSE)
  if (ncol(dat) < 2L) stop("Expression CSV has fewer than two columns: ", path)
  genes <- as.character(dat[[1]])
  mat <- as.matrix(dat[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- genes
  mat
}

compare_matrices <- function(label, rebuilt, current, current_path) {
  same_dim <- identical(dim(rebuilt), dim(current))
  same_genes <- identical(rownames(rebuilt), rownames(current))
  same_samples <- identical(colnames(rebuilt), colnames(current))
  mismatches <- NA_real_
  max_abs <- NA_real_
  if (same_dim) {
    mismatches <- 0
    max_abs <- 0
    chunk <- 500L
    for (lo in seq.int(1L, nrow(rebuilt), by = chunk)) {
      hi <- min(nrow(rebuilt), lo + chunk - 1L)
      a <- rebuilt[lo:hi, , drop = FALSE]
      b <- current[lo:hi, , drop = FALSE]
      bad_na <- xor(is.na(a), is.na(b))
      d <- abs(a - b)
      mismatches <- mismatches + sum(bad_na | (!is.na(d) & d != 0))
      finite_d <- d[is.finite(d)]
      if (length(finite_d)) max_abs <- max(max_abs, max(finite_d))
    }
  }
  hashes <- file_record(label, current_path)
  data.frame(
    cohort = label,
    rebuilt_genes = nrow(rebuilt),
    rebuilt_samples = ncol(rebuilt),
    current_genes = nrow(current),
    current_samples = ncol(current),
    dimensions_exact = same_dim,
    gene_order_exact = same_genes,
    sample_order_exact = same_samples,
    numeric_mismatches = mismatches,
    max_abs_difference = max_abs,
    semantic_exact = isTRUE(same_dim && same_genes && same_samples && mismatches == 0),
    current_md5 = hashes$md5,
    current_sha256 = hashes$sha256,
    stringsAsFactors = FALSE
  )
}

canonical_table <- function(x) {
  out <- lapply(x, function(v) {
    z <- as.character(v)
    z[is.na(v)] <- "<NA>"
    z
  })
  names(out) <- names(x)
  as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}

compare_tables <- function(label, rebuilt, current, current_path) {
  a <- canonical_table(rebuilt)
  b <- canonical_table(current)
  same_dim <- identical(dim(a), dim(b))
  same_columns <- identical(names(a), names(b))
  mismatches <- NA_real_
  if (same_dim && same_columns) {
    mismatches <- sum(as.matrix(a) != as.matrix(b))
  }
  hashes <- file_record(label, current_path)
  data.frame(
    artifact = label,
    rebuilt_rows = nrow(a),
    rebuilt_columns = ncol(a),
    current_rows = nrow(b),
    current_columns = ncol(b),
    dimensions_exact = same_dim,
    column_order_exact = same_columns,
    cell_mismatches = mismatches,
    semantic_exact = isTRUE(same_dim && same_columns && mismatches == 0),
    current_md5 = hashes$md5,
    current_sha256 = hashes$sha256,
    stringsAsFactors = FALSE
  )
}

a_level_mapping_audit <- function(probe_expr, maps, primary, sensitivity, gse_clin) {
  chosen_idx <- match(primary$mapping$PROBEID, rownames(probe_expr))
  chosen_exact <- identical(
    unname(primary$expr),
    unname(probe_expr[chosen_idx, , drop = FALSE])
  )
  checks <- data.frame(
    check_id = c(
      "annotation_has_one_row_per_observed_probe",
      "eligible_probes_have_exactly_one_SYMBOL",
      "primary_and_mean_have_same_gene_order",
      "primary_has_one_fixed_probe_per_gene",
      "primary_values_equal_selected_probe_rows",
      "primary_gene_order_is_lexicographic",
      "sample_order_preserved_in_primary",
      "sample_order_preserved_in_mean",
      "clinical_and_expression_sample_order_match",
      "primary_values_are_finite",
      "mean_sensitivity_values_are_finite"
    ),
    pass = c(
      nrow(maps$audit) == nrow(probe_expr) && identical(maps$audit$PROBEID, rownames(probe_expr)),
      nrow(maps$one) > 0L && all(maps$one$n_symbols == 1L) &&
        all(!is.na(maps$one$SYMBOL) & nzchar(maps$one$SYMBOL)),
      identical(rownames(primary$expr), rownames(sensitivity$expr)),
      nrow(primary$mapping) == nrow(primary$expr) &&
        !anyDuplicated(primary$mapping$SYMBOL) && !anyDuplicated(primary$mapping$PROBEID),
      chosen_exact,
      identical(rownames(primary$expr), sort(rownames(primary$expr))),
      identical(colnames(primary$expr), colnames(probe_expr)),
      identical(colnames(sensitivity$expr), colnames(probe_expr)),
      identical(gse_clin$sample, colnames(probe_expr)),
      all(is.finite(primary$expr)),
      all(is.finite(sensitivity$expr))
    ),
    observed = c(
      sprintf("%d annotation rows / %d probe rows", nrow(maps$audit), nrow(probe_expr)),
      sprintf("%d eligible probes", nrow(maps$one)),
      sprintf("%d / %d genes", nrow(primary$expr), nrow(sensitivity$expr)),
      sprintf("%d selected probes / %d genes", nrow(primary$mapping), nrow(primary$expr)),
      as.character(chosen_exact),
      as.character(identical(rownames(primary$expr), sort(rownames(primary$expr)))),
      sprintf("%d samples", ncol(primary$expr)),
      sprintf("%d samples", ncol(sensitivity$expr)),
      sprintf("%d clinical rows / %d samples", nrow(gse_clin), ncol(probe_expr)),
      sprintf("%d non-finite cells", sum(!is.finite(primary$expr))),
      sprintf("%d non-finite cells", sum(!is.finite(sensitivity$expr)))
    ),
    stringsAsFactors = FALSE
  )
  checks$note <- c(
    "Annotation audit includes unmapped and multi-SYMBOL probes.",
    "Eligibility is annotation-only and may be frozen before resampling.",
    "Both mappings must expose the same gene universe.",
    "The selected probe is fixed across samples for this reference matrix.",
    "Reference values must be copied from the selected probe without per-sample maximization.",
    "Deterministic output ordering.",
    "No expression sample reordering.",
    "No expression sample reordering.",
    "The raw public header anchors outcome rows to expression columns.",
    "Finite matrix required for modelling.",
    "Finite matrix required for modelling."
  )
  checks
}

write_expr_csv <- function(expr, path, id_name = "gene") {
  if (file.exists(path) && !force) stop("Refusing to overwrite: ", path, call. = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".part")
  if (file.exists(tmp)) unlink(tmp)
  con <- gzfile(tmp, open = "wt")
  ok <- FALSE
  on.exit({
    try(close(con), silent = TRUE)
    if (!ok && file.exists(tmp)) unlink(tmp)
  }, add = TRUE)
  ids <- rownames(expr)
  dat <- data.frame(ids, expr, check.names = FALSE)
  names(dat)[1] <- id_name
  write.csv(dat, con, row.names = FALSE)
  close(con)
  ok <- TRUE
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not finalize: ", path, call. = FALSE)
  invisible(path)
}

write_table <- function(x, path, sep = ",") {
  if (file.exists(path) && !force) stop("Refusing to overwrite: ", path, call. = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".part")
  if (file.exists(tmp)) unlink(tmp)
  fwrite(x, tmp, sep = sep, quote = FALSE, na = "")
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not finalize: ", path, call. = FALSE)
  invisible(path)
}

derive_coad_survival <- function(pancan_path, output_path) {
  dat <- fread(pancan_path, sep = "\t", data.table = FALSE)
  cancer_col <- "cancer type abbreviation"
  keep_cols <- c("sample", "_PATIENT", "OS", "OS.time", "DSS", "DSS.time", "DFI", "DFI.time", "PFI", "PFI.time", "Redaction")
  if (!cancer_col %in% names(dat) || !all(keep_cols %in% names(dat))) {
    stop("Unexpected TCGA Pan-Cancer survival schema")
  }
  coad <- dat[dat[[cancer_col]] == "COAD", keep_cols, drop = FALSE]
  fwrite(coad, output_path, sep = "\t", quote = FALSE, na = "")
  invisible(output_path)
}

if (action == "metadata-only") {
  cat("Metadata manifest written to:", file.path(provenance_dir, "source_manifest_runtime.csv"), "\n")
  quit(save = "no", status = 0L)
}

gse_raw <- find_input(
  raw_dir,
  c("GSE39582_series_matrix.txt.gz", "GSE39582_series_matrix.txt"),
  "GSE39582 series matrix"
)
tcga_raw <- find_input(raw_dir, c("COAD_HiSeqV2.gz", "HiSeqV2.gz", "COAD_HiSeqV2.txt"), "TCGA COAD HiSeqV2")
tcga_clinical_raw <- find_input(
  raw_dir,
  c("COAD_clinicalMatrix.tsv", "COAD_clinicalMatrix"),
  "TCGA COAD clinical matrix"
)
survival_path <- file.path(raw_dir, "COAD_survival.txt")
if (!file.exists(survival_path)) {
  pancan_path <- file.path(raw_dir, "TCGA_PANCAN_survival.tsv")
  if (!file.exists(pancan_path)) {
    stop("Missing COAD_survival.txt and TCGA_PANCAN_survival.tsv in ", raw_dir, call. = FALSE)
  }
  derive_coad_survival(pancan_path, survival_path)
}

hash_rows <- rbindlist(list(
  file_record("GSE39582_raw_series_matrix", gse_raw),
  file_record("TCGA_COAD_raw_HiSeqV2", tcga_raw),
  file_record("TCGA_COAD_survival", survival_path),
  file_record("TCGA_COAD_raw_clinical_matrix", tcga_clinical_raw)
), fill = TRUE)

gse_probe <- read_geo_probe_matrix(gse_raw)
gse_clinical_rebuilt <- read_geo_clinical(gse_raw, colnames(gse_probe$expr))
gse_maps <- annotation_maps(rownames(gse_probe$expr))
tcga_rebuilt <- prepare_tcga_current(tcga_raw, survival_path)
tcga_clinical_rebuilt <- prepare_tcga_clinical(tcga_clinical_raw, tcga_rebuilt$common)
legacy_result <- aggregate_geo(gse_probe$expr, "legacy_per_sample_max", maps = gse_maps)

mapping_summary <- data.frame(
  series = gse_probe$meta$series,
  platform = gse_probe$meta$platform,
  source_processing = gse_probe$meta$processing,
  raw_probes = nrow(gse_probe$expr),
  samples = ncol(gse_probe$expr),
  probes_with_nonempty_symbol = length(unique(gse_maps$edges$PROBEID)),
  probes_with_exactly_one_symbol = nrow(gse_maps$one),
  probes_with_multiple_symbols = sum(gse_maps$audit$n_symbols > 1L),
  probes_without_symbol = sum(gse_maps$audit$n_symbols == 0L),
  legacy_genes = nrow(legacy_result$expr),
  unique_probe_genes = length(unique(gse_maps$one$SYMBOL)),
  hgu133plus2_db_version = as.character(packageVersion("hgu133plus2.db")),
  AnnotationDbi_version = as.character(packageVersion("AnnotationDbi")),
  legacy_rule = paste(
    "First nonempty SYMBOL edge returned per PROBEID;",
    "within each SYMBOL take the maximum probe intensity separately in each sample"
  ),
  a_level_annotation_rule = paste(
    "Retain observed probes mapping to exactly one distinct nonempty SYMBOL;",
    "exclude unmapped and multi-SYMBOL probes"
  ),
  a_level_primary_selection_scope = paste(
    "Select the highest-mean eligible probe separately inside each modelling",
    "training partition; tie-break by lexicographic PROBEID"
  ),
  stringsAsFactors = FALSE
)
fwrite(mapping_summary, file.path(provenance_dir, "mapping_summary_runtime.csv"))

if (action %in% c("verify-current", "write-current")) {
  gse_current_rebuilt <- legacy_result$expr

  if (action == "verify-current") {
    gse_current_path <- file.path(derived_dir, "GSE39582_expression.csv.gz")
    tcga_current_path <- file.path(derived_dir, "TCGA_COAD_expression.csv.gz")
    gse_clinical_current_path <- file.path(derived_dir, "GSE39582_clin.csv")
    tcga_clinical_current_path <- file.path(
      current_results_dir, "tcga_clinical_annotations.csv"
    )
    current_gse <- read_derived_expr(gse_current_path)
    current_tcga <- read_derived_expr(tcga_current_path)
    current_gse_clinical <- read.csv(
      gse_clinical_current_path, check.names = FALSE, stringsAsFactors = FALSE,
      na.strings = character()
    )
    current_tcga_clinical <- read.csv(
      tcga_clinical_current_path, check.names = FALSE, stringsAsFactors = FALSE
    )
    audit <- rbind(
      compare_matrices("GSE39582", gse_current_rebuilt, current_gse, gse_current_path),
      compare_matrices("TCGA-COAD", tcga_rebuilt$expr, current_tcga, tcga_current_path)
    )
    clinical_audit <- rbind(
      compare_tables(
        "GSE39582_clinical", gse_clinical_rebuilt, current_gse_clinical,
        gse_clinical_current_path
      ),
      compare_tables(
        "TCGA_COAD_clinical_annotations", tcga_clinical_rebuilt,
        current_tcga_clinical, tcga_clinical_current_path
      )
    )
    fwrite(audit, file.path(provenance_dir, "reconstruction_audit_runtime.csv"))
    fwrite(
      clinical_audit,
      file.path(provenance_dir, "clinical_reconstruction_audit_runtime.csv")
    )
    hash_rows <- rbindlist(list(
      hash_rows,
      file_record("GSE39582_current_derived", gse_current_path),
      file_record("TCGA_COAD_current_derived", tcga_current_path),
      file_record("GSE39582_current_clinical", gse_clinical_current_path),
      file_record("TCGA_COAD_current_clinical_annotations", tcga_clinical_current_path)
    ), fill = TRUE)
    fwrite(hash_rows, file.path(provenance_dir, "local_hash_manifest_runtime.csv"))
    print(audit)
    print(clinical_audit)
    if (!all(audit$semantic_exact)) stop("At least one current matrix failed semantic reconstruction")
    if (!all(clinical_audit$semantic_exact)) {
      stop("At least one current clinical table failed semantic reconstruction")
    }
  } else {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write_expr_csv(gse_current_rebuilt, file.path(output_dir, "GSE39582_expression.csv.gz"))
    write_expr_csv(tcga_rebuilt$expr, file.path(output_dir, "TCGA_COAD_expression.csv.gz"))
  }
}

if (action %in% c("audit-a-level", "write-a-level")) {
  # This full-cohort highest-mean result is an input-contract reference only.
  # The strict primary nested-CV analysis must repeat this selection using only
  # the samples in each relevant training partition.
  full_cohort_reference <- aggregate_geo(
    gse_probe$expr, "unique_highest_mean", maps = gse_maps
  )
  mean_sensitivity <- aggregate_geo(gse_probe$expr, "unique_mean", maps = gse_maps)
  a_level_audit <- a_level_mapping_audit(
    gse_probe$expr, gse_maps, full_cohort_reference, mean_sensitivity,
    gse_clinical_rebuilt
  )
  fwrite(
    gse_maps$audit,
    file.path(provenance_dir, "GSE39582_probe_annotation_runtime.csv")
  )
  fwrite(
    full_cohort_reference$mapping,
    file.path(
      provenance_dir,
      "GSE39582_unique_highest_mean_full_cohort_probe_map_runtime.csv"
    )
  )
  fwrite(
    mean_sensitivity$mapping,
    file.path(provenance_dir, "GSE39582_unique_mean_probe_map_runtime.csv")
  )
  fwrite(
    a_level_audit,
    file.path(provenance_dir, "a_level_mapping_audit_runtime.csv")
  )
  fwrite(
    hash_rows,
    file.path(provenance_dir, "a_level_source_hash_manifest_runtime.csv")
  )
  print(a_level_audit)
  if (!all(a_level_audit$pass)) {
    stop(
      "A-level mapping audit failed: ",
      paste(a_level_audit$check_id[!a_level_audit$pass], collapse = ", ")
    )
  }

  if (action == "write-a-level") {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write_expr_csv(
      full_cohort_reference$expr,
      file.path(
        output_dir,
        "GSE39582_expression_unique_highest_mean_full_cohort_reference.csv.gz"
      )
    )
    write_expr_csv(
      mean_sensitivity$expr,
      file.path(output_dir, "GSE39582_expression_unique_mean.csv.gz")
    )
    write_expr_csv(
      tcga_rebuilt$expr,
      file.path(output_dir, "TCGA_COAD_expression.csv.gz")
    )
    write_table(
      gse_maps$audit,
      file.path(output_dir, "GSE39582_probe_annotation.csv")
    )
    write_table(
      gse_clinical_rebuilt,
      file.path(output_dir, "GSE39582_clin.csv")
    )
    write_table(
      fread(survival_path, sep = "\t", data.table = FALSE),
      file.path(output_dir, "COAD_survival.txt"),
      sep = "\t"
    )
    write_table(
      tcga_clinical_rebuilt,
      file.path(output_dir, "TCGA_COAD_clinical_annotations.csv")
    )
    prepared_paths <- c(
      file.path(
        output_dir,
        "GSE39582_expression_unique_highest_mean_full_cohort_reference.csv.gz"
      ),
      file.path(output_dir, "GSE39582_expression_unique_mean.csv.gz"),
      file.path(output_dir, "TCGA_COAD_expression.csv.gz"),
      file.path(output_dir, "GSE39582_probe_annotation.csv"),
      file.path(output_dir, "GSE39582_clin.csv"),
      file.path(output_dir, "COAD_survival.txt"),
      file.path(output_dir, "TCGA_COAD_clinical_annotations.csv")
    )
    prepared_roles <- c(
      "GSE39582_full_cohort_probe_selection_reference_not_primary_nested_input",
      "GSE39582_primary_unique_probe_mean",
      "TCGA_COAD_gene_expression",
      "GSE39582_frozen_probe_annotation",
      "GSE39582_clinical",
      "TCGA_COAD_survival",
      "TCGA_COAD_clinical_annotations"
    )
    prepared_manifest <- rbindlist(Map(file_record, prepared_roles, prepared_paths))
    fwrite(
      prepared_manifest,
      file.path(provenance_dir, "prepared_input_manifest_runtime.csv")
    )
  }
}

runtime <- data.frame(
  name = c(
    "script_version", "generated_at", "R.version.string", "data.table",
    "AnnotationDbi", "hgu133plus2.db", "digest"
  ),
  value = c(
    SCRIPT_VERSION,
    format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    R.version.string,
    as.character(packageVersion("data.table")),
    as.character(packageVersion("AnnotationDbi")),
    as.character(packageVersion("hgu133plus2.db")),
    as.character(packageVersion("digest"))
  ),
  stringsAsFactors = FALSE
)
fwrite(runtime, file.path(provenance_dir, "runtime_versions.csv"))
cat("Completed action:", action, "\n")
