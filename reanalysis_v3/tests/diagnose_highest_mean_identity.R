#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !nzchar(args[[1L]])) {
  stop("Usage: Rscript diagnose_highest_mean_identity.R <completed-run-directory>")
}
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Cannot locate the diagnostic script")
script_dir <- dirname(sub(
  "/+$", "",
  gsub("\\\\", "/", path.expand(sub("^--file=", "", file_arg)))
))
source(file.path(script_dir, "validate_completed_run.R"))

run_dir <- sub("/+$", "", gsub("\\\\", "/", path.expand(args[[1L]])))
status <- qa_read_csv(run_dir, "results/run_status.csv")
manifest <- qa_validate_manifest(run_dir, status)
params <- qa_validate_parameters(run_dir)
primary_nested <- qa_validate_nested(run_dir, manifest, status, params)
strict_nested <- qa_validate_strict_nested(
  run_dir, manifest, params, primary_nested$analysis_key
)

outer <- qa_read_csv(
  run_dir, "results/highest_mean_outer_probe_maps.csv",
  c(
    "repeat_id", "outer_fold", "SYMBOL", "PROBEID",
    "training_probe_mean", "training_probe_mean_hex", "n_probes_gene",
    "map_hash", "training_sample_hash", "training_n", "analysis_key"
  )
)
rows <- list()
for (repeat_id in seq_len(strict_nested$repeats)) {
  assignment <- strict_nested$outer_assignments[
    strict_nested$outer_assignments$repeat_id == repeat_id, , drop = FALSE
  ]
  for (outer_fold in seq_len(strict_nested$outer_folds)) {
    z <- outer[
      outer$repeat_id == repeat_id & outer$outer_fold == outer_fold,
      , drop = FALSE
    ]
    train_samples <- setdiff(
      strict_nested$sample_ids,
      assignment$sample[assignment$outer_fold == outer_fold]
    )
    expected_map_hash <- qa_hash_probe_map(z)
    expected_sample_hash <- qa_hash_sample_ids(train_samples)
    rows[[length(rows) + 1L]] <- data.frame(
      repeat_id = repeat_id,
      outer_fold = outer_fold,
      rows = nrow(z),
      stored_map_hash = unique(as.character(z$map_hash)),
      expected_map_hash = expected_map_hash,
      map_hash_match = identical(
        as.character(unique(z$map_hash)), expected_map_hash
      ),
      stored_training_n = unique(as.integer(z$training_n)),
      expected_training_n = length(train_samples),
      training_n_match = identical(
        unique(as.integer(z$training_n)), as.integer(length(train_samples))
      ),
      stored_sample_hash = unique(as.character(z$training_sample_hash)),
      expected_sample_hash = expected_sample_hash,
      sample_hash_match = identical(
        as.character(unique(z$training_sample_hash)), expected_sample_hash
      ),
      stringsAsFactors = FALSE
    )
  }
}
diagnosis <- do.call(rbind, rows)
print(diagnosis, row.names = FALSE)
cat("\nMismatch counts:\n")
print(colSums(!diagnosis[, c(
  "map_hash_match", "training_n_match", "sample_hash_match"
)]))

full_csv <- qa_read_csv(
  run_dir, "results/highest_mean_full_development_probe_map.csv"
)
full_model <- readRDS(qa_run_path(
  run_dir, "results/highest_mean_full_development_model.rds"
))
runtime_map <- qa_canonical_probe_map(full_model$probe_map$mapping)
csv_map <- qa_canonical_probe_map(full_csv)
cat("\nFull-map hash comparison:\n")
print(data.frame(
  source = c("stored", "runtime_RDS_rehash", "CSV_rehash"),
  hash = c(
    full_model$probe_map_hash,
    qa_hash_probe_map(runtime_map),
    qa_hash_probe_map(csv_map)
  ),
  stringsAsFactors = FALSE
), row.names = FALSE)
cat("\nCanonical RDS-versus-CSV all.equal:\n")
print(all.equal(runtime_map, csv_map, tolerance = 0, check.attributes = TRUE))
cat("\nCanonical column classes:\n")
print(rbind(
  runtime_RDS = vapply(runtime_map, function(x) paste(class(x), collapse = "/"), character(1)),
  CSV_reload = vapply(csv_map, function(x) paste(class(x), collapse = "/"), character(1))
))

snapshot_dir <- file.path(run_dir, "code_snapshot")
production_env <- new.env(parent = globalenv())
sys.source(file.path(snapshot_dir, "utils.R"), envir = production_env)
sys.source(
  file.path(snapshot_dir, "highest_mean_sensitivity.R"),
  envir = production_env
)
cat("\nProduction-snapshot rehash:\n")
print(data.frame(
  stored = full_model$probe_map_hash,
  production_RDS_rehash = production_env$hash_probe_map(
    full_model$probe_map$mapping
  ),
  validator_CSV_rehash = qa_hash_probe_map(full_csv),
  stringsAsFactors = FALSE
), row.names = FALSE)

original_collate <- Sys.getlocale("LC_COLLATE")
locale_hash_rows <- lapply(c("C", "Chinese_China.utf8"), function(locale_name) {
  locale_result <- suppressWarnings(Sys.setlocale("LC_COLLATE", locale_name))
  data.frame(
    requested_locale = locale_name,
    active_locale = as.character(locale_result),
    hash = production_env$hash_probe_map(full_model$probe_map$mapping),
    stringsAsFactors = FALSE
  )
})
suppressWarnings(Sys.setlocale("LC_COLLATE", original_collate))
cat("\nLocale-specific production rehash:\n")
print(do.call(rbind, locale_hash_rows), row.names = FALSE)

sample_ids <- as.character(full_model$inner_fold_assignments$sample)
probe_expr <- production_env$read_geo_matrix(file.path(
  dirname(dirname(run_dir)), "raw_inputs", "GSE39582_series_matrix.txt.gz"
))
probe_expr <- probe_expr[, sample_ids, drop = FALSE]
probe_annotation <- utils::read.csv(
  file.path(
    dirname(dirname(run_dir)), "prepared_inputs",
    "GSE39582_probe_annotation.csv"
  ),
  stringsAsFactors = FALSE, check.names = FALSE
)
probe_annotation$PROBEID <- as.character(probe_annotation$PROBEID)
probe_annotation$SYMBOL <- as.character(probe_annotation$SYMBOL)
coding <- !is.na(probe_annotation$SYMBOL) &
  nzchar(probe_annotation$SYMBOL) &
  !(grepl(
    "^LINC|^MIR|^SNORA|^SNORD|^RNU|^RN7S|^SCARNA|^RNY|^RPPH|^MIRLET",
    probe_annotation$SYMBOL
  ) |
    grepl("-AS[0-9]*$|-IT[0-9]*$", probe_annotation$SYMBOL) |
    grepl("^AC[0-9]{6}", probe_annotation$SYMBOL))
probe_annotation <- probe_annotation[coding, , drop = FALSE]
relearned <- production_env$learn_unique_highest_probe_map(
  probe_expr, probe_annotation, sample_ids
)
cat("\nRelearned full-development map comparison:\n")
print(data.frame(
  source = c("stored", "relearned", "persisted_rehash"),
  hash = c(
    full_model$probe_map_hash, relearned$map_hash,
    production_env$hash_probe_map(full_model$probe_map$mapping)
  ),
  stringsAsFactors = FALSE
), row.names = FALSE)
print(all.equal(
  production_env$.canonical_probe_mapping(relearned$mapping),
  production_env$.canonical_probe_mapping(full_model$probe_map$mapping),
  tolerance = 0, check.attributes = TRUE
))
