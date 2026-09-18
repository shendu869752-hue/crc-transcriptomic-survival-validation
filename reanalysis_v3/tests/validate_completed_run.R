#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

.qa_sha256_pattern <- "^[0-9a-f]{64}$"
.qa_probe_map_schema <- "unique_highest_mean_probe_map_v4_radix_textsha256"
.qa_probe_model_schema <- "probe_signature_model_v1"
.qa_probe_sample_schema <- "probe_training_sample_set_v1"

qa_fail <- function(...) stop(paste0(...), call. = FALSE)

qa_assert <- function(value, ...) {
  if (length(value) != 1L || is.na(value) || !isTRUE(value)) qa_fail(...)
  invisible(TRUE)
}

qa_sha256_file <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    qa_fail("Package 'digest' is required for completed-run validation")
  }
  qa_assert(file.exists(path) && !dir.exists(path), "Missing file: ", path)
  digest::digest(file = path, algo = "sha256")
}

qa_normalize_relative <- function(path) {
  path <- gsub("\\\\", "/", as.character(path))
  bad <- is.na(path) | !nzchar(path) |
    grepl("^[A-Za-z]:/|^/|^//", path) |
    grepl("(^|/)\\.\\.(/|$)", path)
  if (any(bad)) qa_fail("Unsafe relative artifact path: ", paste(path[bad], collapse = ", "))
  path
}

qa_run_path <- function(run_dir, relative_path) {
  file.path(run_dir, qa_normalize_relative(relative_path))
}

qa_read_csv <- function(run_dir, relative_path, required_columns = character(),
                        allow_empty = FALSE) {
  path <- qa_run_path(run_dir, relative_path)
  qa_assert(file.exists(path) && !dir.exists(path),
            "Required CSV is missing: ", relative_path)
  qa_assert(is.finite(file.info(path)$size) && file.info(path)$size > 0,
            "Required CSV is empty: ", relative_path)
  out <- tryCatch(
    utils::read.csv(
      path, stringsAsFactors = FALSE, check.names = FALSE,
      na.strings = c("NA", "")
    ),
    error = function(e) qa_fail("Cannot parse ", relative_path, ": ", conditionMessage(e))
  )
  missing <- setdiff(required_columns, names(out))
  if (length(missing)) {
    qa_fail(relative_path, " lacks columns: ", paste(missing, collapse = ", "))
  }
  if (!allow_empty && !nrow(out)) qa_fail(relative_path, " has no data rows")
  out
}

qa_character <- function(x, label, allow_blank = FALSE) {
  out <- as.character(x)
  bad <- is.na(out) | (!allow_blank & !nzchar(trimws(out)))
  if (any(bad)) qa_fail(label, " contains missing or blank values")
  out
}

qa_numeric <- function(x, label, finite = TRUE) {
  out <- suppressWarnings(as.numeric(x))
  if (anyNA(out) || (finite && any(!is.finite(out)))) {
    qa_fail(label, " contains invalid numeric values")
  }
  out
}

qa_all_close <- function(observed, expected, tolerance = 1e-10) {
  observed <- as.numeric(observed)
  expected <- as.numeric(expected)
  length(observed) == length(expected) &&
    !anyNA(observed) && !anyNA(expected) &&
    all(is.finite(observed)) && all(is.finite(expected)) &&
    all(abs(observed - expected) <= tolerance * pmax(1, abs(expected)))
}

qa_integer <- function(x, label, minimum = NULL) {
  out <- qa_numeric(x, label)
  if (any(out != as.integer(out))) qa_fail(label, " is not integer-valued")
  out <- as.integer(out)
  if (!is.null(minimum) && any(out < minimum)) qa_fail(label, " is below ", minimum)
  out
}

qa_logical <- function(x, label, allow_na = FALSE) {
  if (is.logical(x)) {
    out <- x
  } else {
    z <- tolower(trimws(as.character(x)))
    out <- rep(NA, length(z))
    out[z %in% c("true", "t", "1")] <- TRUE
    out[z %in% c("false", "f", "0")] <- FALSE
  }
  if (!allow_na && anyNA(out)) qa_fail(label, " contains invalid logical values")
  out
}

qa_require_named_list <- function(x, label, required_names) {
  qa_assert(is.list(x) && !is.data.frame(x), label, " must be a list")
  nm <- names(x)
  qa_assert(!is.null(nm) && !anyNA(nm) && all(nzchar(nm)) && !anyDuplicated(nm),
            label, " must have unique non-empty names")
  missing <- setdiff(required_names, nm)
  if (length(missing)) qa_fail(label, " lacks names: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

qa_split_semicolon <- function(value) {
  if (length(value) != 1L || is.na(value) || !nzchar(trimws(value))) return(character())
  out <- trimws(strsplit(as.character(value), ";", fixed = TRUE)[[1L]])
  out[nzchar(out)]
}

qa_canonical_probe_map <- function(x) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  required <- c(
    "SYMBOL", "PROBEID", "training_probe_mean", "training_probe_mean_hex",
    "n_probes_gene"
  )
  missing <- setdiff(required, names(x))
  if (length(missing)) qa_fail("Probe map lacks columns: ", paste(missing, collapse = ", "))
  x <- x[, required, drop = FALSE]
  x$SYMBOL <- qa_character(x$SYMBOL, "Probe-map SYMBOL")
  x$PROBEID <- qa_character(x$PROBEID, "Probe-map PROBEID")
  reported_mean <- qa_numeric(x$training_probe_mean, "Probe-map training mean")
  x$training_probe_mean_hex <- qa_character(
    x$training_probe_mean_hex, "Probe-map exact training mean"
  )
  valid_hex <- grepl(
    "^hex:[+-]?0x[0-9a-f]+(?:\\.[0-9a-f]*)?p[+-]?[0-9]+$",
    x$training_probe_mean_hex,
    ignore.case = TRUE,
    perl = TRUE
  )
  exact_mean <- suppressWarnings(as.numeric(sub(
    "^hex:", "", x$training_probe_mean_hex
  )))
  qa_assert(
    all(valid_hex) && all(is.finite(exact_mean)) &&
      all(abs(reported_mean - exact_mean) <=
          1e-12 * pmax(1, abs(exact_mean))),
    "Probe-map exact and reported training means disagree"
  )
  x$training_probe_mean <- exact_mean
  x$training_probe_mean_hex <- paste0("hex:", sprintf("%a", exact_mean))
  x$n_probes_gene <- qa_integer(x$n_probes_gene, "Probe-map eligible-probe count", 1L)
  qa_assert(!anyDuplicated(x$SYMBOL) && !anyDuplicated(x$PROBEID),
            "Probe map is not one-to-one")
  x <- x[order(x$SYMBOL, x$PROBEID, method = "radix"), , drop = FALSE]
  rownames(x) <- NULL
  x
}

qa_hash_probe_map <- function(x) {
  mapping <- qa_canonical_probe_map(x)
  qa_assert(
    !any(grepl("[[:cntrl:]]", mapping$SYMBOL)) &&
      !any(grepl("[[:cntrl:]]", mapping$PROBEID)),
    "Probe-map identifiers contain tab or newline characters"
  )
  rows <- paste(
    enc2utf8(mapping$SYMBOL), enc2utf8(mapping$PROBEID),
    mapping$training_probe_mean_hex,
    sprintf("%d", mapping$n_probes_gene),
    sep = "\t"
  )
  payload <- paste(
    c(
      .qa_probe_map_schema,
      "SYMBOL\tPROBEID\ttraining_probe_mean_hex\tn_probes_gene",
      rows
    ),
    collapse = "\n"
  )
  digest::digest(payload, algo = "sha256", serialize = FALSE)
}

qa_hash_sample_ids <- function(sample_ids) {
  sample_ids <- qa_character(sample_ids, "Training sample IDs")
  qa_assert(!anyDuplicated(sample_ids), "Training sample IDs are duplicated")
  digest::digest(
    list(schema = .qa_probe_sample_schema, sample_ids = sort(sample_ids)),
    algo = "sha256", serialize = TRUE
  )
}

qa_hash_gene_model <- function(model) {
  required <- c(
    "genes", "coefs", "center", "scale", "train_lp_mean", "train_lp_sd",
    "lambda_ratio", "lambda"
  )
  missing <- setdiff(required, names(model))
  if (length(missing)) qa_fail("Highest-mean model lacks fields: ", paste(missing, collapse = ", "))
  genes <- qa_character(model$genes, "Highest-mean model genes")
  qa_assert(!anyDuplicated(genes), "Highest-mean model genes are duplicated")
  qa_assert(all(genes %in% names(model$coefs)) &&
              all(genes %in% names(model$center)) &&
              all(genes %in% names(model$scale)),
            "Highest-mean model feature metadata is incomplete")
  digest::digest(list(
    genes = genes,
    coefficients = unname(model$coefs[genes]),
    center = unname(model$center[genes]),
    scale = unname(model$scale[genes]),
    lp_center = model$train_lp_mean,
    lp_scale = model$train_lp_sd,
    lambda_ratio = model$lambda_ratio,
    lambda = model$lambda
  ), algo = "sha256", serialize = TRUE)
}

qa_combine_model_fingerprint <- function(gene_model_fingerprint, probe_map_hash) {
  digest::digest(list(
    schema = .qa_probe_model_schema,
    gene_model_fingerprint = gene_model_fingerprint,
    probe_map_hash = probe_map_hash
  ), algo = "sha256", serialize = TRUE)
}

qa_require_manifested <- function(manifest, relative_paths) {
  relative_paths <- qa_normalize_relative(relative_paths)
  idx <- match(relative_paths, manifest$relative_path)
  if (anyNA(idx)) {
    qa_fail("Required artifacts are absent from output_manifest.csv: ",
            paste(relative_paths[is.na(idx)], collapse = ", "))
  }
  bad <- manifest$status[idx] != "generated" | !manifest$required_value[idx]
  if (any(bad)) {
    qa_fail("Required artifacts are not required/generated in output_manifest.csv: ",
            paste(relative_paths[bad], collapse = ", "))
  }
  invisible(TRUE)
}

qa_validate_manifest <- function(run_dir, status) {
  manifest_rel <- "results/output_manifest.csv"
  manifest <- qa_read_csv(
    run_dir, manifest_rel,
    c("artifact_id", "relative_path", "required", "status", "reason", "bytes", "sha256")
  )
  manifest$artifact_id <- qa_character(manifest$artifact_id, "Manifest artifact_id")
  manifest$relative_path <- qa_normalize_relative(manifest$relative_path)
  manifest$required_value <- qa_logical(manifest$required, "Manifest required")
  manifest$status <- qa_character(manifest$status, "Manifest status")
  qa_assert(!anyDuplicated(manifest$artifact_id), "Manifest artifact_id values are duplicated")
  qa_assert(!anyDuplicated(manifest$relative_path), "Manifest paths are duplicated")
  qa_assert(all(manifest$status %in% c("generated", "not_generated")),
            "Manifest contains unresolved artifact status")
  controls <- c("results/run_status.csv", manifest_rel)
  qa_assert(!any(manifest$relative_path %in% controls),
            "Control files must not be self-listed in output_manifest.csv")

  generated <- manifest$status == "generated"
  qa_assert(!any(manifest$required_value & !generated),
            "A required manifest artifact is not generated")
  generated_hash <- as.character(manifest$sha256[generated])
  qa_assert(all(grepl(.qa_sha256_pattern, generated_hash)),
            "A generated artifact has an invalid SHA-256")
  generated_bytes <- qa_numeric(manifest$bytes[generated], "Generated artifact bytes")
  qa_assert(all(generated_bytes >= 0), "A generated artifact has negative bytes")
  for (i in which(generated)) {
    path <- qa_run_path(run_dir, manifest$relative_path[[i]])
    qa_assert(file.exists(path) && !dir.exists(path),
              "Manifest-generated artifact is missing: ", manifest$relative_path[[i]])
    qa_assert(as.numeric(file.info(path)$size) == as.numeric(manifest$bytes[[i]]),
              "Manifest byte mismatch: ", manifest$relative_path[[i]])
    qa_assert(identical(qa_sha256_file(path), as.character(manifest$sha256[[i]])),
              "Manifest SHA-256 mismatch: ", manifest$relative_path[[i]])
  }
  not_generated <- !generated
  if (any(not_generated)) {
    reasons <- as.character(manifest$reason[not_generated])
    qa_assert(!anyNA(reasons) && all(nzchar(trimws(reasons))),
              "A not_generated artifact lacks a reason")
    for (rel in manifest$relative_path[not_generated]) {
      qa_assert(!file.exists(qa_run_path(run_dir, rel)),
                "A not_generated artifact exists on disk: ", rel)
    }
  }

  actual <- list.files(
    run_dir, recursive = TRUE, full.names = FALSE, all.files = TRUE, no.. = TRUE
  )
  if (length(actual)) {
    is_file <- !file.info(file.path(run_dir, actual))$isdir
    actual <- qa_normalize_relative(actual[is_file])
  }
  expected <- c(manifest$relative_path[generated], controls)
  unexpected <- setdiff(actual, expected)
  missing <- setdiff(expected, actual)
  if (length(unexpected) || length(missing)) {
    qa_fail(
      "Run directory and output manifest differ; unexpected=[",
      paste(unexpected, collapse = ", "), "]; missing=[",
      paste(missing, collapse = ", "), "]"
    )
  }

  chain <- c(
    input_manifest_sha256 = "results/input_manifest.csv",
    parameter_manifest_sha256 = "results/run_parameters.txt",
    code_manifest_sha256 = "results/code_manifest.csv",
    output_manifest_sha256 = manifest_rel
  )
  for (field in names(chain)) {
    recorded <- as.character(status[[field]][[1L]])
    qa_assert(grepl(.qa_sha256_pattern, recorded), "Invalid status hash field: ", field)
    qa_assert(identical(recorded, qa_sha256_file(qa_run_path(run_dir, chain[[field]]))),
              "Broken four-hash chain at ", field)
  }
  manifest
}

qa_validate_code_snapshot <- function(run_dir, manifest) {
  input_manifest <- qa_read_csv(
    run_dir, "results/input_manifest.csv",
    c("source_id", "relative_path", "bytes", "sha256")
  )
  qa_assert(!anyDuplicated(qa_character(input_manifest$source_id, "Input source_id")),
            "Input manifest source_id values are duplicated")
  qa_normalize_relative(input_manifest$relative_path)
  qa_assert(all(qa_numeric(input_manifest$bytes, "Input manifest bytes") >= 0),
            "Input manifest has negative bytes")
  qa_assert(all(grepl(.qa_sha256_pattern, as.character(input_manifest$sha256))),
            "Input manifest contains an invalid SHA-256")

  code <- qa_read_csv(
    run_dir, "results/code_manifest.csv", c("category", "name", "value")
  )
  code$category <- qa_character(code$category, "Code-manifest category")
  code$name <- qa_character(code$name, "Code-manifest name")
  code$value <- qa_character(code$value, "Code-manifest value")
  qa_assert(!anyDuplicated(paste(code$category, code$name, sep = "\r")),
            "Code manifest contains duplicated category/name rows")
  hashes <- code[code$category == "code_sha256", , drop = FALSE]
  paths <- code[code$category == "code_snapshot_relative_path", , drop = FALSE]
  bytes <- code[code$category == "code_snapshot_bytes", , drop = FALSE]
  qa_assert(nrow(hashes) > 0L, "Code manifest contains no code_sha256 rows")
  qa_assert(setequal(hashes$name, paths$name) && setequal(hashes$name, bytes$name),
            "Code snapshot rows do not cover every code_sha256 row")
  qa_assert(all(grepl(.qa_sha256_pattern, hashes$value)),
            "Code manifest contains an invalid code SHA-256")
  for (name in hashes$name) {
    expected_hash <- hashes$value[hashes$name == name]
    relative_path <- qa_normalize_relative(paths$value[paths$name == name])
    expected_bytes <- qa_numeric(bytes$value[bytes$name == name], "Code snapshot bytes")
    qa_assert(length(relative_path) == 1L && startsWith(relative_path, "code_snapshot/"),
              "Invalid code snapshot path for ", name)
    qa_require_manifested(manifest, relative_path)
    path <- qa_run_path(run_dir, relative_path)
    qa_assert(identical(qa_sha256_file(path), expected_hash),
              "Code snapshot differs from code_manifest code_sha256: ", name)
    qa_assert(as.numeric(file.info(path)$size) == expected_bytes,
              "Code snapshot byte mismatch: ", name)
  }
  package_rows <- code[
    code$category == "package_version" & code$name == "org.Hs.eg.db", , drop = FALSE
  ]
  qa_assert(nrow(package_rows) == 1L, "Code manifest lacks one org.Hs.eg.db version")
  as.character(package_rows$value[[1L]])
}

qa_validate_parameters <- function(run_dir) {
  path <- qa_run_path(run_dir, "results/run_parameters.txt")
  params <- tryCatch(dget(path), error = function(e) {
    qa_fail("run_parameters.txt cannot be dget(): ", conditionMessage(e))
  })
  qa_require_named_list(params, "Run parameters", c(
    "output_contract_version", "pipeline_version", "nested_pipeline",
    "highest_mean_nested_pipeline", "bootstrap", "probe_mapping",
    "gene_identity_resolution", "scoring", "time_dependent_auc", "survival",
    "execution"
  ))
  nested_required <- c(
    "pipeline_version", "outer_folds", "outer_repeats", "inner_folds", "base_seed",
    "sd_cutoff", "univ_p_cutoff", "min_candidates", "candidate_fallback_n",
    "coef_cutoff", "min_selected", "selected_fallback_n", "max_selected",
    "lambda_ratio_grid", "lambda_rule", "ties", "development_probe_aggregation"
  )
  qa_require_named_list(params$nested_pipeline, "nested_pipeline", nested_required)
  qa_require_named_list(
    params$highest_mean_nested_pipeline, "highest_mean_nested_pipeline",
    c(nested_required, "external_probe_application")
  )
  qa_require_named_list(params$bootstrap, "bootstrap", c(
    "replicates", "percentile_probabilities", "cohort_seed_base",
    "minimum_valid_fraction", "estimand"
  ))
  qa_require_named_list(params$probe_mapping, "probe_mapping", c(
    "primary", "sensitivity", "eligibility", "primary_multi_probe_rule",
    "sensitivity_tie_break", "annotation", "external_sensitivity_scope"
  ))
  qa_require_named_list(params$gene_identity_resolution, "gene_identity_resolution", c(
    "scope", "rule", "exact_symbol_priority", "alias_requirements", "annotation",
    "annotation_version", "outcome_blind", "expression_values_used",
    "development_candidate_space_filtered_by_TCGA"
  ))
  qa_require_named_list(params$survival, "survival", c(
    "cox_ties", "spline_df", "clinical_min_n", "clinical_min_events",
    "os_meta_method", "os_meta_primary_test"
  ))
  qa_require_named_list(params$execution, "execution", c("force_nested", "entrypoint"))
  pipeline_version <- qa_character(params$pipeline_version, "pipeline_version")
  primary_version <- qa_character(
    params$nested_pipeline$pipeline_version, "nested_pipeline pipeline_version"
  )
  strict_version <- qa_character(
    params$highest_mean_nested_pipeline$pipeline_version,
    "highest_mean_nested_pipeline pipeline_version"
  )
  qa_assert(length(pipeline_version) == 1L && length(primary_version) == 1L &&
              identical(pipeline_version, primary_version),
            "Top-level and nested pipeline_version differ")
  qa_assert(length(strict_version) == 1L && !identical(strict_version, primary_version),
            "Strict highest-mean pipeline_version must differ from primary")
  qa_assert(identical(as.character(params$probe_mapping$primary), "unique_mean"),
            "Primary probe mapping is not unique_mean")
  qa_assert(identical(as.character(params$probe_mapping$sensitivity), "unique_highest_mean"),
            "Prespecified unique_highest_mean sensitivity is absent")
  scope <- tolower(gsub("[ _-]", "", as.character(
    params$probe_mapping$external_sensitivity_scope
  )))
  qa_assert(length(scope) == 1L && grepl("fulldevelopment", scope) && grepl("frozen", scope),
            "External highest-mean scope is not a full-development frozen map")
  qa_assert(identical(as.character(params$survival$os_meta_method), "REML"),
            "OS meta-analysis method is not REML")
  qa_assert(tolower(as.character(params$survival$os_meta_primary_test)) %in%
              c("knha", "knapp-hartung", "knapp_hartung"),
            "OS meta-analysis primary test is not Knapp-Hartung")
  params
}

qa_validate_nested_contract <- function(
    run_dir, manifest, pipeline_params, expected_version, rel,
    expected_analysis_key = NULL,
    forbidden_analysis_key = NULL, label = "Nested CV") {
  qa_assert(length(rel) == 7L, label, " artifact contract must contain seven files")
  qa_require_manifested(manifest, rel)
  summary <- qa_read_csv(run_dir, rel[[1L]], c(
    "pipeline_version", "analysis_key", "repeats", "outer_folds", "inner_folds",
    "completed_folds"
  ))
  qa_assert(nrow(summary) == 1L, label, " summary must have one row")
  nested <- pipeline_params
  repeats <- qa_integer(nested$outer_repeats, "outer_repeats", 1L)
  outer_folds <- qa_integer(nested$outer_folds, "outer_folds", 2L)
  inner_folds <- qa_integer(nested$inner_folds, "inner_folds", 2L)
  analysis_key <- as.character(summary$analysis_key[[1L]])
  qa_assert(grepl(.qa_sha256_pattern, analysis_key), label, " has an invalid analysis_key")
  if (!is.null(expected_analysis_key)) {
    qa_assert(identical(analysis_key, as.character(expected_analysis_key)),
              label, " analysis_key differs from the expected run identity")
  }
  if (!is.null(forbidden_analysis_key)) {
    qa_assert(!identical(analysis_key, as.character(forbidden_analysis_key)),
              "Strict nested analysis_key must be independent from primary nested analysis_key")
  }
  qa_assert(identical(as.character(summary$pipeline_version[[1L]]),
                      as.character(expected_version)),
            label, " summary pipeline_version differs from its pipeline parameters")
  qa_assert(qa_integer(summary$repeats, "Nested summary repeats") == repeats &&
              qa_integer(summary$outer_folds, "Nested summary outer_folds") == outer_folds &&
              qa_integer(summary$inner_folds, "Nested summary inner_folds") == inner_folds &&
              qa_integer(summary$completed_folds, "Nested completed_folds") == repeats * outer_folds,
            "Nested summary dimensions are inconsistent with run parameters")

  oof <- qa_read_csv(run_dir, rel[[2L]], c(
    "sample", "repeat_id", "outer_fold", "analysis_key"
  ))
  outer <- qa_read_csv(run_dir, rel[[3L]], c(
    "sample", "repeat_id", "outer_fold", "analysis_key"
  ))
  inner <- qa_read_csv(run_dir, rel[[4L]], c(
    "sample", "inner_fold", "repeat_id", "outer_fold", "analysis_key"
  ))
  fold_perf <- qa_read_csv(run_dir, rel[[5L]], c(
    "repeat_id", "outer_fold", "analysis_key"
  ))
  repeat_perf <- qa_read_csv(run_dir, rel[[6L]], c("repeat_id", "analysis_key"))
  failures <- qa_read_csv(
    run_dir, rel[[7L]], c("repeat_id", "outer_fold", "stage", "message"),
    allow_empty = TRUE
  )
  qa_assert(nrow(failures) == 0L, "nested_cv_failure_log.csv contains failures")

  for (tab in list(oof, outer, inner, fold_perf, repeat_perf)) {
    qa_assert(all(as.character(tab$analysis_key) == analysis_key),
              "A nested-CV artifact has the wrong analysis_key")
  }
  oof$repeat_id <- qa_integer(oof$repeat_id, "OOF repeat_id", 1L)
  oof$outer_fold <- qa_integer(oof$outer_fold, "OOF outer_fold", 1L)
  outer$repeat_id <- qa_integer(outer$repeat_id, "Outer assignment repeat_id", 1L)
  outer$outer_fold <- qa_integer(outer$outer_fold, "Outer assignment fold", 1L)
  inner$repeat_id <- qa_integer(inner$repeat_id, "Inner assignment repeat_id", 1L)
  inner$outer_fold <- qa_integer(inner$outer_fold, "Inner assignment outer_fold", 1L)
  inner$inner_fold <- qa_integer(inner$inner_fold, "Inner assignment inner_fold", 1L)
  oof$sample <- qa_character(oof$sample, "OOF sample")
  outer$sample <- qa_character(outer$sample, "Outer assignment sample")
  inner$sample <- qa_character(inner$sample, "Inner assignment sample")

  expected_repeats <- seq_len(repeats)
  qa_assert(identical(sort(unique(oof$repeat_id)), expected_repeats),
            "OOF repeat IDs are incomplete")
  qa_assert(nrow(oof) == 573L * repeats,
            "OOF row count is not 573 x outer_repeats")
  qa_assert(!anyDuplicated(paste(oof$repeat_id, oof$sample, sep = "\r")),
            "OOF sample/repeat coverage is not unique")
  sample_ids <- sort(unique(oof$sample[oof$repeat_id == expected_repeats[[1L]]]))
  qa_assert(length(sample_ids) == 573L, "Development OOF sample count is not 573")
  for (repeat_id in expected_repeats) {
    z <- oof[oof$repeat_id == repeat_id, , drop = FALSE]
    qa_assert(setequal(z$sample, sample_ids), "OOF sample coverage differs by repeat")
    qa_assert(identical(sort(unique(z$outer_fold)), seq_len(outer_folds)),
              "OOF outer folds are incomplete")
  }

  qa_assert(nrow(outer) == 573L * repeats &&
              !anyDuplicated(paste(outer$repeat_id, outer$sample, sep = "\r")),
            "Outer-fold assignments do not uniquely cover 573 samples per repeat")
  oof_key <- paste(oof$repeat_id, oof$sample, sep = "\r")
  outer_key <- paste(outer$repeat_id, outer$sample, sep = "\r")
  matched <- match(oof_key, outer_key)
  qa_assert(!anyNA(matched) && all(oof$outer_fold == outer$outer_fold[matched]),
            "OOF outer folds differ from the assignment table")
  for (repeat_id in expected_repeats) {
    z <- outer[outer$repeat_id == repeat_id, , drop = FALSE]
    qa_assert(setequal(z$sample, sample_ids) &&
                identical(sort(unique(z$outer_fold)), seq_len(outer_folds)),
              "Outer-fold assignment coverage is incomplete")
    for (outer_fold in seq_len(outer_folds)) {
      expected_train <- setdiff(sample_ids, z$sample[z$outer_fold == outer_fold])
      current <- inner[
        inner$repeat_id == repeat_id & inner$outer_fold == outer_fold, , drop = FALSE
      ]
      qa_assert(!anyDuplicated(current$sample) && setequal(current$sample, expected_train),
                "Inner assignments do not exactly cover the outer-training samples")
      qa_assert(identical(sort(unique(current$inner_fold)), seq_len(inner_folds)),
                "Inner folds are incomplete within an outer split")
    }
  }
  qa_assert(!anyDuplicated(paste(
    inner$repeat_id, inner$outer_fold, inner$sample, sep = "\r"
  )), "Inner assignment rows are duplicated")

  fold_perf$repeat_id <- qa_integer(fold_perf$repeat_id, "Fold performance repeat_id", 1L)
  fold_perf$outer_fold <- qa_integer(fold_perf$outer_fold, "Fold performance outer_fold", 1L)
  qa_assert(nrow(fold_perf) == repeats * outer_folds &&
              !anyDuplicated(paste(fold_perf$repeat_id, fold_perf$outer_fold, sep = "\r")),
            "Nested fold-performance coverage is incomplete")
  repeat_perf$repeat_id <- qa_integer(repeat_perf$repeat_id, "Repeat performance repeat_id", 1L)
  qa_assert(nrow(repeat_perf) == repeats &&
              identical(sort(repeat_perf$repeat_id), expected_repeats),
            "Nested repeat-performance coverage is incomplete")
  list(
    analysis_key = analysis_key, sample_ids = sample_ids, repeats = repeats,
    outer_folds = outer_folds, inner_folds = inner_folds,
    outer_assignments = outer, inner_assignments = inner
  )
}

qa_validate_nested <- function(run_dir, manifest, status, params) {
  qa_validate_nested_contract(
    run_dir, manifest,
    pipeline_params = params$nested_pipeline,
    expected_version = params$nested_pipeline$pipeline_version,
    rel = c(
      "results/nested_cv_summary.csv",
      "results/nested_cv_oof_predictions.csv",
      "results/nested_cv_outer_fold_assignments.csv",
      "results/nested_cv_inner_fold_assignments.csv",
      "results/nested_cv_fold_performance.csv",
      "results/nested_cv_repeat_performance.csv",
      "results/nested_cv_failure_log.csv"
    ),
    expected_analysis_key = as.character(status$nested_analysis_key[[1L]]),
    label = "Primary nested CV"
  )
}

qa_validate_strict_nested <- function(run_dir, manifest, params,
                                      primary_analysis_key) {
  qa_validate_nested_contract(
    run_dir, manifest,
    pipeline_params = params$highest_mean_nested_pipeline,
    expected_version = params$highest_mean_nested_pipeline$pipeline_version,
    rel = c(
      "results/highest_mean_nested_cv_summary.csv",
      "results/highest_mean_nested_cv_oof_predictions.csv",
      "results/highest_mean_nested_cv_outer_fold_assignments.csv",
      "results/highest_mean_nested_cv_inner_fold_assignments.csv",
      "results/highest_mean_nested_cv_fold_performance.csv",
      "results/highest_mean_nested_cv_repeat_performance.csv",
      "results/highest_mean_nested_cv_failure_log.csv"
    ),
    forbidden_analysis_key = primary_analysis_key,
    label = "Strict highest-mean nested CV"
  )
}

qa_validate_mapping_audit <- function(run_dir, manifest, relative_path,
                                      coefficient_path, annotation_version) {
  qa_require_manifested(manifest, c(relative_path, coefficient_path))
  coefficients <- qa_read_csv(run_dir, coefficient_path, c("gene", "coef"))
  coefficients$gene <- qa_character(coefficients$gene, paste(coefficient_path, "gene"))
  qa_assert(!anyDuplicated(coefficients$gene), coefficient_path, " has duplicated genes")
  audit <- qa_read_csv(run_dir, relative_path, c(
    "required_gene", "required_entrez_id", "source_row_id", "resolution_method",
    "exact_symbol_present", "expression_alias_candidates",
    "expression_alias_candidate_count", "reverse_unique", "symbol_consistent",
    "entrez_consistent", "reverse_symbols", "reverse_entrez_ids",
    "reverse_symbol_count", "reverse_entrez_count", "annotation_package",
    "annotation_version", "outcome_or_expression_values_used_for_resolution"
  ))
  qa_assert(identical(as.character(audit$required_gene), coefficients$gene),
            relative_path, " does not preserve coefficient-gene order")
  qa_assert(!anyDuplicated(audit$required_gene) && !anyDuplicated(audit$source_row_id),
            relative_path, " reuses a required gene or source row")
  exact <- qa_logical(audit$exact_symbol_present, paste(relative_path, "exact_symbol_present"))
  used_values <- qa_logical(
    audit$outcome_or_expression_values_used_for_resolution,
    paste(relative_path, "outcome_or_expression_values_used_for_resolution")
  )
  entrez_ok <- qa_logical(audit$entrez_consistent, paste(relative_path, "entrez_consistent"))
  reverse_unique <- qa_logical(
    audit$reverse_unique, paste(relative_path, "reverse_unique"), allow_na = TRUE
  )
  symbol_ok <- qa_logical(
    audit$symbol_consistent, paste(relative_path, "symbol_consistent"), allow_na = TRUE
  )
  qa_assert(!any(used_values), relative_path, " used outcomes or expression values for resolution")
  qa_assert(all(as.character(audit$annotation_package) == "org.Hs.eg.db") &&
              all(as.character(audit$annotation_version) == annotation_version),
            relative_path, " annotation identity differs from code_manifest.csv")
  qa_assert(all(audit$resolution_method %in% c(
    "exact_SYMBOL", "unique_reverse_unique_ENTREZ_alias"
  )), relative_path, " has an unsupported resolution_method")
  counts <- qa_integer(
    audit$expression_alias_candidate_count,
    paste(relative_path, "expression_alias_candidate_count"), 0L
  )
  for (i in seq_len(nrow(audit))) {
    candidates <- qa_split_semicolon(audit$expression_alias_candidates[[i]])
    qa_assert(length(unique(candidates)) == counts[[i]],
              relative_path, " alias-candidate field/count mismatch at row ", i)
    if (exact[[i]]) {
      qa_assert(audit$resolution_method[[i]] == "exact_SYMBOL" &&
                  audit$source_row_id[[i]] == audit$required_gene[[i]] && entrez_ok[[i]],
                relative_path, " exact-SYMBOL row is internally inconsistent")
    } else {
      reverse_symbols <- qa_split_semicolon(audit$reverse_symbols[[i]])
      reverse_entrez <- qa_split_semicolon(audit$reverse_entrez_ids[[i]])
      qa_assert(audit$resolution_method[[i]] == "unique_reverse_unique_ENTREZ_alias" &&
                  counts[[i]] == 1L &&
                  identical(candidates, as.character(audit$source_row_id[[i]])) &&
                  isTRUE(reverse_unique[[i]]) && isTRUE(symbol_ok[[i]]) &&
                  isTRUE(entrez_ok[[i]]) &&
                  identical(reverse_symbols, as.character(audit$required_gene[[i]])) &&
                  identical(reverse_entrez, as.character(audit$required_entrez_id[[i]])) &&
                  qa_integer(audit$reverse_symbol_count[[i]], "reverse_symbol_count") == 1L &&
                  qa_integer(audit$reverse_entrez_count[[i]], "reverse_entrez_count") == 1L,
                relative_path, " alias-resolution row is not reverse-unique/Entrez-consistent")
    }
  }
  invisible(audit)
}

qa_validate_meta <- function(
    run_dir, manifest,
    rel = c(
      "results/os_validation_meta_inputs.csv",
      "results/os_validation_meta_analysis.csv",
      "results/os_validation_meta_leave_one_out.csv"
    ),
    label = "OS meta-analysis", identity = NULL, performance = NULL) {
  qa_assert(length(rel) == 3L, label, " artifact contract must contain three files")
  identity_fields <- character()
  if (!is.null(identity)) {
    qa_assert(!is.null(names(identity)) && !anyDuplicated(names(identity)) &&
                all(nzchar(names(identity))), label, " identity must be named")
    identity_fields <- names(identity)
  }
  qa_require_manifested(manifest, rel)
  input_required <- c("cohort", "endpoint", "n", "events", "yi", "sei")
  if (!is.null(performance)) {
    input_required <- c(
      input_required, "mapping_method", "HR_per_SD", "lower95", "upper95", "cox_p"
    )
  }
  inputs <- qa_read_csv(run_dir, rel[[1L]], c(input_required, identity_fields))
  cohorts <- sort(qa_character(inputs$cohort, "OS meta input cohort"))
  expected_cohorts <- sort(c("TCGA-COAD", "GSE17536", "GSE17537"))
  qa_assert(identical(cohorts, expected_cohorts) && !anyDuplicated(inputs$cohort),
            "OS meta inputs are not the three prespecified external OS cohorts")
  qa_assert(all(qa_numeric(inputs$sei, "OS meta input SE") > 0),
            "OS meta input SE must be positive")
  qa_numeric(inputs$yi, "OS meta input log-HR")
  cohort_inputs <- paste(expected_cohorts, collapse = ";")
  required <- c(
    "endpoint", "cohorts", "pooled_HR", "lower95", "upper95", "p",
    "prediction_lower95", "prediction_upper95", "Q", "Q_df", "Q_p", "I2",
    "tau2", "test_statistic", "inference_df", "inference", "primary",
    "cohort_inputs", "analysis", identity_fields
  )
  meta <- qa_read_csv(run_dir, rel[[2L]], required)
  qa_assert(nrow(meta) == 2L, "Main OS meta-analysis must contain KH and Wald rows")
  primary <- qa_logical(meta$primary, "OS meta primary")
  qa_assert(sum(primary) == 1L &&
              identical(as.character(meta$inference[primary]), "Knapp-Hartung/t") &&
              identical(as.character(meta$inference[!primary]), "normal/Wald"),
            "Main OS meta-analysis does not use KH/t primary plus Wald sensitivity")
  qa_assert(all(qa_integer(meta$cohorts, "OS meta cohort count") == 3L) &&
              all(qa_integer(meta$Q_df, "OS meta Q_df") == 2L) &&
              qa_numeric(meta$inference_df[primary], "OS meta KH df") == 2 &&
              all(as.character(meta$cohort_inputs) == cohort_inputs),
            "Main OS meta-analysis cohort/df audit fields are inconsistent")
  for (field in c("Q", "Q_p", "tau2", "I2")) qa_numeric(meta[[field]], paste("OS meta", field))
  qa_assert(all(meta$Q >= 0) && all(meta$Q_p >= 0 & meta$Q_p <= 1) &&
              all(meta$tau2 >= 0) && all(meta$I2 >= 0),
            "Main OS meta heterogeneity fields are outside valid ranges")

  loo <- qa_read_csv(run_dir, rel[[3L]], c("dropped_cohort", required))
  qa_assert(nrow(loo) == 2L * length(expected_cohorts),
            "Leave-one-out OS meta-analysis does not have two rows per omission")
  qa_assert(setequal(unique(loo$dropped_cohort), expected_cohorts),
            "Leave-one-out omissions differ from OS meta inputs")
  for (drop in expected_cohorts) {
    z <- loo[loo$dropped_cohort == drop, , drop = FALSE]
    z_primary <- qa_logical(z$primary, paste("LOO primary", drop))
    expected_inputs <- paste(setdiff(expected_cohorts, drop), collapse = ";")
    qa_assert(nrow(z) == 2L && sum(z_primary) == 1L &&
                identical(as.character(z$inference[z_primary]), "Knapp-Hartung/t") &&
                identical(as.character(z$inference[!z_primary]), "normal/Wald") &&
                all(qa_integer(z$cohorts, "LOO cohort count") == 2L) &&
                all(qa_integer(z$Q_df, "LOO Q_df") == 1L) &&
                qa_numeric(z$inference_df[z_primary], "LOO KH df") == 1 &&
                all(as.character(z$cohort_inputs) == expected_inputs),
              "Leave-one-out dual inference is inconsistent for ", drop)
    for (field in c("Q", "Q_p", "tau2", "I2")) {
      qa_numeric(z[[field]], paste("LOO", drop, field))
    }
  }
  if (length(identity_fields)) {
    for (tab in list(inputs = inputs, meta = meta, leave_one_out = loo)) {
      for (field in identity_fields) {
        qa_assert(all(as.character(tab[[field]]) == as.character(identity[[field]])),
                  label, " ", field, " differs from the certified strict model")
      }
    }
  }
  if (!is.null(performance)) {
    idx <- match(as.character(inputs$cohort), as.character(performance$cohort))
    qa_assert(!anyNA(idx), label, " inputs include a cohort absent from performance")
    matched <- performance[idx, , drop = FALSE]
    qa_assert(all(as.character(inputs$endpoint) == as.character(matched$endpoint)) &&
                all(as.character(inputs$mapping_method) ==
                      as.character(matched$mapping_method)) &&
                all(qa_integer(inputs$n, paste(label, "input n"), 1L) ==
                      qa_integer(matched$n, paste(label, "performance n"), 1L)) &&
                all(qa_integer(inputs$events, paste(label, "input events"), 1L) ==
                      qa_integer(matched$events, paste(label, "performance events"), 1L)),
              label, " inputs disagree with cohort performance identities/counts")
    for (field in c("HR_per_SD", "lower95", "upper95", "cox_p")) {
      qa_assert(qa_all_close(
        qa_numeric(inputs[[field]], paste(label, "input", field)),
        qa_numeric(matched[[field]], paste(label, "performance", field)),
        tolerance = 1e-9
      ), label, " input ", field, " disagrees with cohort performance")
    }
    expected_yi <- log(qa_numeric(inputs$HR_per_SD, paste(label, "input HR")))
    expected_sei <- (
      log(qa_numeric(inputs$upper95, paste(label, "input upper95"))) -
        log(qa_numeric(inputs$lower95, paste(label, "input lower95")))
    ) / (2 * stats::qnorm(0.975))
    qa_assert(qa_all_close(inputs$yi, expected_yi, tolerance = 1e-9) &&
                qa_all_close(inputs$sei, expected_sei, tolerance = 1e-9),
              label, " yi/sei do not reconstruct from cohort effects")
  }
  invisible(TRUE)
}

qa_validate_clinical_and_diagnostics <- function(run_dir, manifest, params) {
  rel <- c(
    "results/clinical_adjustment_eligibility.csv",
    "results/proportional_hazards_test.csv",
    "results/risk_score_proportional_hazards_test.csv",
    "results/clinical_model_proportional_hazards_test.csv",
    "results/risk_score_linearity_test.csv"
  )
  qa_require_manifested(manifest, rel)
  eligibility <- qa_read_csv(run_dir, rel[[1L]], c(
    "cohort", "endpoint", "mapping_method", "source_n", "source_events",
    "complete_case_n", "complete_case_events", "excluded_from_adjustment_n",
    "excluded_events", "minimum_n_required", "minimum_events_required",
    "eligible", "reason"
  ))
  expected <- c(
    "GSE39582 (training)", "TCGA-COAD", "GSE14333", "GSE17536", "GSE17537"
  )
  qa_assert(nrow(eligibility) == 5L && setequal(eligibility$cohort, expected) &&
              !anyDuplicated(eligibility$cohort),
            "Clinical eligibility does not contain exactly the five prespecified cohorts")
  eligible <- qa_logical(eligibility$eligible, "Clinical eligible")
  reason <- qa_character(eligibility$reason, "Clinical eligibility reason")
  qa_assert(all(reason[eligible] == "eligible") && all(reason[!eligible] != "eligible"),
            "Clinical eligibility reasons do not preserve inclusion/exclusion decisions")
  source_n <- qa_integer(eligibility$source_n, "Clinical source_n", 0L)
  source_events <- qa_integer(eligibility$source_events, "Clinical source_events", 0L)
  complete_n <- qa_integer(eligibility$complete_case_n, "Clinical complete_case_n", 0L)
  complete_events <- qa_integer(
    eligibility$complete_case_events, "Clinical complete_case_events", 0L
  )
  qa_assert(all(source_events <= source_n) && all(complete_events <= complete_n) &&
              all(complete_n <= source_n) &&
              all(qa_integer(eligibility$excluded_from_adjustment_n,
                             "Clinical excluded n", 0L) == source_n - complete_n) &&
              all(qa_integer(eligibility$excluded_events,
                             "Clinical excluded events", 0L) == source_events - complete_events),
            "Clinical eligibility counts do not reconcile")
  min_n <- qa_integer(params$survival$clinical_min_n, "clinical_min_n", 1L)
  min_events <- qa_integer(params$survival$clinical_min_events, "clinical_min_events", 1L)
  qa_assert(all(qa_integer(eligibility$minimum_n_required,
                           "Clinical minimum_n_required") == min_n) &&
              all(qa_integer(eligibility$minimum_events_required,
                             "Clinical minimum_events_required") == min_events),
            "Clinical eligibility thresholds differ from run parameters")
  qa_assert(all(complete_n[eligible] >= min_n) &&
              all(complete_events[eligible] >= min_events),
            "An eligible clinical cohort is below a prespecified threshold")

  qa_read_csv(run_dir, rel[[2L]], c("term", "chisq", "df", "p"))
  qa_read_csv(run_dir, rel[[3L]], c(
    "cohort", "endpoint", "mapping_method", "n", "events", "chisq", "p"
  ))
  qa_read_csv(run_dir, rel[[4L]], c("cohort", "endpoint", "term", "chisq", "p"))
  linearity <- qa_read_csv(run_dir, rel[[5L]], c(
    "cohort", "endpoint", "n", "events", "linear_loglik", "spline_loglik",
    "p_nonlinearity", "method"
  ))
  qa_assert(all(nzchar(trimws(as.character(linearity$method)))),
            "Nonlinearity diagnostic lacks named methods")
  invisible(TRUE)
}

qa_validate_probe_frequency <- function(tab, denominator, analysis_key, label) {
  required <- c(
    "SYMBOL", "PROBEID", "selected_splits", "denominator", "frequency", "analysis_key"
  )
  missing <- setdiff(required, names(tab))
  if (length(missing)) qa_fail(label, " lacks columns: ", paste(missing, collapse = ", "))
  tab$SYMBOL <- qa_character(tab$SYMBOL, paste(label, "SYMBOL"))
  tab$PROBEID <- qa_character(tab$PROBEID, paste(label, "PROBEID"))
  qa_assert(!anyDuplicated(paste(tab$SYMBOL, tab$PROBEID, sep = "\r")),
            label, " has duplicated SYMBOL/PROBEID rows")
  selected <- qa_integer(tab$selected_splits, paste(label, "selected_splits"), 1L)
  den <- qa_integer(tab$denominator, paste(label, "denominator"), 1L)
  frequency <- qa_numeric(tab$frequency, paste(label, "frequency"))
  qa_assert(all(den == denominator) && all(selected <= denominator) &&
              all(abs(frequency - selected / denominator) <= 1e-12) &&
              all(as.character(tab$analysis_key) == analysis_key),
            label, " has invalid denominator/frequency/key")
  totals <- stats::aggregate(selected, list(SYMBOL = tab$SYMBOL), sum)
  qa_assert(all(totals$x == denominator),
            label, " does not conserve one selected probe per gene per split")
  invisible(tab)
}

qa_recompute_external_performance <- function(scores, ties) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    qa_fail("Package 'survival' is required to validate external performance")
  }
  rows <- lapply(split(scores, as.character(scores$cohort)), function(z) {
    z <- as.data.frame(z, stringsAsFactors = FALSE)
    qa_assert(length(unique(z$endpoint)) == 1L &&
                length(unique(z$mapping_method)) == 1L,
              "Highest-mean score cohort has mixed endpoint or mapping labels")
    z$group <- factor(as.character(z$group), levels = c("Low", "High"))
    qa_assert(!anyNA(z$group) && nlevels(droplevels(z$group)) == 2L,
              "Highest-mean score cohort lacks both Low and High groups")
    fit <- survival::coxph(
      survival::Surv(time, status) ~ risk_sd,
      data = z, ties = ties, x = TRUE, y = TRUE, na.action = stats::na.fail
    )
    fit_summary <- summary(fit)
    ci <- stats::confint(fit)
    km <- survival::survdiff(survival::Surv(time, status) ~ group, data = z)
    c_obj <- survival::concordance(
      survival::Surv(time, status) ~ risk_sd,
      data = z, reverse = TRUE, timewt = "n"
    )
    data.frame(
      cohort = unique(as.character(z$cohort)),
      endpoint = unique(as.character(z$endpoint)),
      mapping_method = unique(as.character(z$mapping_method)),
      n = nrow(z), events = sum(z$status),
      HR_per_SD = unname(exp(stats::coef(fit))),
      lower95 = unname(exp(ci[1L])), upper95 = unname(exp(ci[2L])),
      cox_p = fit_summary$coefficients[1L, "Pr(>|z|)"],
      logrank_p = stats::pchisq(km$chisq, df = 1, lower.tail = FALSE),
      c_index = unname(c_obj$concordance),
      c_index_method = "Harrell C; timewt=n; higher score=higher hazard",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

qa_validate_highest_mean <- function(run_dir, manifest, params, nested,
                                     annotation_version) {
  files <- c(
    outer_maps = "results/highest_mean_outer_probe_maps.csv",
    inner_hashes = "results/highest_mean_inner_probe_map_hashes.csv",
    outer_frequency = "results/highest_mean_outer_probe_frequency.csv",
    inner_frequency = "results/highest_mean_inner_probe_frequency.csv",
    full_map = "results/highest_mean_full_development_probe_map.csv",
    full_model = "results/highest_mean_full_development_model.rds",
    coefficients = "results/highest_mean_final_model_coefficients.csv",
    external_scores = "results/highest_mean_external_scores.csv",
    external_performance = "results/highest_mean_external_cohort_performance.csv",
    tcga_mapping = "results/highest_mean_TCGA_COAD_signature_gene_mapping.csv",
    meta_inputs = "results/highest_mean_os_validation_meta_inputs.csv",
    meta_analysis = "results/highest_mean_os_validation_meta_analysis.csv",
    meta_loo = "results/highest_mean_os_validation_meta_leave_one_out.csv"
  )
  qa_require_manifested(manifest, unname(files))
  analysis_key <- nested$analysis_key
  total_outer <- nested$repeats * nested$outer_folds
  total_inner <- total_outer * nested$inner_folds

  outer <- qa_read_csv(run_dir, files[["outer_maps"]], c(
    "repeat_id", "outer_fold", "SYMBOL", "PROBEID", "training_probe_mean",
    "training_probe_mean_hex",
    "n_probes_gene", "map_hash", "training_sample_hash", "training_n", "analysis_key"
  ))
  outer$repeat_id <- qa_integer(outer$repeat_id, "Highest-mean outer repeat_id", 1L)
  outer$outer_fold <- qa_integer(outer$outer_fold, "Highest-mean outer fold", 1L)
  qa_assert(all(as.character(outer$analysis_key) == analysis_key),
            "Highest-mean outer maps have the wrong analysis_key")
  groups <- unique(paste(outer$repeat_id, outer$outer_fold, sep = "\r"))
  qa_assert(length(groups) == total_outer, "Highest-mean outer maps omit model splits")
  for (repeat_id in seq_len(nested$repeats)) {
    assignment <- nested$outer_assignments[
      nested$outer_assignments$repeat_id == repeat_id, , drop = FALSE
    ]
    for (outer_fold in seq_len(nested$outer_folds)) {
      z <- outer[
        outer$repeat_id == repeat_id & outer$outer_fold == outer_fold, , drop = FALSE
      ]
      qa_assert(nrow(z) > 0L && !anyDuplicated(z$SYMBOL) &&
                  length(unique(z$map_hash)) == 1L &&
                  length(unique(z$training_sample_hash)) == 1L &&
                  length(unique(z$training_n)) == 1L,
                "Highest-mean outer map is not unique within a split")
      map_hash <- qa_hash_probe_map(z)
      train_samples <- setdiff(
        nested$sample_ids, assignment$sample[assignment$outer_fold == outer_fold]
      )
      qa_assert(identical(as.character(unique(z$map_hash)), map_hash) &&
                  unique(qa_integer(z$training_n, "Highest-mean outer training_n")) ==
                    length(train_samples) &&
                  identical(as.character(unique(z$training_sample_hash)),
                            qa_hash_sample_ids(train_samples)),
                "Highest-mean outer map hash/training identity mismatch")
    }
  }

  inner <- qa_read_csv(run_dir, files[["inner_hashes"]], c(
    "repeat_id", "outer_fold", "inner_fold", "probe_map_hash",
    "training_sample_hash", "training_n", "analysis_key"
  ))
  inner$repeat_id <- qa_integer(inner$repeat_id, "Highest-mean inner repeat_id", 1L)
  inner$outer_fold <- qa_integer(inner$outer_fold, "Highest-mean inner outer_fold", 1L)
  inner$inner_fold <- qa_integer(inner$inner_fold, "Highest-mean inner fold", 1L)
  qa_assert(nrow(inner) == total_inner &&
              !anyDuplicated(paste(
                inner$repeat_id, inner$outer_fold, inner$inner_fold, sep = "\r"
              )) && all(grepl(.qa_sha256_pattern, as.character(inner$probe_map_hash))) &&
              all(as.character(inner$analysis_key) == analysis_key),
            "Highest-mean inner map-hash table is incomplete")
  for (i in seq_len(nrow(inner))) {
    row <- inner[i, , drop = FALSE]
    assignments <- nested$inner_assignments[
      nested$inner_assignments$repeat_id == row$repeat_id &
        nested$inner_assignments$outer_fold == row$outer_fold, , drop = FALSE
    ]
    train_samples <- assignments$sample[assignments$inner_fold != row$inner_fold]
    qa_assert(qa_integer(row$training_n, "Highest-mean inner training_n") ==
                  length(train_samples) &&
                identical(as.character(row$training_sample_hash),
                          qa_hash_sample_ids(train_samples)),
              "Highest-mean inner map training identity mismatch")
  }

  outer_frequency <- qa_read_csv(run_dir, files[["outer_frequency"]])
  inner_frequency <- qa_read_csv(run_dir, files[["inner_frequency"]])
  qa_validate_probe_frequency(
    outer_frequency, total_outer, analysis_key, "Highest-mean outer frequency"
  )
  qa_validate_probe_frequency(
    inner_frequency, total_inner, analysis_key, "Highest-mean inner frequency"
  )
  observed <- stats::aggregate(
    rep.int(1L, nrow(outer)),
    list(SYMBOL = outer$SYMBOL, PROBEID = outer$PROBEID), sum
  )
  names(observed)[[3L]] <- "selected_splits"
  freq_key <- paste(outer_frequency$SYMBOL, outer_frequency$PROBEID, sep = "\r")
  observed_key <- paste(observed$SYMBOL, observed$PROBEID, sep = "\r")
  idx <- match(freq_key, observed_key)
  qa_assert(!anyNA(idx) && nrow(observed) == nrow(outer_frequency) &&
              all(qa_integer(outer_frequency$selected_splits,
                             "Highest-mean outer selected_splits") ==
                    observed$selected_splits[idx]),
            "Highest-mean outer frequencies do not reconstruct from outer maps")

  full_map <- qa_read_csv(run_dir, files[["full_map"]], c(
    "SYMBOL", "PROBEID", "training_probe_mean", "training_probe_mean_hex",
    "n_probes_gene", "map_hash",
    "training_sample_hash", "training_n", "analysis_key"
  ))
  full_hash <- qa_hash_probe_map(full_map)
  qa_assert(length(unique(full_map$map_hash)) == 1L &&
              identical(as.character(unique(full_map$map_hash)), full_hash) &&
              length(unique(full_map$training_sample_hash)) == 1L &&
              identical(as.character(unique(full_map$training_sample_hash)),
                        qa_hash_sample_ids(nested$sample_ids)) &&
              length(unique(full_map$training_n)) == 1L &&
              qa_integer(unique(full_map$training_n),
                         "Highest-mean full-development training_n") ==
                length(nested$sample_ids) &&
              all(as.character(full_map$analysis_key) == analysis_key),
            "Highest-mean full-development map hash/training identity mismatch")

  model_path <- qa_run_path(run_dir, files[["full_model"]])
  model <- tryCatch(readRDS(model_path), error = function(e) {
    qa_fail("Cannot read highest-mean full-development model: ", conditionMessage(e))
  })
  qa_assert(is.list(model) && all(c(
    "probe_map", "probe_map_hash", "gene_model_fingerprint", "model_fingerprint"
  ) %in% names(model)), "Highest-mean full-development model lacks frozen-map metadata")
  model_map <- qa_canonical_probe_map(model$probe_map$mapping)
  csv_map <- qa_canonical_probe_map(full_map)
  qa_assert(isTRUE(all.equal(model_map, csv_map, check.attributes = FALSE, tolerance = 0)),
            "Full-development map CSV differs from model RDS")
  qa_assert(identical(as.character(model$probe_map_hash), full_hash) &&
              identical(as.character(model$probe_map$map_hash), full_hash) &&
              identical(as.character(model$probe_map$training_sample_hash),
                        qa_hash_sample_ids(nested$sample_ids)) &&
              as.integer(model$probe_map$training_n) == length(nested$sample_ids),
            "Full-development model does not carry the certified map/training set")
  gene_hash <- qa_hash_gene_model(model)
  expected_model_hash <- qa_combine_model_fingerprint(gene_hash, full_hash)
  qa_assert(identical(as.character(model$gene_model_fingerprint), gene_hash) &&
              identical(as.character(model$model_fingerprint), expected_model_hash),
            "Full-development model fingerprint is invalid")

  coefficients <- qa_read_csv(run_dir, files[["coefficients"]], c(
    "gene", "coef", "probe_map_hash", "model_fingerprint", "analysis_key"
  ))
  coefficients$gene <- qa_character(coefficients$gene, "Highest-mean coefficient gene")
  qa_assert(identical(coefficients$gene, as.character(model$genes)) &&
              isTRUE(all.equal(
                qa_numeric(coefficients$coef, "Highest-mean coefficient"),
                as.numeric(model$coefs[model$genes]), tolerance = 1e-12,
                check.attributes = FALSE
              )) &&
              all(as.character(coefficients$probe_map_hash) == full_hash) &&
              all(as.character(coefficients$model_fingerprint) == expected_model_hash) &&
              all(as.character(coefficients$analysis_key) == analysis_key),
            "Highest-mean coefficient table is not bound to the certified model/map")

  scores <- qa_read_csv(run_dir, files[["external_scores"]], c(
    "cohort", "sample", "time", "status", "endpoint", "lp_raw", "risk_score",
    "risk_sd", "group", "mapping_method", "probe_map_hash", "model_fingerprint",
    "analysis_key"
  ))
  expected_cohorts <- c("TCGA-COAD", "GSE14333", "GSE17536", "GSE17537")
  qa_assert(setequal(unique(scores$cohort), expected_cohorts),
            "Highest-mean external scores do not cover all four external cohorts")
  qa_assert(!anyDuplicated(paste(scores$cohort, scores$sample, sep = "\r")),
            "Highest-mean external score sample IDs are duplicated within cohort")
  time <- qa_numeric(scores$time, "Highest-mean external time")
  status <- qa_integer(scores$status, "Highest-mean external status", 0L)
  lp_raw <- qa_numeric(scores$lp_raw, "Highest-mean external lp_raw")
  risk_score <- qa_numeric(scores$risk_score, "Highest-mean external risk_score")
  risk_sd <- qa_numeric(scores$risk_sd, "Highest-mean external risk_sd")
  scores$time <- time
  scores$status <- status
  scores$lp_raw <- lp_raw
  scores$risk_score <- risk_score
  scores$risk_sd <- risk_sd
  qa_assert(all(time > 0), "Highest-mean external time must be positive")
  qa_assert(all(status %in% 0:1), "Highest-mean external status must be 0 or 1")
  qa_assert(isTRUE(all.equal(
    lp_raw, risk_score, tolerance = 1e-12, check.attributes = FALSE
  )), "Highest-mean external lp_raw and risk_score differ")
  score_groups <- split(seq_len(nrow(scores)), as.character(scores$cohort))
  for (cohort in names(score_groups)) {
    idx <- score_groups[[cohort]]
    qa_assert(length(idx) >= 2L &&
                abs(mean(risk_sd[idx])) <= 1e-10 &&
                abs(stats::sd(risk_sd[idx]) - 1) <= 1e-10 &&
                is.finite(stats::sd(lp_raw[idx])) && stats::sd(lp_raw[idx]) > 0,
              "Highest-mean external score scaling is invalid for ", cohort)
  }
  is_tcga <- scores$cohort == "TCGA-COAD"
  qa_assert(all(as.character(scores$mapping_method[is_tcga]) == "gene_level_input") &&
              all(as.character(scores$mapping_method[!is_tcga]) ==
                    "full_development_frozen_unique_highest_mean") &&
              all(as.character(scores$probe_map_hash) == full_hash) &&
              all(as.character(scores$model_fingerprint) == expected_model_hash) &&
              all(as.character(scores$analysis_key) == analysis_key),
            "External scores are not bound to the same frozen full-development model/map")

  performance <- qa_read_csv(run_dir, files[["external_performance"]], c(
    "cohort", "endpoint", "mapping_method", "n", "events", "HR_per_SD",
    "lower95", "upper95", "cox_p", "logrank_p", "c_index", "c_index_method",
    "probe_map_hash", "model_fingerprint", "analysis_key"
  ))
  qa_assert(nrow(performance) == length(expected_cohorts) &&
              setequal(performance$cohort, expected_cohorts) &&
              !anyDuplicated(performance$cohort),
            "Highest-mean external performance must contain one row per cohort")
  expected_performance <- qa_recompute_external_performance(
    scores, ties = as.character(params$survival$cox_ties)
  )
  perf_idx <- match(as.character(performance$cohort), expected_performance$cohort)
  qa_assert(!anyNA(perf_idx), "Highest-mean performance cohort matching failed")
  expected_performance <- expected_performance[perf_idx, , drop = FALSE]
  qa_assert(all(as.character(performance$endpoint) == expected_performance$endpoint) &&
              all(as.character(performance$mapping_method) ==
                    expected_performance$mapping_method) &&
              all(qa_integer(performance$n, "Highest-mean performance n", 1L) ==
                    expected_performance$n) &&
              all(qa_integer(performance$events, "Highest-mean performance events", 1L) ==
                    expected_performance$events) &&
              all(as.character(performance$c_index_method) ==
                    expected_performance$c_index_method),
            "Highest-mean performance labels/counts differ from external scores")
  for (field in c(
      "HR_per_SD", "lower95", "upper95", "cox_p", "logrank_p", "c_index")) {
    qa_assert(qa_all_close(
      qa_numeric(performance[[field]], paste("Highest-mean performance", field)),
      expected_performance[[field]], tolerance = 1e-8
    ), "Highest-mean performance ", field, " does not reconstruct from scores")
  }
  qa_assert(all(as.character(performance$probe_map_hash) == full_hash) &&
              all(as.character(performance$model_fingerprint) == expected_model_hash) &&
              all(as.character(performance$analysis_key) == analysis_key),
            "Highest-mean performance is not bound to the certified model/map")

  tcga_audit <- qa_validate_mapping_audit(
    run_dir, manifest, files[["tcga_mapping"]], files[["coefficients"]],
    annotation_version
  )
  identity_fields <- c("probe_map_hash", "model_fingerprint", "analysis_key")
  missing_identity <- setdiff(identity_fields, names(tcga_audit))
  if (length(missing_identity)) {
    qa_fail(
      files[["tcga_mapping"]], " lacks certified-model identity fields: ",
      paste(missing_identity, collapse = ", ")
    )
  }
  expected_identity <- c(
    probe_map_hash = full_hash,
    model_fingerprint = expected_model_hash,
    analysis_key = analysis_key
  )
  for (field in identity_fields) {
    qa_assert(all(as.character(tcga_audit[[field]]) == expected_identity[[field]]),
              "Highest-mean TCGA mapping audit has the wrong ", field)
  }
  qa_validate_meta(
    run_dir, manifest,
    rel = unname(files[c("meta_inputs", "meta_analysis", "meta_loo")]),
    label = "Strict highest-mean OS meta-analysis",
    identity = expected_identity,
    performance = performance
  )
  invisible(TRUE)
}

validate_completed_run <- function(run_dir, quiet = FALSE) {
  qa_assert(length(run_dir) == 1L && !is.na(run_dir) && nzchar(run_dir),
            "Provide exactly one completed run directory")
  qa_assert(dir.exists(run_dir), "Run directory does not exist: ", run_dir)
  # Do not call normalizePath() here on Windows: it resolves the project's
  # ASCII junction to the Chinese target path, which R cannot reopen under the
  # active C locale. Keeping the caller's existing junction path is deliberate.
  run_dir <- sub("/+$", "", gsub("\\\\", "/", path.expand(run_dir)))
  directory_key <- basename(run_dir)
  qa_assert(grepl(.qa_sha256_pattern, directory_key),
            "Completed run directory name is not a SHA-256 run_key")
  qa_assert(!endsWith(directory_key, ".tmp"), "A temporary run cannot be certified")

  status <- qa_read_csv(run_dir, "results/run_status.csv", c(
    "status", "run_key", "nested_analysis_key", "input_manifest_sha256",
    "parameter_manifest_sha256", "code_manifest_sha256",
    "output_manifest_sha256", "finished_at_utc"
  ))
  qa_assert(nrow(status) == 1L && identical(as.character(status$status[[1L]]), "complete"),
            "run_status.csv is not exactly one complete row")
  qa_assert(identical(as.character(status$run_key[[1L]]), directory_key),
            "run_status run_key differs from the directory name")
  qa_assert(!is.na(status$finished_at_utc[[1L]]) &&
              nzchar(trimws(as.character(status$finished_at_utc[[1L]]))),
            "Complete status lacks finished_at_utc")

  manifest <- qa_validate_manifest(run_dir, status)
  qa_require_manifested(manifest, c(
    "results/input_manifest.csv", "results/code_manifest.csv",
    "results/run_parameters.txt"
  ))
  annotation_version <- qa_validate_code_snapshot(run_dir, manifest)
  params <- qa_validate_parameters(run_dir)
  primary_nested <- qa_validate_nested(run_dir, manifest, status, params)
  strict_nested <- qa_validate_strict_nested(
    run_dir, manifest, params, primary_nested$analysis_key
  )
  qa_validate_mapping_audit(
    run_dir, manifest, "results/TCGA_COAD_signature_gene_mapping.csv",
    "results/final_model_coefficients.csv", annotation_version
  )
  qa_validate_meta(run_dir, manifest)
  qa_validate_clinical_and_diagnostics(run_dir, manifest, params)
  qa_validate_highest_mean(
    run_dir, manifest, params, strict_nested, annotation_version
  )

  if (!quiet) {
    cat("COMPLETED-RUN QA: PASS\n")
    cat("run_key:", directory_key, "\n")
    cat("pipeline_version:", as.character(params$pipeline_version), "\n")
    cat("development_samples:", length(primary_nested$sample_ids), "\n")
    cat("strict_analysis_key:", strict_nested$analysis_key, "\n")
  }
  invisible(list(
    run_key = directory_key,
    pipeline_version = as.character(params$pipeline_version),
    development_samples = length(primary_nested$sample_ids),
    strict_analysis_key = strict_nested$analysis_key
  ))
}

qa_cli_main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L || !nzchar(args[[1L]])) {
    cat("Usage: Rscript validate_completed_run.R <completed-run-directory>\n", file = stderr())
    quit(save = "no", status = 2L)
  }
  error <- NULL
  tryCatch(
    validate_completed_run(args[[1L]], quiet = FALSE),
    error = function(e) error <<- conditionMessage(e)
  )
  if (!is.null(error)) {
    cat("COMPLETED-RUN QA: FAIL\n", file = stderr())
    cat(error, "\n", file = stderr())
    quit(save = "no", status = 1L)
  }
  quit(save = "no", status = 0L)
}

if (sys.nframe() == 0L) qa_cli_main()
