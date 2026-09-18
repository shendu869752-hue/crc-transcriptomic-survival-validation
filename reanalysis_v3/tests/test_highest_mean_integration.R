#!/usr/bin/env Rscript

# Static production-integration contract for the strict highest-mean
# sensitivity. The core script is parsed but never sourced or executed.

options(stringsAsFactors = FALSE, warn = 1)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript")
this_file <- gsub("\\\\", "/", sub("^--file=", "", script_arg[[1L]]))
test_dir <- dirname(this_file)
v3_root <- dirname(test_dir)
core_file <- file.path(v3_root, "scripts", "01_core_reanalysis_v3.R")
highest_file <- file.path(v3_root, "scripts", "highest_mean_sensitivity.R")
if (!file.exists(core_file) || !file.exists(highest_file)) {
  stop("Cannot locate production integration files")
}

core_lines <- readLines(core_file, warn = FALSE, encoding = "UTF-8")
core_text <- paste(core_lines, collapse = "\n")
results <- list()

assert_true <- function(value, message) {
  if (length(value) != 1L || is.na(value) || !isTRUE(value)) stop(message)
  invisible(TRUE)
}

assert_contains <- function(pattern, message, fixed = TRUE) {
  assert_true(grepl(pattern, core_text, fixed = fixed), message)
}

run_test <- function(name, expr) {
  test_expr <- substitute(expr)
  status <- "PASS"
  detail <- ""
  started <- proc.time()[["elapsed"]]
  tryCatch(
    eval(test_expr, envir = parent.frame()),
    error = function(e) {
      status <<- "FAIL"
      detail <<- conditionMessage(e)
    }
  )
  elapsed <- proc.time()[["elapsed"]] - started
  results[[length(results) + 1L]] <<- data.frame(
    test = name, status = status, seconds = round(elapsed, 3), detail = detail,
    stringsAsFactors = FALSE
  )
  cat(sprintf("[%s] %s%s\n", status, name,
              if (nzchar(detail)) paste0(" -- ", detail) else ""))
}

run_test("production scripts parse", {
  invisible(parse(file = core_file))
  invisible(parse(file = highest_file))
})

run_test("raw probe input and sensitivity code enter run identity", {
  assert_contains(
    "GSE39582_probe_expression = gse_probe_expression_path",
    "raw GSE39582 probe matrix is absent from required_inputs"
  )
  assert_contains(
    "highest_mean_sensitivity = file.path(",
    "highest_mean_sensitivity.R is absent from code_files"
  )
  assert_contains(
    "highest_mean_sensitivity.R\"\n))",
    "highest-mean implementation is not sourced"
  )
  assert_contains(
    "make_run_key(\n  input_paths = required_inputs,\n  code_paths = code_files",
    "run key is not derived from the expanded input/code sets"
  )
})

run_test("mapping scope is training learned and external frozen", {
  assert_contains(
    "training_partition_learned_full_development_frozen_external",
    "locked highest-mean scope is missing from version/parameters"
  )
  assert_true(
    !grepl("cohort_adaptive_highest_mean", core_text, fixed = TRUE),
    "legacy cohort-adaptive highest-mean claim remains in production core"
  )
})

run_test("strict nested result is rerun or fail-closed validated", {
  for (symbol in c(
      "run_probe_nested_cv(", "validate_probe_nested_result(",
      "highest_probe_annotation", "highest_analysis_identity$manifest",
      "highest_mean_nested_cv_")) {
    assert_contains(symbol, paste("missing strict nested integration:", symbol))
  }
})

run_test("all mandatory strict artifacts are allow-listed", {
  required <- c(
    "highest_mean_nested_cv_outer_fold_assignments.csv",
    "highest_mean_nested_cv_inner_fold_assignments.csv",
    "highest_mean_nested_cv_oof_predictions.csv",
    "highest_mean_nested_cv_inner_tuning_performance.csv",
    "highest_mean_outer_probe_maps.csv",
    "highest_mean_inner_probe_map_hashes.csv",
    "highest_mean_outer_probe_frequency.csv",
    "highest_mean_inner_probe_frequency.csv",
    "highest_mean_full_development_probe_map.csv",
    "highest_mean_full_development_model.rds",
    "highest_mean_final_model_coefficients.csv",
    "highest_mean_external_scores.csv",
    "highest_mean_external_cohort_performance.csv",
    "highest_mean_TCGA_COAD_signature_gene_mapping.csv",
    "highest_mean_os_validation_meta_inputs.csv",
    "highest_mean_os_validation_meta_analysis.csv",
    "highest_mean_os_validation_meta_leave_one_out.csv"
  )
  absent <- required[!vapply(
    required, function(x) grepl(paste0("\"", x, "\""), core_text, fixed = TRUE),
    logical(1)
  )]
  assert_true(!length(absent), paste(
    "missing mandatory artifact literals:", paste(absent, collapse = ", ")
  ))
})

run_test("external scoring freezes probes but preserves cohort scaling", {
  assert_contains(
    "full_development_frozen_unique_highest_mean",
    "external GPL570 sensitivity label is not explicit"
  )
  assert_contains(
    "apply_frozen_probe_map(\n    probe_mat, highest_full_model$probe_map",
    "GPL570 external scoring does not apply the frozen full-development map"
  )
  assert_contains(
    "highest_geo_data[[nm]] <- apply_signature(",
    "frozen GEO gene matrices do not enter apply_signature"
  )
  assert_contains(
    "highest_tcga_d <- apply_signature(",
    "TCGA sensitivity rows do not enter apply_signature"
  )
  assert_true(
    !grepl("predict_probe_signature_pipeline(", core_text, fixed = TRUE),
    "production external effects incorrectly use training-center prediction"
  )
  for (field in c(
      "highest_external_d$lp_raw",
      "highest_external_d$risk_score",
      "highest_external_d$probe_map_hash",
      "highest_external_d$model_fingerprint",
      "highest_external_d$analysis_key")) {
    assert_contains(field, paste("external score identity missing:", field))
  }
  assert_contains(
    '"TCGA-COAD",\n  "gene_level_input"',
    "TCGA strict sensitivity is not labeled as gene-level input"
  )
})

run_test("strict OS meta-analysis uses the shared REML dual-inference helpers", {
  assert_contains(
    "highest_meta_summary <- fit_reml_meta_dual(",
    "strict OS meta-analysis does not use fit_reml_meta_dual"
  )
  assert_contains(
    "highest_leave_one_out <- leave_one_out_reml_meta_dual(",
    "strict leave-one-out meta-analysis is missing"
  )
  assert_contains(
    "c(\"TCGA-COAD\", \"GSE17536\", \"GSE17537\")",
    "strict OS meta-analysis cohort set is not explicit"
  )
})

run_test("Fig4 compares primary and strict full pipelines", {
  assert_contains(
    "primary_unique_mean_full_pipeline",
    "Fig4/CSV lacks primary full-pipeline label"
  )
  assert_contains(
    "strict_highest_mean_full_pipeline",
    "Fig4/CSV lacks strict full-pipeline label"
  )
  assert_true(
    !grepl("Cohort-adaptive highest-mean probe", core_text, fixed = TRUE),
    "Fig4 retains the invalid cohort-adaptive sensitivity label"
  )
})

run_test("KM censor coordinates retain the corrected survival summary", {
  assert_contains(
    "summary(fit, censored = TRUE)",
    "KM curve no longer requests censored observations"
  )
})

summary_df <- do.call(rbind, results)
cat("\nHighest-mean production integration test summary\n")
print(summary_df, row.names = FALSE)
cat(sprintf("\nPASS=%d FAIL=%d\n",
            sum(summary_df$status == "PASS"),
            sum(summary_df$status == "FAIL")))
if (any(summary_df$status != "PASS")) quit(save = "no", status = 1L)
quit(save = "no", status = 0L)
