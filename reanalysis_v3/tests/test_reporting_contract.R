#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript")
this_file <- gsub("\\\\", "/", sub("^--file=", "", script_arg[[1L]]))
test_dir <- dirname(this_file)
v3_root <- dirname(test_dir)
work_root <- dirname(v3_root)
utils_file <- file.path(v3_root, "scripts", "utils.R")
core_file <- file.path(v3_root, "scripts", "01_core_reanalysis_v3.R")

Sys.setenv(CRC_WORK_ROOT = work_root)
if (!nzchar(Sys.getenv("CRC_SOURCE_ROOT"))) {
  Sys.setenv(CRC_SOURCE_ROOT = work_root)
}
source(utils_file)

results <- list()

assert_true <- function(value, message = "assertion failed") {
  if (length(value) != 1L || is.na(value) || !isTRUE(value)) stop(message)
  invisible(TRUE)
}

expect_error <- function(expr, pattern = NULL) {
  captured <- NULL
  tryCatch(force(expr), error = function(e) captured <<- e)
  if (is.null(captured)) stop("expected an error, but expression returned normally")
  if (!is.null(pattern) && !grepl(pattern, conditionMessage(captured), ignore.case = TRUE)) {
    stop("error did not match /", pattern, "/: ", conditionMessage(captured))
  }
  invisible(captured)
}

run_test <- function(name, code) {
  error <- NULL
  tryCatch(force(code), error = function(e) error <<- conditionMessage(e))
  if (is.null(error)) {
    results[[length(results) + 1L]] <<- data.frame(
      test = name, status = "PASS", detail = "", stringsAsFactors = FALSE
    )
    cat("PASS:", name, "\n")
  } else {
    results[[length(results) + 1L]] <<- data.frame(
      test = name, status = "FAIL", detail = error, stringsAsFactors = FALSE
    )
    cat("FAIL:", name, "-", error, "\n")
  }
  invisible(NULL)
}

run_test("core and utility scripts parse", {
  parse(file = utils_file)
  parse(file = core_file)
})

run_test("parameter manifest preserves every list name", {
  td <- tempfile("crc9-parameters-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE, force = TRUE), add = TRUE)
  parameters <- list(
    output_contract_version = "test-v1",
    nested_pipeline = list(
      pipeline_version = "synthetic-v1",
      outer_folds = 5L,
      lambda_ratio_grid = c(
        0.562341325190349, 0.431933279405154, 0.331767112784286,
        0.254829674797935, 0.195734178148766, 0.150343041978733,
        0.115478198468946, 0.0886985799018192, 0.0681292069057962,
        0.0523299114681495, 0.0401945033361513, 0.0308733199257026,
        0.0237137370566166, 0.0182144753639595, 0.0139905031413729,
        0.0107460782832132, 0.00825404185268019,
        0.00633991351172485, 0.00486967525165863,
        0.00374038810036779, 0.00287298483335367,
        0.00220673406908459, 0.00169498815139035,
        0.00130191710619008, 0.001
      )
    ),
    bootstrap = list(replicates = 100L, probabilities = c(0.025, 0.975))
  )
  path <- file.path(td, "run_parameters.txt")
  write_dput_manifest(parameters, path)
  restored <- dget(path)
  assert_true(identical(restored, parameters), "parameter round-trip changed values")
  assert_true(identical(names(restored), names(parameters)), "top-level names were lost")
  assert_true(
    identical(names(restored$nested_pipeline), names(parameters$nested_pipeline)),
    "nested parameter names were lost"
  )
})

run_test("code snapshot preserves exact bytes and hashes", {
  td <- tempfile("crc9-code-snapshot-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE, force = TRUE), add = TRUE)
  source_dir <- file.path(td, "source")
  run_dir <- file.path(td, "run")
  dir.create(source_dir)
  dir.create(run_dir)
  writeLines("alpha <- 1L", file.path(source_dir, "alpha.R"), useBytes = TRUE)
  writeLines("beta <- 2L", file.path(source_dir, "beta.R"), useBytes = TRUE)
  paths <- c(
    alpha = file.path(source_dir, "alpha.R"),
    beta = file.path(source_dir, "beta.R")
  )
  snapshot <- snapshot_code_files(paths, run_dir)
  assert_true(nrow(snapshot) == 2L, "snapshot row count differs from code inputs")
  assert_true(all(file.exists(file.path(run_dir, snapshot$relative_path))),
              "snapshot file is missing")
  source_hashes <- vapply(paths[sort(names(paths))], sha256_file, character(1))
  assert_true(identical(unname(source_hashes), snapshot$sha256),
              "snapshot hashes differ from the source hashes")
  expect_error(snapshot_code_files(paths, run_dir), "already exists")
})

run_test("dual REML meta-analysis reports KH primary and Wald sensitivity", {
  d <- data.frame(
    cohort = c("Cohort-C", "Cohort-A", "Cohort-B"),
    yi = log(c(1.48, 1.12, 1.31)),
    sei = c(0.13, 0.10, 0.12),
    stringsAsFactors = FALSE
  )
  out <- fit_reml_meta_dual(d, analysis = "synthetic contract test")
  required <- c(
    "pooled_HR", "lower95", "upper95", "p", "prediction_lower95",
    "prediction_upper95", "Q", "Q_df", "Q_p", "I2", "tau2",
    "test_statistic", "inference_df", "inference", "primary",
    "cohort_inputs", "analysis"
  )
  assert_true(nrow(out) == 2L, "dual inference did not return two rows")
  assert_true(all(required %in% names(out)), "meta summary lacks required audit fields")
  assert_true(sum(out$primary) == 1L, "meta summary does not have exactly one primary row")
  assert_true(out$inference[out$primary] == "Knapp-Hartung/t",
              "KH/t is not the primary inference")
  assert_true(out$inference[!out$primary] == "normal/Wald",
              "Wald sensitivity row is missing")
  assert_true(out$inference_df[out$primary] == 2,
              "KH degrees of freedom are not k-1")
  assert_true(is.na(out$inference_df[!out$primary]),
              "Wald row unexpectedly reports t degrees of freedom")
  assert_true(all(out$Q_df == 2L), "Cochran Q degrees of freedom are incorrect")
  assert_true(all(is.finite(out$Q)) && all(is.finite(out$Q_p)),
              "Cochran Q or Q-test p is not finite")
  assert_true(all(out$cohort_inputs == "Cohort-A;Cohort-B;Cohort-C"),
              "cohort input identifiers are not exact and deterministic")
  assert_true(out$lower95[out$primary] != out$lower95[!out$primary],
              "KH and Wald intervals were not computed separately")
})

run_test("leave-one-out retains KH primary and Wald sensitivity for every omission", {
  d <- data.frame(
    cohort = c("TCGA", "GSE17536", "GSE17537"),
    yi = log(c(1.28, 1.34, 1.21)),
    sei = c(0.10, 0.12, 0.11),
    stringsAsFactors = FALSE
  )
  out <- leave_one_out_reml_meta_dual(d)
  assert_true(nrow(out) == 6L, "three-cohort LOO did not return six rows")
  for (cohort in sort(d$cohort)) {
    z <- out[out$dropped_cohort == cohort, , drop = FALSE]
    assert_true(nrow(z) == 2L, paste("LOO row count failed for", cohort))
    assert_true(sum(z$primary) == 1L, paste("LOO primary count failed for", cohort))
    assert_true(z$inference[z$primary] == "Knapp-Hartung/t",
                paste("LOO KH primary failed for", cohort))
    assert_true(z$inference_df[z$primary] == 1,
                paste("LOO KH df failed for", cohort))
    assert_true(!grepl(cohort, z$cohort_inputs[[1L]], fixed = TRUE),
                paste("dropped cohort remains in cohort_inputs for", cohort))
  }
})

run_test("KM censor marks retain the actual censoring times", {
  d <- data.frame(
    time = c(1, 2, 3, 4, 1.5, 2.5, 3.5, 4.5),
    status = c(1, 0, 1, 0, 1, 0, 1, 0),
    group = factor(rep(c("Low", "High"), each = 4L),
                   levels = c("Low", "High"))
  )
  fit <- survival::survfit(survival::Surv(time, status) ~ group, data = d)
  event_only <- summary(fit)
  with_censor <- summary(fit, censored = TRUE)
  assert_true(length(with_censor$time) > length(event_only$time),
              "censored=TRUE did not retain the censor-only time points")
  assert_true(any(with_censor$n.censor > 0L),
              "survfit summary has no auditable censor markers")
  censor_times <- with_censor$time[with_censor$n.censor > 0L]
  assert_true(all(c(2, 4, 2.5, 4.5) %in% censor_times),
              "one or more actual censoring times were lost")
})

run_test("core declares the revised reporting and reproducibility artifacts", {
  core <- readLines(core_file, warn = FALSE)
  required_literals <- c(
    '"clinical_adjustment_eligibility.csv"',
    '"os_validation_meta_inputs.csv"',
    '"os_adjusted_meta_inputs.csv"',
    '"km_risk_table.csv"',
    'file.path("code_snapshot", basename(code_files))',
    "write_dput_manifest(run_parameters, parameter_manifest_path)",
    "leave_one_out_reml_meta_dual(",
    "summary(fit, censored = TRUE)"
  )
  for (literal in required_literals) {
    assert_true(any(grepl(literal, core, fixed = TRUE)),
                paste("core contract is missing", literal))
  }
  assert_true(!any(grepl("dput(run_parameters, control", core, fixed = TRUE)),
              "name-dropping dput control is still present")
  assert_true(!any(grepl("lambda_ratio_grid = 10^seq", core, fixed = TRUE)),
              "runtime-computed lambda grid can fail exact text round-trip")
  assert_true(any(grepl("0.562341325190349", core, fixed = TRUE)) &&
                any(grepl("0.001", core, fixed = TRUE)),
              "frozen decimal lambda-ratio grid is missing")
})

run_test("core declares the A-level manuscript figure contract", {
  core <- readLines(core_file, warn = FALSE)
  core_text <- paste(core, collapse = "\n")
  required_literals <- c(
    '"Fig1_nested_validation_stability.pdf"',
    '"Supplementary_FigS1_LASSO_diagnostic.pdf"',
    "grid::grid.layout(nrow = 2L, ncol = 2L)",
    "pooled_primary$prediction_lower95",
    "shape = 23",
    "KH 95%% prediction interval",
    "km_external_order",
    "km_x_breaks",
    "time_years = 0, survival = 1",
    "scale_linetype_manual",
    "summary(fit, censored = TRUE)",
    "identical(primary_identity$sample, strict_identity$sample)",
    "Strict highest-mean development pipeline"
  )
  for (literal in required_literals) {
    assert_true(grepl(literal, core_text, fixed = TRUE),
                paste("A-level figure contract is missing", literal))
  }
  assert_true(!grepl("Fig1_LASSO_diagnostic", core_text, fixed = TRUE),
              "LASSO diagnostic remains incorrectly registered as main Fig1")
  nested_write <- grep(
    'fwrite(nested_summary, file.path(RESULT_DIR, "nested_cv_summary.csv"))',
    core, fixed = TRUE
  )
  fig1_draw <- grep("draw_nested_validation_figure()", core, fixed = TRUE)
  assert_true(length(nested_write) == 1L && length(fig1_draw) == 2L &&
                max(fig1_draw) > nested_write,
              "main Fig1 is not generated after the primary nested result")
  assert_true(grepl(
    'save_pub_plot(p_forest, file.path(FIGURE_DIR, "Fig2_cohort_forest"), width = 7.2',
    core_text, fixed = TRUE
  ), "Fig2 exceeds the 7.2-inch manuscript width")
  assert_true(grepl(
    'file.path(FIGURE_DIR, "Fig3_multicohort_KM.pdf"),\n  width = 7.2',
    core_text, fixed = TRUE
  ), "Fig3 exceeds the 7.2-inch manuscript width")
})

report <- do.call(rbind, results)
cat("\nSummary:\n")
print(table(report$status))
if (any(report$status != "PASS")) quit(save = "no", status = 1L)
quit(save = "no", status = 0L)
