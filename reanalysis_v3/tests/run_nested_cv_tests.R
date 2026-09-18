#!/usr/bin/env Rscript

# Contract and regression tests for the strict CRC9 nested-CV implementation.
#
# The suite intentionally uses base-R assertions instead of testthat so it can
# run in a clean R installation containing only the analysis dependencies.  A
# missing required interface is BLOCKED (and yields a non-zero exit status),
# never silently counted as a pass.

options(stringsAsFactors = FALSE, warn = 1)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript")
# Keep a relative invocation relative.  R 4.6.1 on this Windows host starts
# with a C.UTF-8 locale warning and can corrupt non-ASCII components returned
# by normalizePath(); relative paths avoid that platform issue.
this_file <- gsub("\\\\", "/", sub("^--file=", "", script_arg[[1L]]))
test_dir <- dirname(this_file)
v3_root <- dirname(test_dir)
work_root <- dirname(v3_root)
utils_file <- file.path(v3_root, "scripts", "utils.R")
core_file <- file.path(v3_root, "scripts", "01_core_reanalysis_v3.R")

if (!file.exists(utils_file) || !file.exists(core_file)) {
  stop("Cannot locate v3 scripts relative to test file")
}

Sys.setenv(CRC_WORK_ROOT = work_root)
if (!nzchar(Sys.getenv("CRC_SOURCE_ROOT"))) {
  Sys.setenv(CRC_SOURCE_ROOT = work_root)
}

results <- list()
artifacts <- new.env(parent = emptyenv())
pipeline_env <- new.env(parent = globalenv())
source_error <- NULL

new_blocked <- function(message) {
  structure(list(message = message, call = NULL),
            class = c("blocked_test", "error", "condition"))
}

block_if_missing <- function(names) {
  absent <- names[!vapply(names, exists, logical(1), envir = pipeline_env,
                          mode = "function", inherits = FALSE)]
  if (length(absent)) {
    stop(new_blocked(paste0("required interface not implemented: ",
                            paste(absent, collapse = ", "))))
  }
  invisible(TRUE)
}

assert_true <- function(value, message = "assertion failed") {
  if (length(value) != 1L || is.na(value) || !isTRUE(value)) stop(message)
  invisible(TRUE)
}

assert_equal <- function(actual, expected, tolerance = 1e-10,
                         message = NULL) {
  ok <- isTRUE(all.equal(actual, expected, tolerance = tolerance,
                         check.attributes = TRUE))
  if (!ok) {
    if (is.null(message)) {
      message <- paste0("objects differ: ",
                        paste(all.equal(actual, expected, tolerance = tolerance,
                                        check.attributes = TRUE), collapse = "; "))
    }
    stop(message)
  }
  invisible(TRUE)
}

assert_set_equal <- function(actual, expected, message = "sets differ") {
  if (!setequal(actual, expected)) stop(message)
  invisible(TRUE)
}

expect_error <- function(expr, pattern = NULL) {
  captured <- NULL
  tryCatch(
    force(expr),
    error = function(e) captured <<- e
  )
  if (is.null(captured)) stop("expected an error, but expression returned normally")
  if (!is.null(pattern) && !grepl(pattern, conditionMessage(captured),
                                  ignore.case = TRUE)) {
    stop("error did not match /", pattern, "/: ", conditionMessage(captured))
  }
  invisible(captured)
}

run_test <- function(name, expr) {
  test_expr <- substitute(expr)
  test_env <- new.env(parent = parent.frame())
  started <- proc.time()[["elapsed"]]
  status <- "PASS"
  detail <- ""
  tryCatch(
    eval(test_expr, envir = test_env),
    blocked_test = function(e) {
      status <<- "BLOCKED"
      detail <<- conditionMessage(e)
    },
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
  invisible(status)
}

function_definitions <- function(node, target = NULL) {
  found <- list()
  walk <- function(x) {
    if (!is.call(x)) return(invisible(NULL))
    head <- if (is.symbol(x[[1L]])) as.character(x[[1L]]) else ""
    if (head %in% c("<-", "=", "assign") && length(x) >= 3L) {
      lhs <- x[[2L]]
      rhs <- x[[3L]]
      lhs_name <- if (is.symbol(lhs)) as.character(lhs) else NA_character_
      if (!is.na(lhs_name) && is.call(rhs) && identical(as.character(rhs[[1L]]), "function") &&
          (is.null(target) || identical(lhs_name, target))) {
        found[[length(found) + 1L]] <<- list(name = lhs_name, definition = rhs)
      }
    }
    for (i in seq_along(x)[-1L]) walk(x[[i]])
    invisible(NULL)
  }
  if (is.expression(node)) {
    for (x in node) walk(x)
  } else {
    walk(node)
  }
  found
}

called_symbols <- function(node) {
  out <- character()
  walk <- function(x) {
    if (!is.call(x)) return(invisible(NULL))
    head <- x[[1L]]
    if (is.symbol(head)) out <<- c(out, as.character(head))
    for (i in seq_along(x)[-1L]) walk(x[[i]])
    invisible(NULL)
  }
  walk(node)
  unique(out)
}

synthetic_params <- function() {
  list(
    pipeline_version = "synthetic-nested-v1",
    outer_folds = 3L,
    outer_repeats = 1L,
    inner_folds = 2L,
    base_seed = 20240916L,
    sd_cutoff = 0.05,
    univ_p_cutoff = 0.20,
    min_candidates = 3L,
    candidate_fallback_n = 12L,
    coef_cutoff = 0.05,
    min_selected = 2L,
    selected_fallback_n = 4L,
    max_selected = 6L,
    lambda_ratio_grid = c(0.50, 0.15, 0.05),
    lambda_rule = "one_se_larger_penalty",
    ties = "efron"
  )
}

make_synthetic_survival <- function(n = 72L, p = 24L, seed = 7001L) {
  set.seed(seed)
  expr <- matrix(stats::rnorm(p * n), nrow = p, ncol = n)
  rownames(expr) <- sprintf("G%03d", seq_len(p))
  sample_ids <- sprintf("S%03d", seq_len(n))
  colnames(expr) <- sample_ids
  linear_predictor <- 0.90 * expr[1L, ] - 0.70 * expr[2L, ] +
    0.45 * expr[3L, ]
  event_time <- stats::rexp(n, rate = exp(linear_predictor) / 12)
  censor_time <- stats::rexp(n, rate = 1 / 20)
  time <- pmax(pmin(event_time, censor_time) * 30.4375, 0.01)
  status <- as.integer(event_time <= censor_time)
  if (sum(status) < 12L || sum(status == 0L) < 12L) {
    stop("synthetic-data generator produced inadequate event balance")
  }
  list(expr = expr, time = time, status = status, sample_ids = sample_ids)
}

assignment_table <- function(x, sample_ids, fold_name, repeat_default = 1L) {
  if (is.atomic(x) && is.null(dim(x))) {
    if (length(x) != length(sample_ids)) stop("assignment vector has wrong length")
    return(data.frame(sample = sample_ids, repeat_id = repeat_default,
                      fold_value = as.integer(x), stringsAsFactors = FALSE))
  }
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  sample_col <- intersect(c("sample", "sample_id"), names(x))
  repeat_col <- intersect(c("repeat_id", "repeat"), names(x))
  fold_col <- intersect(c(fold_name, "fold", "fold_id"), names(x))
  if (!length(sample_col) || !length(fold_col)) {
    stop("assignment table must expose sample and fold columns")
  }
  data.frame(
    sample = as.character(x[[sample_col[[1L]]]]),
    repeat_id = if (length(repeat_col)) as.integer(x[[repeat_col[[1L]]]]) else repeat_default,
    fold_value = as.integer(x[[fold_col[[1L]]]]),
    stringsAsFactors = FALSE
  )
}

model_fingerprint <- function(model) {
  fp <- model$model_fingerprint
  if (is.null(fp) || length(fp) != 1L || !nzchar(fp)) {
    stop("fit_signature_pipeline() must return a scalar model_fingerprint")
  }
  as.character(fp)
}

key_value <- function(x) {
  if (is.character(x) && length(x) == 1L) return(x)
  if (is.list(x) && is.character(x$analysis_key) && length(x$analysis_key) == 1L) {
    return(x$analysis_key)
  }
  stop("make_analysis_key() must return a scalar key or list$analysis_key")
}

cat("Strict nested-CV tests\n")
cat("R:", R.version.string, "\n")
cat("utils:", utils_file, "\n")
cat("core:", core_file, "\n\n")

run_test("R parser accepts utils.R and 01_core_reanalysis_v3.R", {
  parse(file = utils_file, keep.source = TRUE)
  parse(file = core_file, keep.source = TRUE)
})

run_test("utils.R sources without reading project data", {
  tryCatch(
    sys.source(utils_file, envir = pipeline_env),
    error = function(e) {
      source_error <<- conditionMessage(e)
      stop(e)
    }
  )
})

required_functions <- c(
  "make_stratified_folds", "prepare_training", "select_refit_genes",
  "fit_at_ratio", "fit_signature_pipeline", "predict_signature_pipeline",
  "run_nested_cv", "summarize_repeated_oof", "make_analysis_key",
  "make_run_key"
)
run_test("strict public test interface is complete", {
  if (!is.null(source_error)) stop(new_blocked(paste("utils source failed:", source_error)))
  block_if_missing(required_functions)
  forbidden <- grep("test|valid|new", names(formals(get("fit_signature_pipeline",
                                                        envir = pipeline_env))),
                    ignore.case = TRUE, value = TRUE)
  assert_true(!length(forbidden), paste0(
    "training interface leaks validation/test inputs through formals: ",
    paste(forbidden, collapse = ", ")
  ))
})

run_test("stratified folds are deterministic, exhaustive, and class-balanced", {
  block_if_missing("make_stratified_folds")
  f <- get("make_stratified_folds", envir = pipeline_env)
  status <- rep(c(0L, 1L), each = 18L)
  fold_1 <- f(status, k = 3L, seed = 91L)
  fold_2 <- f(status, k = 3L, seed = 91L)
  assert_equal(fold_1, fold_2, message = "same seed did not reproduce folds")
  assert_set_equal(unique(fold_1), 1:3, "fold ids are incomplete")
  assert_true(all(tabulate(fold_1, nbins = 3L) > 0L), "an empty fold was generated")
  for (j in 1:3) {
    assert_set_equal(unique(status[fold_1 == j]), 0:1,
                     paste("fold", j, "does not contain both outcome classes"))
  }
  expect_error(f(c(rep(0L, 10L), 1L), k = 3L, seed = 91L))
})

run_test("feature-selection fallbacks are deterministic and fully specified", {
  block_if_missing(c("select_refit_genes", "prepare_training"))
  select_fun <- get("select_refit_genes", envir = pipeline_env)
  params <- synthetic_params()
  params$min_selected <- 5L
  params$selected_fallback_n <- 8L
  params$max_selected <- 30L
  genes <- sprintf("G%02d", seq_len(40L))
  candidate_rank <- rev(genes)

  beta_zero <- stats::setNames(rep(0, length(genes)), genes)
  zero_1 <- select_fun(beta_zero, candidate_rank, params)
  zero_2 <- select_fun(beta_zero, candidate_rank, params)
  assert_equal(as.character(zero_1), as.character(zero_2),
               message = "all-zero fallback is not deterministic")
  assert_equal(as.character(zero_1), head(candidate_rank, 8L),
               message = "all-zero fallback did not use the locked univariate order")

  beta_sparse <- beta_zero
  beta_sparse[c("G04", "G03", "G02", "G01")] <- c(1, -1, 0.5, -0.5)
  sparse <- select_fun(beta_sparse, candidate_rank, params)
  expected_prefix <- c("G03", "G04", "G01", "G02")
  assert_equal(as.character(head(sparse, 4L)), expected_prefix,
               message = "absolute-coefficient ties were not broken by gene name")
  assert_true(length(sparse) == 8L && !anyDuplicated(sparse),
              "1-4 nonzero coefficients did not deterministically fill to 8 genes")

  beta_many <- stats::setNames(rep(0.2, length(genes)), rev(genes))
  many <- select_fun(beta_many, candidate_rank, params)
  assert_equal(as.character(many), head(sort(names(beta_many)), 30L),
               message = ">30-gene cap/tie rule is not deterministic")

  dat <- make_synthetic_survival(n = 48L, p = 16L, seed = 8001L)
  fallback_params <- synthetic_params()
  fallback_params$univ_p_cutoff <- 0
  fallback_params$candidate_fallback_n <- 10L
  prep_1 <- get("prepare_training", envir = pipeline_env)(
    dat$expr, dat$time, dat$status, fallback_params
  )
  prep_2 <- get("prepare_training", envir = pipeline_env)(
    dat$expr, dat$time, dat$status, fallback_params
  )
  assert_equal(prep_1$candidates, prep_2$candidates,
               message = "candidate fallback changed across identical calls")
  assert_equal(prep_1$candidates, head(prep_1$univ_rank, 10L),
               message = "candidate fallback did not follow locked univariate rank")
})

run_test("pipeline version is mandatory and propagated into nested metadata", {
  block_if_missing(c("validate_pipeline_params", "run_nested_cv", "validate_nested_result"))
  validate_params <- get("validate_pipeline_params", envir = pipeline_env)
  base <- synthetic_params()
  validate_params(base)
  invalid <- list(
    within(base, rm(pipeline_version)),
    modifyList(base, list(pipeline_version = NA_character_)),
    modifyList(base, list(pipeline_version = "")),
    modifyList(base, list(pipeline_version = c("v1", "v2")))
  )
  for (candidate in invalid) {
    expect_error(validate_params(candidate), "pipeline_version")
  }
  dat <- make_synthetic_survival(n = 48L, p = 16L, seed = 7033L)
  nested <- get("run_nested_cv", envir = pipeline_env)(
    dat$expr, dat$time, dat$status, dat$sample_ids, base,
    analysis_key = "pipeline-version-contract-key"
  )
  assert_true(identical(nested$pipeline_version, base$pipeline_version),
              "top-level pipeline version was not propagated")
  assert_true(identical(nested$summary$pipeline_version, base$pipeline_version),
              "summary pipeline version was not propagated")
  tampered <- nested
  tampered$pipeline_version <- "tampered-version"
  expect_error(
    get("validate_nested_result", envir = pipeline_env)(
      tampered, dat$expr, dat$time, dat$status, dat$sample_ids,
      base, nested$analysis_key
    ),
    "pipeline_version"
  )
})

run_test("Cox refit preserves non-syntactic gene SYMBOLs", {
  block_if_missing("fit_signature_pipeline")
  dat <- make_synthetic_survival(n = 60L, p = 18L, seed = 8123L)
  rownames(dat$expr) <- sprintf("GENE-%03d", seq_len(nrow(dat$expr)))
  params <- synthetic_params()
  fit_fun <- get("fit_signature_pipeline", envir = pipeline_env)

  model <- fit_fun(
    dat$expr, dat$time, dat$status, seed = 444L, params = params
  )
  assert_true(any(make.names(model$genes) != model$genes),
              "fixture did not retain a non-syntactic selected SYMBOL")
  assert_equal(names(model$coefs), model$genes,
               message = "refit coefficient names did not preserve raw SYMBOLs")
  assert_equal(names(stats::coef(model$fit)), model$genes,
               message = "coxph object did not preserve raw SYMBOLs")
  prediction <- get("predict_signature_pipeline", envir = pipeline_env)(
    model, dat$expr
  )
  assert_true(all(is.finite(prediction$lp)) &&
                all(is.finite(prediction$lp_train_sd)),
              "prediction failed after restoring non-syntactic SYMBOLs")
})

run_test("fit is reproducible and inner validation is stratified", {
  block_if_missing("fit_signature_pipeline")
  dat <- make_synthetic_survival()
  params <- synthetic_params()
  fit_fun <- get("fit_signature_pipeline", envir = pipeline_env)
  model_1 <- fit_fun(dat$expr, dat$time, dat$status, seed = 333L,
                     params = params)
  model_2 <- fit_fun(dat$expr, dat$time, dat$status, seed = 333L,
                     params = params)
  assert_equal(model_fingerprint(model_1), model_fingerprint(model_2),
               message = "model fingerprint changed for identical seed/input")
  assert_equal(model_1$selected, model_2$selected,
               message = "selected genes changed for identical seed/input")
  assert_equal(model_1$coefs, model_2$coefs, tolerance = 1e-9,
               message = "coefficients changed for identical seed/input")
  assert_true(grepl("^[0-9a-f]{64}$", model_fingerprint(model_1)),
              "model_fingerprint is not a SHA-256 hex digest")

  fingerprint_payload <- list(
    genes = model_1$genes,
    coefficients = unname(model_1$coefs),
    center = unname(model_1$center[model_1$genes]),
    scale = unname(model_1$scale[model_1$genes]),
    lp_center = model_1$train_lp_mean,
    lp_scale = model_1$train_lp_sd,
    lambda_ratio = model_1$lambda_ratio,
    lambda = model_1$lambda
  )
  expected_fingerprint <- digest::digest(
    fingerprint_payload, algo = "sha256", serialize = TRUE
  )
  assert_equal(model_fingerprint(model_1), expected_fingerprint,
               message = "model fingerprint does not match predictive model fields")

  tuning <- as.data.frame(model_1$tuning_curve, stringsAsFactors = FALSE)
  required_tuning <- c(
    "lambda_ratio", "mean_c_index", "se_c_index", "best_mean",
    "one_se_eligible", "selected", "one_se_threshold"
  )
  assert_true(all(required_tuning %in% names(tuning)),
              "tuning curve lacks one-SE audit columns")
  assert_true(sum(tuning$selected) == 1L, "one-SE rule selected != 1 ratio")
  best_mean <- max(tuning$mean_c_index)
  best_rows <- tuning[tuning$mean_c_index == best_mean, , drop = FALSE]
  best_row <- best_rows[which.max(best_rows$lambda_ratio), , drop = FALSE]
  threshold <- best_row$mean_c_index - best_row$se_c_index
  expected_ratio <- max(tuning$lambda_ratio[tuning$mean_c_index >= threshold])
  assert_equal(unique(tuning$one_se_threshold), threshold, tolerance = 1e-12,
               message = "stored one-SE threshold is incorrect")
  assert_equal(tuning$lambda_ratio[tuning$selected], expected_ratio,
               tolerance = 1e-12,
               message = "one-SE rule did not choose the largest eligible penalty")

  inner <- assignment_table(model_1$inner_fold_assignments, dat$sample_ids,
                            "inner_fold")
  assert_true(nrow(inner) == length(dat$sample_ids),
              "inner assignments do not cover all training samples")
  assert_true(!anyDuplicated(inner$sample),
              "a sample occurs more than once in inner assignments")
  status_lookup <- stats::setNames(dat$status, dat$sample_ids)
  for (j in sort(unique(inner$fold_value))) {
    y <- unname(status_lookup[inner$sample[inner$fold_value == j]])
    assert_set_equal(unique(y), 0:1,
                     paste("inner validation fold", j, "is not class-balanced"))
  }
  artifacts$model <- model_1
  artifacts$synthetic <- dat
})

run_test("leakage sentinel and training-only standardization hold", {
  block_if_missing(c("make_stratified_folds", "fit_signature_pipeline",
                     "predict_signature_pipeline"))
  dat <- if (exists("synthetic", envir = artifacts, inherits = FALSE)) {
    artifacts$synthetic
  } else {
    make_synthetic_survival()
  }
  params <- synthetic_params()
  folds <- get("make_stratified_folds", envir = pipeline_env)(
    dat$status, k = 3L, seed = 444L
  )
  train_idx <- which(folds != 1L)
  test_idx <- which(folds == 1L)
  fit_fun <- get("fit_signature_pipeline", envir = pipeline_env)
  pred_fun <- get("predict_signature_pipeline", envir = pipeline_env)

  model_1 <- fit_fun(dat$expr[, train_idx, drop = FALSE], dat$time[train_idx],
                     dat$status[train_idx], seed = 445L, params = params)
  status_sentinel <- dat$status
  status_sentinel[test_idx] <- 1L - status_sentinel[test_idx]
  model_2 <- fit_fun(dat$expr[, train_idx, drop = FALSE], dat$time[train_idx],
                     status_sentinel[train_idx], seed = 445L, params = params)
  assert_equal(model_fingerprint(model_1), model_fingerprint(model_2),
               message = "changing outer-test outcomes changed the trained model")

  selected <- names(model_1$coefs)
  assert_true(length(selected) > 0L, "model returned no named coefficients")
  direct_center <- rowMeans(dat$expr[selected, train_idx, drop = FALSE])
  direct_scale <- apply(dat$expr[selected, train_idx, drop = FALSE], 1L, stats::sd)
  assert_equal(unname(model_1$center[selected]), unname(direct_center), tolerance = 1e-10,
               message = "stored centers are not outer-training means")
  assert_equal(unname(model_1$scale[selected]), unname(direct_scale), tolerance = 1e-10,
               message = "stored scales are not outer-training SDs")
  test_center <- rowMeans(dat$expr[selected, test_idx, drop = FALSE])
  assert_true(any(abs(unname(model_1$center[selected]) - unname(test_center)) > 1e-6),
              "sentinel is uninformative because train and test centers coincide")

  pred_1 <- pred_fun(model_1, dat$expr[, test_idx, drop = FALSE])
  changed_expr <- dat$expr[, test_idx, drop = FALSE]
  changed_expr[selected[[1L]], ] <- changed_expr[selected[[1L]], ] +
    seq_len(ncol(changed_expr)) / 3
  pred_2 <- pred_fun(model_1, changed_expr)
  assert_true(any(abs(pred_1$lp_train_sd - pred_2$lp_train_sd) > 1e-8),
              "changing outer-test expression did not change predictions")
  assert_equal(model_fingerprint(model_1), model_fingerprint(model_2),
               message = "test-expression sentinel changed the model fingerprint")
})

run_test("one-repeat synthetic nested CV has exact OOF coverage and disjoint folds", {
  block_if_missing("run_nested_cv")
  dat <- make_synthetic_survival()
  params <- synthetic_params()
  nested <- get("run_nested_cv", envir = pipeline_env)(
    dat$expr, dat$time, dat$status, dat$sample_ids, params,
    analysis_key = "synthetic-contract-key"
  )
  expected_components <- c(
    "outer_fold_assignments", "inner_fold_assignments", "oof", "fold_performance",
    "repeat_performance", "selection", "selection_frequency", "jaccard",
    "summary", "failure_log"
  )
  assert_true(all(expected_components %in% names(nested)),
              paste("run_nested_cv() missing components:",
                    paste(setdiff(expected_components, names(nested)), collapse = ", ")))

  oof <- as.data.frame(nested$oof, stringsAsFactors = FALSE)
  required_oof <- c(
    "sample", "repeat_id", "outer_fold", "time", "status", "lp_oof",
    "n_train", "events_train", "n_test", "events_test", "n_genes",
    "analysis_key"
  )
  assert_true(all(required_oof %in% names(oof)),
              paste("OOF table missing columns:",
                    paste(setdiff(required_oof, names(oof)), collapse = ", ")))
  assert_true(nrow(oof) == length(dat$sample_ids) * params$outer_repeats,
              "OOF row count is not n_samples * outer_repeats")
  assert_true(!anyDuplicated(oof[c("sample", "repeat_id")]),
              "a patient has more than one OOF prediction in a repeat")
  assert_set_equal(oof$sample, dat$sample_ids, "OOF table omits samples")
  assert_true(all(is.finite(oof$lp_oof)), "OOF predictions contain non-finite values")
  assert_true(all(oof$analysis_key == "synthetic-contract-key"),
              "analysis_key was not propagated to every OOF row")
  assert_true(nrow(as.data.frame(nested$fold_performance)) ==
                params$outer_folds * params$outer_repeats,
              "completed fold count is not repeats * folds")

  assignments <- assignment_table(nested$outer_fold_assignments,
                                  dat$sample_ids, "outer_fold")
  status_lookup <- stats::setNames(dat$status, dat$sample_ids)
  for (r in sort(unique(assignments$repeat_id))) {
    ar <- assignments[assignments$repeat_id == r, , drop = FALSE]
    assert_true(nrow(ar) == length(dat$sample_ids) && !anyDuplicated(ar$sample),
                paste("repeat", r, "does not assign every sample exactly once"))
    for (j in sort(unique(ar$fold_value))) {
      test_samples <- ar$sample[ar$fold_value == j]
      train_samples <- setdiff(dat$sample_ids, test_samples)
      assert_true(!length(intersect(train_samples, test_samples)),
                  paste("repeat", r, "fold", j, "train/test overlap"))
      assert_set_equal(union(train_samples, test_samples), dat$sample_ids,
                       paste("repeat", r, "fold", j, "does not partition samples"))
      y <- unname(status_lookup[test_samples])
      assert_set_equal(unique(y), 0:1,
                       paste("repeat", r, "outer fold", j,
                             "does not contain events and censoring"))
      oo <- oof[oof$repeat_id == r & oof$outer_fold == j, , drop = FALSE]
      assert_true(all(oo$n_test == length(test_samples)) &&
                    all(oo$n_train == length(train_samples)),
                  "OOF fold sizes disagree with assignments")
      assert_true(all(oo$events_test == sum(status_lookup[test_samples])) &&
                    all(oo$events_train == sum(status_lookup[train_samples])),
                  "OOF event counts disagree with assignments")

      inner <- as.data.frame(nested$inner_fold_assignments,
                             stringsAsFactors = FALSE)
      required_inner <- c("sample", "repeat_id", "outer_fold", "inner_fold")
      assert_true(all(required_inner %in% names(inner)),
                  "inner-fold audit table is missing identifiers")
      ii <- inner[inner$repeat_id == r & inner$outer_fold == j, , drop = FALSE]
      assert_set_equal(ii$sample, train_samples,
                       "inner-fold assignments do not cover the outer training set")
      assert_true(!anyDuplicated(ii$sample),
                  "an outer-training sample has multiple inner-fold assignments")
      for (inner_j in sort(unique(ii$inner_fold))) {
        inner_y <- unname(status_lookup[ii$sample[ii$inner_fold == inner_j]])
        assert_set_equal(unique(inner_y), 0:1,
                         paste("repeat", r, "outer fold", j, "inner fold",
                               inner_j, "lacks events or censoring"))
      }
    }
  }
  assert_true(!nrow(as.data.frame(nested$failure_log)),
              "successful synthetic run contains failure-log rows")
  artifacts$nested_1 <- nested
})

run_test("repeated OOF aggregation rejects inconsistent patient outcomes", {
  block_if_missing("summarize_repeated_oof")
  summarize_oof <- get("summarize_repeated_oof", envir = pipeline_env)
  oof <- data.frame(
    sample = rep(c("P1", "P2"), each = 2L),
    repeat_id = rep(1:2, times = 2L),
    time = rep(c(100, 200), each = 2L),
    status = rep(c(1L, 0L), each = 2L),
    lp_oof = c(-0.2, -0.1, 0.3, 0.4),
    stringsAsFactors = FALSE
  )
  ensemble <- summarize_oof(oof, expected_repeats = 2L)
  assert_true(nrow(ensemble) == 2L && !anyDuplicated(ensemble$sample),
              "valid repeated OOF rows were not reduced one-per-patient")

  inconsistent_time <- oof
  inconsistent_time$time[2L] <- inconsistent_time$time[2L] + 1
  expect_error(
    summarize_oof(inconsistent_time, expected_repeats = 2L),
    "inconsistent within patient"
  )

  duplicated_repeat <- oof
  duplicated_repeat$repeat_id[2L] <- duplicated_repeat$repeat_id[1L]
  expect_error(
    summarize_oof(duplicated_repeat, expected_repeats = 2L),
    "inconsistent within patient"
  )
})

run_test("Jaccard stability is calculated within each repeat", {
  block_if_missing("run_nested_cv")
  dat <- make_synthetic_survival(n = 66L, p = 30L, seed = 7027L)
  params <- synthetic_params()
  params$outer_repeats <- 2L
  params$lambda_ratio_grid <- c(0.30, 0.08)
  nested <- get("run_nested_cv", envir = pipeline_env)(
    dat$expr, dat$time, dat$status, dat$sample_ids, params,
    analysis_key = "synthetic-jaccard-key"
  )
  selection <- as.data.frame(nested$selection, stringsAsFactors = FALSE)
  jaccard <- as.data.frame(nested$jaccard, stringsAsFactors = FALSE)
  required <- c("repeat_id", "fold_a", "fold_b", "intersection_n", "union_n",
                "jaccard")
  assert_true(all(required %in% names(jaccard)),
              "Jaccard table lacks repeat/fold audit columns")
  assert_true(nrow(jaccard) == params$outer_repeats * base::choose(params$outer_folds, 2L),
              "Jaccard table has the wrong number of within-repeat fold pairs")

  repeat_gene_sets <- lapply(seq_len(params$outer_repeats), function(r) {
    lapply(seq_len(params$outer_folds), function(j) {
      sort(selection$gene[selection$repeat_id == r & selection$outer_fold == j])
    })
  })
  assert_true(!identical(repeat_gene_sets[[1L]], repeat_gene_sets[[2L]]),
              "Jaccard regression sentinel is uninformative: repeat selections coincide")

  for (i in seq_len(nrow(jaccard))) {
    r <- jaccard$repeat_id[[i]]
    a <- selection$gene[
      selection$repeat_id == r & selection$outer_fold == jaccard$fold_a[[i]]
    ]
    b <- selection$gene[
      selection$repeat_id == r & selection$outer_fold == jaccard$fold_b[[i]]
    ]
    expected_intersection <- length(intersect(a, b))
    expected_union <- length(union(a, b))
    expected_jaccard <- expected_intersection / expected_union
    assert_true(jaccard$intersection_n[[i]] == expected_intersection &&
                  jaccard$union_n[[i]] == expected_union,
                paste("Jaccard counts use the wrong repeat for row", i))
    assert_equal(jaccard$jaccard[[i]], expected_jaccard, tolerance = 1e-12,
                 message = paste("Jaccard value uses the wrong repeat for row", i))
  }
})

run_test("nested CV is reproducible for identical input and seed", {
  block_if_missing("run_nested_cv")
  dat <- make_synthetic_survival()
  params <- synthetic_params()
  nested_1 <- if (exists("nested_1", envir = artifacts, inherits = FALSE)) {
    artifacts$nested_1
  } else {
    get("run_nested_cv", envir = pipeline_env)(
      dat$expr, dat$time, dat$status, dat$sample_ids, params,
      analysis_key = "synthetic-contract-key"
    )
  }
  nested_2 <- get("run_nested_cv", envir = pipeline_env)(
    dat$expr, dat$time, dat$status, dat$sample_ids, params,
    analysis_key = "synthetic-contract-key"
  )
  a1 <- as.data.frame(nested_1$outer_fold_assignments)
  a2 <- as.data.frame(nested_2$outer_fold_assignments)
  assert_equal(a1, a2, message = "outer assignments are not reproducible")
  o1 <- as.data.frame(nested_1$oof)
  o2 <- as.data.frame(nested_2$oof)
  ord1 <- order(o1$repeat_id, o1$sample)
  ord2 <- order(o2$repeat_id, o2$sample)
  assert_equal(o1[ord1, c("sample", "repeat_id", "outer_fold")],
               o2[ord2, c("sample", "repeat_id", "outer_fold")],
               message = "OOF identities/folds are not reproducible")
  assert_equal(o1$lp_oof[ord1], o2$lp_oof[ord2], tolerance = 1e-9,
               message = "OOF predictions are not reproducible")
  s1 <- as.data.frame(nested_1$selection)
  s2 <- as.data.frame(nested_2$selection)
  assert_equal(s1, s2, message = "fold-level gene selections are not reproducible")
})

run_test("cache key invalidates on parameter, input, and code changes", {
  block_if_missing("make_analysis_key")
  make_key <- get("make_analysis_key", envir = pipeline_env)
  td <- tempfile("crc9-cache-key-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE, force = TRUE), add = TRUE)
  input_path <- file.path(td, "input.txt")
  code_path <- file.path(td, "code.R")
  writeLines(c("sample,value", "S1,1"), input_path, useBytes = TRUE)
  writeLines("identity <- function(x) x", code_path, useBytes = TRUE)
  params <- synthetic_params()
  identity <- make_key(input_path, code_path, params, return_manifest = TRUE)
  assert_true(is.list(identity) && all(c("key", "manifest") %in% names(identity)),
              "return_manifest=TRUE lacks key/manifest")
  assert_true(all(c("inputs", "code", "params", "r_version", "packages") %in%
                    names(identity$manifest)),
              "cache manifest is incomplete")
  assert_true(grepl("^[0-9a-f]{64}$", identity$key),
              "analysis key is not a SHA-256 hex digest")
  assert_true(all(grepl("^[0-9a-f]{64}$", c(identity$manifest$inputs,
                                             identity$manifest$code))),
              "manifest contains a non-SHA-256 file hash")
  base_1 <- identity$key
  base_2 <- key_value(make_key(input_path, code_path, params))
  assert_equal(base_1, base_2, message = "identical cache inputs changed the key")

  params_changed <- params
  params_changed$coef_cutoff <- params_changed$coef_cutoff + 0.001
  param_key <- key_value(make_key(input_path, code_path, params_changed))
  assert_true(!identical(base_1, param_key), "parameter change did not invalidate key")

  writeLines(c("sample,value", "S1,2"), input_path, useBytes = TRUE)
  input_key <- key_value(make_key(input_path, code_path, params))
  assert_true(!identical(base_1, input_key), "input-content change did not invalidate key")
  writeLines(c("sample,value", "S1,1"), input_path, useBytes = TRUE)

  writeLines("identity <- function(x) x + 0", code_path, useBytes = TRUE)
  code_key <- key_value(make_key(input_path, code_path, params))
  assert_true(!identical(base_1, code_key), "code-content change did not invalidate key")
})

run_test("run key covers all inputs and reporting parameters", {
  block_if_missing("make_run_key")
  make_key <- get("make_run_key", envir = pipeline_env)
  td <- tempfile("crc9-run-key-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE, force = TRUE), add = TRUE)
  training_path <- file.path(td, "training.csv")
  external_path <- file.path(td, "external.csv")
  code_path <- file.path(td, "pipeline.R")
  writeLines("training", training_path, useBytes = TRUE)
  writeLines("external", external_path, useBytes = TRUE)
  writeLines("pipeline", code_path, useBytes = TRUE)
  params <- list(nested_repeats = 1L, bootstrap_reps = 100L,
                 primary_mapping = "unique_mean")
  base <- make_key(
    c(training_path, external_path), code_path, params,
    return_manifest = TRUE
  )
  assert_true(grepl("^[0-9a-f]{64}$", base$key),
              "run key is not a SHA-256 hex digest")

  changed_bootstrap <- params
  changed_bootstrap$bootstrap_reps <- 200L
  assert_true(!identical(
    base$key,
    make_key(c(training_path, external_path), code_path, changed_bootstrap)
  ), "bootstrap parameter change did not invalidate run key")

  writeLines("external changed", external_path, useBytes = TRUE)
  assert_true(!identical(
    base$key,
    make_key(c(training_path, external_path), code_path, params)
  ), "external input change did not invalidate run key")
})

run_test("cached nested result is revalidated against source identities", {
  block_if_missing(c("run_nested_cv", "validate_nested_result"))
  dat <- make_synthetic_survival()
  params <- synthetic_params()
  analysis_key <- "synthetic-cache-validation-key"
  nested <- get("run_nested_cv", envir = pipeline_env)(
    dat$expr, dat$time, dat$status, dat$sample_ids, params, analysis_key
  )
  validate <- get("validate_nested_result", envir = pipeline_env)
  cache_path <- tempfile("crc9-nested-cache-", fileext = ".rds")
  on.exit(unlink(cache_path, force = TRUE), add = TRUE)
  saveRDS(nested, cache_path, version = 3)
  cached <- readRDS(cache_path)
  validate(cached, dat$expr, dat$time, dat$status, dat$sample_ids,
           params, analysis_key)

  wrong_sample <- cached
  wrong_sample$oof$sample[[1L]] <- "NOT_A_SOURCE_SAMPLE"
  expect_error(validate(
    wrong_sample, dat$expr, dat$time, dat$status, dat$sample_ids,
    params, analysis_key
  ))

  wrong_outcome <- cached
  wrong_outcome$oof$time[[1L]] <- wrong_outcome$oof$time[[1L]] + 1000
  wrong_outcome$oof$status[[2L]] <- 1L - wrong_outcome$oof$status[[2L]]
  expect_error(validate(
    wrong_outcome, dat$expr, dat$time, dat$status, dat$sample_ids,
    params, analysis_key
  ))

  wrong_fold <- cached
  wrong_fold$oof$outer_fold[[1L]] <-
    (wrong_fold$oof$outer_fold[[1L]] %% params$outer_folds) + 1L
  expect_error(validate(
    wrong_fold, dat$expr, dat$time, dat$status, dat$sample_ids,
    params, analysis_key
  ))
})

run_test("cached nested summary tampering is rejected field by field", {
  block_if_missing(c("run_nested_cv", "validate_nested_result"))
  dat <- if (exists("synthetic", envir = artifacts, inherits = FALSE)) {
    artifacts$synthetic
  } else {
    make_synthetic_survival()
  }
  params <- synthetic_params()
  cached <- if (exists("nested_1", envir = artifacts, inherits = FALSE)) {
    artifacts$nested_1
  } else {
    get("run_nested_cv", envir = pipeline_env)(
      dat$expr, dat$time, dat$status, dat$sample_ids, params,
      analysis_key = "synthetic-contract-key"
    )
  }
  validate <- get("validate_nested_result", envir = pipeline_env)
  summary_fields <- c(
    "pipeline_version", "analysis_key", "repeats", "outer_folds", "inner_folds",
    "completed_folds", "median_repeat_oof_c_index",
    "repeat_oof_c_index_q1", "repeat_oof_c_index_q3",
    "min_repeat_oof_c_index", "max_repeat_oof_c_index",
    "median_within_repeat_jaccard", "jaccard_q1", "jaccard_q3"
  )
  for (field in summary_fields) {
    tampered <- cached
    original <- tampered$summary[[field]]
    tampered$summary[[field]] <- if (is.character(original)) {
      paste0(original, "-tampered")
    } else if (is.integer(original)) {
      original + 1L
    } else {
      original + 0.01
    }
    expect_error(
      validate(
        tampered, dat$expr, dat$time, dat$status, dat$sample_ids,
        params, cached$analysis_key
      ),
      if (identical(field, "pipeline_version")) {
        "pipeline_version"
      } else {
        paste0("summary field mismatch: ", field)
      }
    )
  }
})

run_test("an unfit outer fold aborts the entire nested run", {
  block_if_missing("run_nested_cv")
  params <- synthetic_params()
  n <- 30L
  expr <- matrix(1, nrow = 12L, ncol = n,
                 dimnames = list(sprintf("CONST%02d", 1:12),
                                 sprintf("BAD%02d", seq_len(n))))
  time <- seq_len(n) + 1
  status <- rep(c(0L, 1L), length.out = n)
  err <- expect_error(get("run_nested_cv", envir = pipeline_env)(
    expr, time, status, colnames(expr), params,
    analysis_key = "must-fail-contract-key"
  ))
  assert_true(nzchar(conditionMessage(err)), "failure was raised without a diagnostic")
})

run_test("final fit and outer fits share one training implementation", {
  parsed_utils <- parse(file = utils_file, keep.source = TRUE)
  parsed_core <- parse(file = core_file, keep.source = TRUE)
  fit_defs <- c(function_definitions(parsed_utils, "fit_signature_pipeline"),
                function_definitions(parsed_core, "fit_signature_pipeline"))
  assert_true(length(fit_defs) == 1L,
              paste("expected exactly one fit_signature_pipeline definition; found",
                    length(fit_defs)))
  run_defs <- c(function_definitions(parsed_utils, "run_nested_cv"),
                function_definitions(parsed_core, "run_nested_cv"))
  assert_true(length(run_defs) == 1L,
              paste("expected exactly one run_nested_cv definition; found",
                    length(run_defs)))
  run_calls <- called_symbols(run_defs[[1L]]$definition)
  assert_true("fit_signature_pipeline" %in% run_calls,
              "run_nested_cv() does not call fit_signature_pipeline()")
  assert_true(!"next" %in% run_calls,
              "run_nested_cv() contains a silent next path")
  core_calls <- called_symbols(as.call(c(quote(`{`), as.list(parsed_core))))
  assert_true("fit_signature_pipeline" %in% core_calls,
              "core script does not use fit_signature_pipeline() for the final model")
  assert_true("run_nested_cv" %in% core_calls,
              "core script does not use run_nested_cv() for outer validation")
})

summary_df <- do.call(rbind, results)
cat("\nTest summary\n")
print(summary_df, row.names = FALSE)
cat(sprintf("\nPASS=%d FAIL=%d BLOCKED=%d\n",
            sum(summary_df$status == "PASS"),
            sum(summary_df$status == "FAIL"),
            sum(summary_df$status == "BLOCKED")))

if (any(summary_df$status != "PASS")) quit(save = "no", status = 1L)
quit(save = "no", status = 0L)
