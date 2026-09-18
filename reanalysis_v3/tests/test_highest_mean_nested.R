#!/usr/bin/env Rscript

# Synthetic regression tests for strict training-partition-only highest-mean
# GPL570 probe selection. No public cohort or previously learned map is used.

options(stringsAsFactors = FALSE, warn = 1)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript")
this_file <- gsub("\\\\", "/", sub("^--file=", "", script_arg[[1L]]))
test_dir <- dirname(this_file)
v3_root <- dirname(test_dir)
work_root <- dirname(v3_root)
utils_file <- file.path(v3_root, "scripts", "utils.R")
probe_file <- file.path(v3_root, "scripts", "highest_mean_sensitivity.R")
if (!file.exists(utils_file) || !file.exists(probe_file)) {
  stop("Cannot locate strict probe-sensitivity scripts relative to test file")
}

Sys.setenv(CRC_WORK_ROOT = work_root)
if (!nzchar(Sys.getenv("CRC_SOURCE_ROOT"))) {
  Sys.setenv(CRC_SOURCE_ROOT = work_root)
}
source(utils_file)
source(probe_file)
validator_env <- new.env(parent = globalenv())
sys.source(
  file.path(test_dir, "validate_completed_run.R"),
  envir = validator_env
)

results <- list()
artifacts <- new.env(parent = emptyenv())

assert_true <- function(value, message = "assertion failed") {
  if (length(value) != 1L || is.na(value) || !isTRUE(value)) stop(message)
  invisible(TRUE)
}

assert_equal <- function(actual, expected, tolerance = 1e-10,
                         message = "objects differ") {
  if (!isTRUE(all.equal(actual, expected, tolerance = tolerance,
                        check.attributes = TRUE))) {
    stop(message, ": ", paste(
      all.equal(actual, expected, tolerance = tolerance,
                check.attributes = TRUE), collapse = "; "
    ))
  }
  invisible(TRUE)
}

expect_error <- function(expr, pattern = NULL) {
  captured <- NULL
  tryCatch(force(expr), error = function(e) captured <<- e)
  if (is.null(captured)) stop("expected an error, but expression returned normally")
  if (!is.null(pattern) && !grepl(
    pattern, conditionMessage(captured), ignore.case = TRUE
  )) {
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

probe_params <- function() {
  list(
    pipeline_version = "probe-highest-mean-synthetic-v1",
    outer_folds = 3L,
    outer_repeats = 1L,
    inner_folds = 2L,
    base_seed = 20260917L,
    sd_cutoff = 0.01,
    univ_p_cutoff = 0.30,
    min_candidates = 3L,
    candidate_fallback_n = 10L,
    coef_cutoff = 0.03,
    min_selected = 2L,
    selected_fallback_n = 4L,
    max_selected = 5L,
    lambda_ratio_grid = c(0.35, 0.08),
    lambda_rule = "one_se_larger_penalty",
    ties = "efron"
  )
}

make_probe_survival <- function(n = 66L, genes = 16L, seed = 91701L) {
  set.seed(seed)
  gene_expr <- matrix(stats::rnorm(genes * n), nrow = genes, ncol = n)
  gene_names <- sprintf("GENE%02d", seq_len(genes))
  sample_ids <- sprintf("P%03d", seq_len(n))
  rownames(gene_expr) <- gene_names
  colnames(gene_expr) <- sample_ids

  probe_ids <- as.vector(rbind(
    sprintf("G%02d_P1", seq_len(genes)),
    sprintf("G%02d_P2", seq_len(genes))
  ))
  probe_expr <- matrix(
    NA_real_, nrow = length(probe_ids) + 2L, ncol = n,
    dimnames = list(c(probe_ids, "MULTI_1", "NO_SYMBOL_1"), sample_ids)
  )
  for (g in seq_len(genes)) {
    probe_expr[sprintf("G%02d_P1", g), ] <-
      gene_expr[g, ] + stats::rnorm(n, sd = 0.06) + 0.20
    probe_expr[sprintf("G%02d_P2", g), ] <-
      gene_expr[g, ] + stats::rnorm(n, sd = 0.06)
  }
  probe_expr["MULTI_1", ] <- gene_expr[1L, ] + gene_expr[2L, ]
  probe_expr["NO_SYMBOL_1", ] <- stats::rnorm(n)

  annotation <- data.frame(
    PROBEID = c(
      probe_ids, "MULTI_1", "MULTI_1", "NO_SYMBOL_1"
    ),
    SYMBOL = c(
      rep(gene_names, each = 2L), gene_names[1L], gene_names[2L], ""
    ),
    stringsAsFactors = FALSE
  )
  lp <- 0.95 * gene_expr[1L, ] - 0.75 * gene_expr[2L, ] +
    0.55 * gene_expr[3L, ]
  event_time <- stats::rexp(n, rate = exp(lp) / 14)
  censor_time <- stats::rexp(n, rate = 1 / 21)
  time <- pmax(pmin(event_time, censor_time) * 30.4375, 0.01)
  status <- as.integer(event_time <= censor_time)
  if (sum(status) < 15L || sum(status == 0L) < 15L) {
    stop("Synthetic probe data have inadequate outcome balance")
  }
  list(
    probe_expr = probe_expr,
    annotation = annotation,
    time = time,
    status = status,
    sample_ids = sample_ids
  )
}

ordered_rows <- function(x, columns) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  ord <- do.call(order, x[columns])
  out <- x[ord, , drop = FALSE]
  rownames(out) <- NULL
  out
}

required_functions <- c(
  "learn_unique_highest_probe_map", "apply_frozen_probe_map",
  "hash_probe_map", "fit_probe_signature_pipeline",
  "predict_probe_signature_pipeline", "run_probe_nested_cv",
  "validate_probe_nested_result"
)

run_test("all strict probe-sensitivity interfaces are present", {
  missing <- required_functions[!vapply(
    required_functions, exists, logical(1), mode = "function", inherits = TRUE
  )]
  assert_true(!length(missing), paste("missing:", paste(missing, collapse = ", ")))
})

run_test("exact mean ties use lexicographic PROBEID", {
  x <- rbind(
    "100_at" = c(2, 3, 4, 5),
    "200_at" = c(5, 4, 3, 2),
    "300_at" = c(1, 2, 1, 2),
    "MULTI" = c(9, 9, 9, 9)
  )
  colnames(x) <- paste0("S", seq_len(ncol(x)))
  ann <- data.frame(
    PROBEID = c("100_at", "200_at", "300_at", "MULTI", "MULTI"),
    SYMBOL = c("A", "A", "B", "A", "B"),
    stringsAsFactors = FALSE
  )
  map <- learn_unique_highest_probe_map(x, ann)
  assert_true(
    map$mapping$PROBEID[map$mapping$SYMBOL == "A"] == "100_at",
    "lexicographic tie-break did not select 100_at"
  )
  assert_true(
    map$mapping$n_probes_gene[map$mapping$SYMBOL == "A"] == 2L,
    "eligible probe count for A is wrong"
  )
  assert_true(!"MULTI" %in% map$mapping$PROBEID,
              "multi-SYMBOL probe was not excluded")
  assert_true(identical(map$map_hash, hash_probe_map(map)),
              "stored map hash is not reproducible")
})

run_test("probe-map hash survives RDS and CSV round trips", {
  set.seed(91709L)
  probe_ids <- sprintf("HASH_P%02d_at", seq_len(12L))
  sample_ids <- sprintf("HASH_S%02d", seq_len(30L))
  probe_expr <- matrix(
    stats::rnorm(length(probe_ids) * length(sample_ids)),
    nrow = length(probe_ids),
    dimnames = list(probe_ids, sample_ids)
  )
  annotation <- data.frame(
    PROBEID = probe_ids,
    SYMBOL = rep(c("A-1", "A1", "A_1", "AA", "A.B", "AB"), each = 2L),
    stringsAsFactors = FALSE
  )
  map <- learn_unique_highest_probe_map(
    probe_expr, annotation, sample_ids
  )
  expected <- map$map_hash

  rds_path <- tempfile(fileext = ".rds")
  csv_path <- tempfile(fileext = ".csv")
  on.exit(unlink(c(rds_path, csv_path)), add = TRUE)
  saveRDS(map$mapping, rds_path, version = 3)
  utils::write.csv(map$mapping, csv_path, row.names = FALSE, quote = TRUE)

  from_rds <- readRDS(rds_path)
  from_csv <- utils::read.csv(
    csv_path, stringsAsFactors = FALSE, check.names = FALSE
  )
  assert_true(identical(hash_probe_map(from_rds), expected),
              "RDS round trip changed the probe-map hash")
  assert_true(identical(hash_probe_map(from_csv), expected),
              "CSV round trip changed the probe-map hash")
  assert_true(
    identical(validator_env$qa_hash_probe_map(from_csv), expected),
    "production and completed-run validator map hashes disagree"
  )

  original_collate <- Sys.getlocale("LC_COLLATE")
  on.exit(suppressWarnings(Sys.setlocale("LC_COLLATE", original_collate)),
          add = TRUE)
  locale_names <- c("C", "Chinese_China.utf8")
  default_orders <- locale_hashes <- lapply(locale_names, function(locale_name) {
    active <- suppressWarnings(Sys.setlocale("LC_COLLATE", locale_name))
    assert_true(nzchar(active), paste("unavailable test locale:", locale_name))
    default_order <- paste(
      from_csv$SYMBOL[order(from_csv$SYMBOL, from_csv$PROBEID)],
      collapse = "\r"
    )
    hashes <- c(
      production = hash_probe_map(from_csv),
      validator = validator_env$qa_hash_probe_map(from_csv)
    )
    list(default_order = default_order, hashes = hashes)
  })
  default_order_values <- vapply(
    default_orders, `[[`, character(1), "default_order"
  )
  assert_true(
    length(unique(default_order_values)) > 1L,
    "locale fixture does not exercise locale-sensitive default ordering"
  )
  locale_hash_values <- unlist(lapply(locale_hashes, `[[`, "hashes"),
                               use.names = FALSE)
  assert_true(
    length(unique(locale_hash_values)) == 1L,
    "probe-map hash depends on LC_COLLATE"
  )
})

run_test("frozen application never relearns from validation expression", {
  train <- rbind(
    "A_LOW_ID" = c(4, 4, 4, 4),
    "A_HIGH_ID" = c(1, 1, 1, 1),
    "B_ONLY" = c(2, 3, 2, 3)
  )
  colnames(train) <- paste0("T", seq_len(ncol(train)))
  ann <- data.frame(
    PROBEID = rownames(train), SYMBOL = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )
  map <- learn_unique_highest_probe_map(train, ann)
  valid <- train[, 1:2, drop = FALSE]
  colnames(valid) <- c("V1", "V2")
  valid["A_HIGH_ID", ] <- 1000
  frozen <- apply_frozen_probe_map(valid, map)
  assert_equal(
    as.numeric(frozen["A", ]), as.numeric(valid["A_LOW_ID", ]),
    message = "validation expression changed the frozen probe choice"
  )
  expect_error(
    apply_frozen_probe_map(valid[setdiff(rownames(valid), "A_LOW_ID"), ], map),
    "missing frozen probes"
  )
})

run_test("strict nested probe run has complete auditable OOF coverage", {
  dat <- make_probe_survival()
  params <- probe_params()
  nested <- run_probe_nested_cv(
    dat$probe_expr, dat$time, dat$status, dat$sample_ids, params,
    analysis_key = "probe-synthetic-contract", probe_annotation = dat$annotation
  )
  validate_probe_nested_result(
    nested, dat$probe_expr, dat$time, dat$status, dat$sample_ids, params,
    "probe-synthetic-contract", dat$annotation
  )
  assert_true(
    nrow(nested$oof) == length(dat$sample_ids) * params$outer_repeats,
    "OOF row count is incomplete"
  )
  assert_true(!anyDuplicated(paste(nested$oof$repeat_id, nested$oof$sample)),
              "a patient has multiple OOF rows in one repeat")
  assert_true(
    nrow(nested$fold_performance) ==
      params$outer_folds * params$outer_repeats,
    "outer model count is wrong"
  )
  artifacts$dat <- dat
  artifacts$params <- params
  artifacts$nested <- nested
})

run_test("outer and inner probe frequencies conserve every split", {
  nested <- artifacts$nested
  params <- artifacts$params
  outer_denominator <- params$outer_repeats * params$outer_folds
  inner_denominator <- outer_denominator * params$inner_folds
  outer_totals <- stats::aggregate(
    selected_splits ~ SYMBOL, nested$outer_probe_frequency, sum
  )
  inner_totals <- stats::aggregate(
    selected_splits ~ SYMBOL, nested$inner_probe_frequency, sum
  )
  assert_true(all(outer_totals$selected_splits == outer_denominator),
              "outer probe frequencies do not sum to one per model")
  assert_true(all(inner_totals$selected_splits == inner_denominator),
              "inner probe frequencies do not sum to one per split")
  assert_true(all(abs(
    nested$outer_probe_frequency$frequency -
      nested$outer_probe_frequency$selected_splits / outer_denominator
  ) < 1e-12), "outer frequencies use the wrong denominator")
  assert_true(all(abs(
    nested$inner_probe_frequency$frequency -
      nested$inner_probe_frequency$selected_splits / inner_denominator
  ) < 1e-12), "inner frequencies use the wrong denominator")
})

run_test("identical seed reproduces OOF maps hashes and full model", {
  dat <- artifacts$dat
  params <- artifacts$params
  first <- artifacts$nested
  second <- run_probe_nested_cv(
    dat$probe_expr, dat$time, dat$status, dat$sample_ids, params,
    analysis_key = "probe-synthetic-contract", probe_annotation = dat$annotation
  )
  oof_cols <- c("repeat_id", "outer_fold", "sample")
  first_oof <- ordered_rows(first$oof, oof_cols)
  second_oof <- ordered_rows(second$oof, oof_cols)
  assert_equal(first_oof, second_oof, tolerance = 1e-9,
               message = "OOF predictions are not reproducible")
  map_cols <- c("repeat_id", "outer_fold", "SYMBOL", "PROBEID")
  assert_equal(
    ordered_rows(first$outer_probe_maps, map_cols),
    ordered_rows(second$outer_probe_maps, map_cols),
    message = "outer maps are not reproducible"
  )
  hash_cols <- c("repeat_id", "outer_fold", "inner_fold")
  assert_equal(
    ordered_rows(first$inner_probe_map_hashes, hash_cols),
    ordered_rows(second$inner_probe_map_hashes, hash_cols),
    message = "inner map hashes are not reproducible"
  )
  assert_equal(first$outer_probe_frequency, second$outer_probe_frequency,
               message = "outer probe frequencies are not reproducible")
  assert_equal(first$inner_probe_frequency, second$inner_probe_frequency,
               message = "inner probe frequencies are not reproducible")
  assert_true(
    identical(
      first$full_development_model$model_fingerprint,
      second$full_development_model$model_fingerprint
    ),
    "full-development fingerprint is not reproducible"
  )
  artifacts$nested_repeat <- second
})

run_test("inner-validation expression cannot alter its training map or models", {
  dat <- artifacts$dat
  params <- artifacts$params
  baseline <- artifacts$nested$full_development_model
  target_fold <- 1L
  valid_samples <- baseline$inner_fold_assignments$sample[
    baseline$inner_fold_assignments$inner_fold == target_fold
  ]
  changed_expr <- dat$probe_expr
  changed_expr["G01_P2", valid_samples] <-
    changed_expr["G01_P2", valid_samples] + 50
  changed <- fit_probe_signature_pipeline(
    changed_expr, dat$time, dat$status,
    seed = baseline$seed, params = params, probe_annotation = dat$annotation
  )
  base_hash <- baseline$inner_probe_map_hashes$probe_map_hash[
    baseline$inner_probe_map_hashes$inner_fold == target_fold
  ]
  changed_hash <- changed$inner_probe_map_hashes$probe_map_hash[
    changed$inner_probe_map_hashes$inner_fold == target_fold
  ]
  assert_true(identical(base_hash, changed_hash),
              "inner-validation values altered the training-only probe map")
  base_models <- baseline$inner_fold_performance[
    baseline$inner_fold_performance$inner_fold == target_fold,
    c("lambda_ratio", "model_fingerprint"), drop = FALSE
  ]
  changed_models <- changed$inner_fold_performance[
    changed$inner_fold_performance$inner_fold == target_fold,
    c("lambda_ratio", "model_fingerprint"), drop = FALSE
  ]
  base_models <- base_models[order(base_models$lambda_ratio), , drop = FALSE]
  changed_models <- changed_models[order(changed_models$lambda_ratio), , drop = FALSE]
  rownames(base_models) <- rownames(changed_models) <- NULL
  assert_equal(base_models, changed_models,
               message = "inner-validation values altered candidate models")
})

run_test("outer-test expression cannot alter its training map or model", {
  dat <- artifacts$dat
  params <- artifacts$params
  baseline <- artifacts$nested
  target_repeat <- 1L
  target_fold <- 1L
  test_samples <- baseline$outer_fold_assignments$sample[
    baseline$outer_fold_assignments$repeat_id == target_repeat &
      baseline$outer_fold_assignments$outer_fold == target_fold
  ]
  changed_expr <- dat$probe_expr
  changed_expr["G01_P2", test_samples] <-
    changed_expr["G01_P2", test_samples] + 50
  changed <- run_probe_nested_cv(
    changed_expr, dat$time, dat$status, dat$sample_ids, params,
    analysis_key = "probe-outer-sentinel", probe_annotation = dat$annotation
  )
  base_map <- baseline$outer_probe_maps[
    baseline$outer_probe_maps$repeat_id == target_repeat &
      baseline$outer_probe_maps$outer_fold == target_fold,
    c("SYMBOL", "PROBEID", "training_probe_mean", "training_probe_mean_hex",
      "n_probes_gene",
      "map_hash", "training_sample_hash", "training_n"),
    drop = FALSE
  ]
  changed_map <- changed$outer_probe_maps[
    changed$outer_probe_maps$repeat_id == target_repeat &
      changed$outer_probe_maps$outer_fold == target_fold,
    names(base_map), drop = FALSE
  ]
  base_map <- base_map[order(base_map$SYMBOL), , drop = FALSE]
  changed_map <- changed_map[order(changed_map$SYMBOL), , drop = FALSE]
  rownames(base_map) <- rownames(changed_map) <- NULL
  assert_equal(base_map, changed_map,
               message = "outer-test values altered the outer-training map")
  base_fp <- baseline$model_fingerprints$model_fingerprint[
    baseline$model_fingerprints$repeat_id == target_repeat &
      baseline$model_fingerprints$outer_fold == target_fold
  ]
  changed_fp <- changed$model_fingerprints$model_fingerprint[
    changed$model_fingerprints$repeat_id == target_repeat &
      changed$model_fingerprints$outer_fold == target_fold
  ]
  assert_true(identical(base_fp, changed_fp),
              "outer-test values altered the outer-training model")
})

run_test("validator rejects probe-map and fingerprint tampering", {
  dat <- artifacts$dat
  params <- artifacts$params
  baseline <- artifacts$nested
  validate <- function(x) validate_probe_nested_result(
    x, dat$probe_expr, dat$time, dat$status, dat$sample_ids, params,
    "probe-synthetic-contract", dat$annotation
  )

  wrong_probe <- baseline
  wrong_probe$outer_probe_maps$PROBEID[[1L]] <- "TAMPERED_PROBE"
  expect_error(validate(wrong_probe), "hash")

  wrong_map_hash <- baseline
  first_group <- wrong_map_hash$outer_probe_maps$repeat_id == 1L &
    wrong_map_hash$outer_probe_maps$outer_fold == 1L
  wrong_map_hash$outer_probe_maps$map_hash[first_group] <-
    paste(rep("0", 64L), collapse = "")
  expect_error(validate(wrong_map_hash), "hash")

  wrong_model_hash <- baseline
  wrong_model_hash$model_fingerprints$model_fingerprint[[1L]] <-
    paste(rep("f", 64L), collapse = "")
  expect_error(validate(wrong_model_hash), "fingerprint")

  coordinated <- baseline
  first_group <- coordinated$outer_probe_maps$repeat_id == 1L &
    coordinated$outer_probe_maps$outer_fold == 1L
  victim_row <- which(first_group)[[1L]]
  victim_symbol <- coordinated$outer_probe_maps$SYMBOL[[victim_row]]
  original_probe <- coordinated$outer_probe_maps$PROBEID[[victim_row]]
  alternative_probe <- setdiff(
    dat$annotation$PROBEID[dat$annotation$SYMBOL == victim_symbol],
    original_probe
  )[[1L]]
  outer_test_samples <- coordinated$outer_fold_assignments$sample[
    coordinated$outer_fold_assignments$repeat_id == 1L &
      coordinated$outer_fold_assignments$outer_fold == 1L
  ]
  outer_train_samples <- setdiff(dat$sample_ids, outer_test_samples)
  coordinated$outer_probe_maps$PROBEID[[victim_row]] <- alternative_probe
  coordinated$outer_probe_maps$training_probe_mean[[victim_row]] <- mean(
    dat$probe_expr[alternative_probe, outer_train_samples]
  )
  coordinated$outer_probe_maps$training_probe_mean_hex[[victim_row]] <-
    paste0(
      "hex:",
      sprintf(
        "%a", coordinated$outer_probe_maps$training_probe_mean[[victim_row]]
      )
    )
  coordinated_hash <- hash_probe_map(
    coordinated$outer_probe_maps[first_group, , drop = FALSE]
  )
  coordinated$outer_probe_maps$map_hash[first_group] <- coordinated_hash
  coordinated_fp_row <- coordinated$model_fingerprints$repeat_id == 1L &
    coordinated$model_fingerprints$outer_fold == 1L
  coordinated$model_fingerprints$probe_map_hash[coordinated_fp_row] <-
    coordinated_hash
  coordinated_model_hash <- .combine_probe_model_fingerprint(
    coordinated$model_fingerprints$gene_model_fingerprint[coordinated_fp_row],
    coordinated_hash
  )
  coordinated$model_fingerprints$model_fingerprint[coordinated_fp_row] <-
    coordinated_model_hash
  coordinated_fold_row <- coordinated$fold_performance$repeat_id == 1L &
    coordinated$fold_performance$outer_fold == 1L
  coordinated$fold_performance$model_fingerprint[coordinated_fold_row] <-
    coordinated_model_hash
  if ("probe_map_hash" %in% names(coordinated$fold_performance)) {
    coordinated$fold_performance$probe_map_hash[coordinated_fold_row] <-
      coordinated_hash
  }
  expect_error(validate(coordinated), "source expression|training partition")

  coordinated_full_inner <- baseline
  full_model <- coordinated_full_inner$full_development_model
  target_inner_fold <- full_model$inner_probe_map_hashes$inner_fold[[1L]]
  coordinated_inner_hash <- paste(rep("b", 64L), collapse = "")
  full_model$inner_probe_map_hashes$probe_map_hash[[1L]] <-
    coordinated_inner_hash
  inner_tuning_idx <- full_model$inner_fold_performance$inner_fold ==
    target_inner_fold
  full_model$inner_fold_performance$probe_map_hash[inner_tuning_idx] <-
    coordinated_inner_hash
  full_model$inner_fold_performance$model_fingerprint[inner_tuning_idx] <-
    mapply(
      .combine_probe_model_fingerprint,
      full_model$inner_fold_performance$gene_model_fingerprint[
        inner_tuning_idx
      ],
      full_model$inner_fold_performance$probe_map_hash[inner_tuning_idx],
      USE.NAMES = FALSE
    )
  frequency <- full_model$inner_probe_frequency
  victim_frequency_row <- 1L
  victim_frequency_symbol <- frequency$SYMBOL[[victim_frequency_row]]
  victim_frequency_probe <- frequency$PROBEID[[victim_frequency_row]]
  alternative_frequency_probe <- setdiff(
    dat$annotation$PROBEID[
      dat$annotation$SYMBOL == victim_frequency_symbol
    ],
    victim_frequency_probe
  )[[1L]]
  frequency$selected_splits[[victim_frequency_row]] <-
    frequency$selected_splits[[victim_frequency_row]] - 1L
  frequency$frequency[[victim_frequency_row]] <-
    frequency$selected_splits[[victim_frequency_row]] /
    frequency$denominator[[victim_frequency_row]]
  alternative_row <- frequency[victim_frequency_row, , drop = FALSE]
  alternative_row$PROBEID <- alternative_frequency_probe
  alternative_row$selected_splits <- 1L
  alternative_row$frequency <- 1 / alternative_row$denominator
  frequency <- rbind(frequency, alternative_row)
  frequency <- frequency[frequency$selected_splits > 0L, , drop = FALSE]
  full_model$inner_probe_frequency <- frequency
  coordinated_full_inner$full_development_model <- full_model
  expect_error(
    validate(coordinated_full_inner),
    "Full-development inner probe-map|source-derived"
  )

  wrong_full_hash <- baseline
  wrong_full_hash$full_development_model$model_fingerprint <-
    paste(rep("a", 64L), collapse = "")
  expect_error(validate(wrong_full_hash), "fingerprint")
})

summary_df <- do.call(rbind, results)
cat("\nHighest-mean nested probe test summary\n")
print(summary_df, row.names = FALSE)
cat(sprintf("\nPASS=%d FAIL=%d\n",
            sum(summary_df$status == "PASS"),
            sum(summary_df$status == "FAIL")))
if (any(summary_df$status != "PASS")) quit(save = "no", status = 1L)
quit(save = "no", status = 0L)
