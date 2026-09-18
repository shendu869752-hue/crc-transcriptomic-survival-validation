#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript")
this_file <- gsub("\\\\", "/", sub("^--file=", "", script_arg[[1L]]))
test_dir <- dirname(this_file)
validator_file <- file.path(test_dir, "validate_completed_run.R")
if (!file.exists(validator_file)) stop("Cannot locate validate_completed_run.R")
source(validator_file)

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, quote = TRUE, na = "")
  invisible(path)
}

write_dput <- function(x, path) {
  writeLines(capture.output(dput(x)), path, useBytes = TRUE)
  invisible(path)
}

make_mapping_audit <- function(genes, annotation_version) {
  stopifnot(length(genes) == 2L)
  data.frame(
    required_gene = genes,
    required_entrez_id = c("101", "202"),
    source_row_id = c(genes[[1L]], "OLD2"),
    resolution_method = c("exact_SYMBOL", "unique_reverse_unique_ENTREZ_alias"),
    exact_symbol_present = c(TRUE, FALSE),
    expression_alias_candidates = c("", "OLD2"),
    expression_alias_candidate_count = c(0L, 1L),
    reverse_unique = c(NA, TRUE),
    symbol_consistent = c(NA, TRUE),
    entrez_consistent = c(TRUE, TRUE),
    reverse_symbols = c("", genes[[2L]]),
    reverse_entrez_ids = c("", "202"),
    reverse_symbol_count = c(NA_integer_, 1L),
    reverse_entrez_count = c(NA_integer_, 1L),
    annotation_package = "org.Hs.eg.db",
    annotation_version = annotation_version,
    outcome_or_expression_values_used_for_resolution = FALSE,
    stringsAsFactors = FALSE
  )
}

make_meta_rows <- function(cohort_inputs, cohorts, inference_df) {
  data.frame(
    endpoint = "OS",
    cohorts = cohorts,
    pooled_HR = 1.25,
    lower95 = c(1.05, 1.10),
    upper95 = c(1.49, 1.42),
    p = c(0.03, 0.002),
    prediction_lower95 = c(0.90, 0.95),
    prediction_upper95 = c(1.75, 1.65),
    Q = 1.2,
    Q_df = cohorts - 1L,
    Q_p = 0.55,
    I2 = 0,
    tau2 = 0,
    test_statistic = c(3.1, 3.1),
    inference_df = c(inference_df, NA_real_),
    inference = c("Knapp-Hartung/t", "normal/Wald"),
    primary = c(TRUE, FALSE),
    method = "REML random-effects inverse-variance meta-analysis",
    cohort_inputs = cohort_inputs,
    analysis = "synthetic completed-run fixture",
    stringsAsFactors = FALSE
  )
}

make_external_performance <- function(scores) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required for the completed-run fixture")
  }
  rows <- lapply(split(scores, scores$cohort), function(z) {
    fit <- survival::coxph(
      survival::Surv(time, status) ~ risk_sd,
      data = z, ties = "efron", x = TRUE, y = TRUE, na.action = stats::na.fail
    )
    fit_summary <- summary(fit)
    ci <- stats::confint(fit)
    km <- survival::survdiff(survival::Surv(time, status) ~ group, data = z)
    c_obj <- survival::concordance(
      survival::Surv(time, status) ~ risk_sd,
      data = z, reverse = TRUE, timewt = "n"
    )
    data.frame(
      cohort = unique(z$cohort), endpoint = unique(z$endpoint),
      mapping_method = unique(z$mapping_method), n = nrow(z),
      events = sum(z$status), HR_per_SD = unname(exp(stats::coef(fit))),
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

build_fixture <- function(parent_dir) {
  run_key <- paste(rep("b", 64L), collapse = "")
  analysis_key <- paste(rep("c", 64L), collapse = "")
  strict_analysis_key <- paste(rep("e", 64L), collapse = "")
  pipeline_version <- "synthetic-primary-nested-v1"
  strict_pipeline_version <- "synthetic-highest-mean-nested-v1"
  annotation_version <- "3.23.1"
  run_dir <- file.path(parent_dir, run_key)
  results_dir <- file.path(run_dir, "results")
  snapshot_dir <- file.path(run_dir, "code_snapshot")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(snapshot_dir, recursive = TRUE, showWarnings = FALSE)

  samples <- sprintf("S%03d", seq_len(573L))
  outer_fold <- ((seq_along(samples) - 1L) %% 5L) + 1L
  outer <- data.frame(
    sample = samples, repeat_id = 1L, outer_fold = outer_fold,
    analysis_key = analysis_key, stringsAsFactors = FALSE
  )
  inner <- do.call(rbind, lapply(seq_len(5L), function(fold) {
    train <- samples[outer_fold != fold]
    data.frame(
      sample = train,
      inner_fold = ((seq_along(train) - 1L) %% 5L) + 1L,
      repeat_id = 1L,
      outer_fold = fold,
      analysis_key = analysis_key,
      stringsAsFactors = FALSE
    )
  }))
  oof <- data.frame(
    sample = samples, repeat_id = 1L, outer_fold = outer_fold,
    time = seq_len(573L) + 10, status = rep(c(0L, 1L), length.out = 573L),
    lp_raw = seq(-1, 1, length.out = 573L),
    lp_oof = seq(-1, 1, length.out = 573L),
    analysis_key = analysis_key, stringsAsFactors = FALSE
  )
  fold_performance <- data.frame(
    repeat_id = 1L, outer_fold = seq_len(5L), c_index = 0.6,
    analysis_key = analysis_key, stringsAsFactors = FALSE
  )
  repeat_performance <- data.frame(
    repeat_id = 1L, n = 573L, events = sum(oof$status), c_index = 0.6,
    analysis_key = analysis_key, stringsAsFactors = FALSE
  )
  nested_summary <- data.frame(
    pipeline_version = pipeline_version, analysis_key = analysis_key,
    repeats = 1L, outer_folds = 5L, inner_folds = 5L, completed_folds = 5L,
    median_repeat_oof_c_index = 0.6, repeat_oof_c_index_q1 = 0.6,
    repeat_oof_c_index_q3 = 0.6, min_repeat_oof_c_index = 0.6,
    max_repeat_oof_c_index = 0.6, median_within_repeat_jaccard = 0.5,
    jaccard_q1 = 0.4, jaccard_q3 = 0.6, stringsAsFactors = FALSE
  )
  failure_log <- data.frame(
    repeat_id = integer(), outer_fold = integer(), stage = character(),
    message = character(), stringsAsFactors = FALSE
  )

  write_csv(nested_summary, file.path(results_dir, "nested_cv_summary.csv"))
  write_csv(oof, file.path(results_dir, "nested_cv_oof_predictions.csv"))
  write_csv(outer, file.path(results_dir, "nested_cv_outer_fold_assignments.csv"))
  write_csv(inner, file.path(results_dir, "nested_cv_inner_fold_assignments.csv"))
  write_csv(fold_performance, file.path(results_dir, "nested_cv_fold_performance.csv"))
  write_csv(repeat_performance, file.path(results_dir, "nested_cv_repeat_performance.csv"))
  write_csv(failure_log, file.path(results_dir, "nested_cv_failure_log.csv"))

  # The strict highest-mean sensitivity has its own content identity and its
  # own folds/OOF audit trail. Shift the synthetic folds so a validator that
  # accidentally borrows the primary assignments will fail the map hashes.
  strict_outer_fold <- (outer_fold %% 5L) + 1L
  strict_outer <- data.frame(
    sample = samples, repeat_id = 1L, outer_fold = strict_outer_fold,
    analysis_key = strict_analysis_key, stringsAsFactors = FALSE
  )
  strict_inner <- do.call(rbind, lapply(seq_len(5L), function(fold) {
    train <- samples[strict_outer_fold != fold]
    data.frame(
      sample = train,
      inner_fold = (seq_along(train) %% 5L) + 1L,
      repeat_id = 1L,
      outer_fold = fold,
      analysis_key = strict_analysis_key,
      stringsAsFactors = FALSE
    )
  }))
  strict_oof <- transform(
    oof,
    outer_fold = strict_outer_fold,
    lp_raw = lp_raw + 0.01,
    lp_oof = lp_oof + 0.01,
    analysis_key = strict_analysis_key
  )
  strict_fold_performance <- transform(
    fold_performance, c_index = 0.59, analysis_key = strict_analysis_key
  )
  strict_repeat_performance <- transform(
    repeat_performance, c_index = 0.59, analysis_key = strict_analysis_key
  )
  strict_summary <- transform(
    nested_summary,
    pipeline_version = strict_pipeline_version,
    analysis_key = strict_analysis_key,
    median_repeat_oof_c_index = 0.59,
    repeat_oof_c_index_q1 = 0.59,
    repeat_oof_c_index_q3 = 0.59,
    min_repeat_oof_c_index = 0.59,
    max_repeat_oof_c_index = 0.59
  )
  write_csv(
    strict_summary,
    file.path(results_dir, "highest_mean_nested_cv_summary.csv")
  )
  write_csv(
    strict_oof,
    file.path(results_dir, "highest_mean_nested_cv_oof_predictions.csv")
  )
  write_csv(
    strict_outer,
    file.path(results_dir, "highest_mean_nested_cv_outer_fold_assignments.csv")
  )
  write_csv(
    strict_inner,
    file.path(results_dir, "highest_mean_nested_cv_inner_fold_assignments.csv")
  )
  write_csv(
    strict_fold_performance,
    file.path(results_dir, "highest_mean_nested_cv_fold_performance.csv")
  )
  write_csv(
    strict_repeat_performance,
    file.path(results_dir, "highest_mean_nested_cv_repeat_performance.csv")
  )
  write_csv(
    failure_log,
    file.path(results_dir, "highest_mean_nested_cv_failure_log.csv")
  )

  input_manifest <- data.frame(
    source_id = c("development", "external"),
    relative_path = c("inputs/development.csv", "inputs/external.csv"),
    bytes = c(123, 456),
    sha256 = c(paste(rep("a", 64L), collapse = ""),
               paste(rep("d", 64L), collapse = "")),
    stringsAsFactors = FALSE
  )
  write_csv(input_manifest, file.path(results_dir, "input_manifest.csv"))

  code_paths <- c(
    core_analysis = file.path(snapshot_dir, "01_core_reanalysis_v3.R"),
    local_launcher = file.path(snapshot_dir, "run_v3_local.ps1")
  )
  writeLines(c("x <- 1L", "stopifnot(x == 1L)"), code_paths[[1L]], useBytes = TRUE)
  writeLines("Write-Output 'fixture'", code_paths[[2L]], useBytes = TRUE)
  code_hashes <- vapply(code_paths, qa_sha256_file, character(1))
  code_rel <- gsub("\\\\", "/", file.path("code_snapshot", basename(code_paths)))
  code_manifest <- rbind(
    data.frame(category = "code_sha256", name = names(code_paths),
               value = unname(code_hashes), stringsAsFactors = FALSE),
    data.frame(category = "code_snapshot_relative_path", name = names(code_paths),
               value = code_rel, stringsAsFactors = FALSE),
    data.frame(category = "code_snapshot_bytes", name = names(code_paths),
               value = as.character(file.info(code_paths)$size), stringsAsFactors = FALSE),
    data.frame(category = "runtime", name = c("R.version.string", "R.platform"),
               value = c(R.version.string, R.version$platform), stringsAsFactors = FALSE),
    data.frame(category = "package_version", name = "org.Hs.eg.db",
               value = annotation_version, stringsAsFactors = FALSE)
  )
  write_csv(code_manifest, file.path(results_dir, "code_manifest.csv"))

  params <- list(
    output_contract_version = "crc9-completed-run-fixture-v1",
    pipeline_version = pipeline_version,
    nested_pipeline = list(
      pipeline_version = pipeline_version,
      outer_folds = 5L, outer_repeats = 1L, inner_folds = 5L,
      base_seed = 20240916L, sd_cutoff = 0.2, univ_p_cutoff = 1e-4,
      min_candidates = 5L, candidate_fallback_n = 50L, coef_cutoff = 0.05,
      min_selected = 2L, selected_fallback_n = 2L, max_selected = 30L,
      lambda_ratio_grid = c(1, 0.1), lambda_rule = "one_se_larger_penalty",
      ties = "efron", development_probe_aggregation = "unique_mean"
    ),
    highest_mean_nested_pipeline = list(
      pipeline_version = strict_pipeline_version,
      outer_folds = 5L, outer_repeats = 1L, inner_folds = 5L,
      base_seed = 20240916L, sd_cutoff = 0.2, univ_p_cutoff = 1e-4,
      min_candidates = 5L, candidate_fallback_n = 50L, coef_cutoff = 0.05,
      min_selected = 2L, selected_fallback_n = 2L, max_selected = 30L,
      lambda_ratio_grid = c(1, 0.1), lambda_rule = "one_se_larger_penalty",
      ties = "efron",
      development_probe_aggregation =
        "unique_highest_mean_training_partition_learned",
      external_probe_application = "full_development_frozen"
    ),
    bootstrap = list(
      replicates = 100L, percentile_probabilities = c(0.025, 0.975),
      cohort_seed_base = 20240917L, ensemble_seed = 20240930L,
      minimum_valid_fraction = 0.95,
      estimand = "fixed-score patient percentile bootstrap"
    ),
    probe_mapping = list(
      primary = "unique_mean", sensitivity = "unique_highest_mean",
      eligibility = "exactly_one_nonempty_SYMBOL",
      primary_multi_probe_rule = "arithmetic_mean_on_log2_scale",
      sensitivity_tie_break = "PROBEID_ascending", annotation = "hgu133plus2.db",
      external_sensitivity_scope =
        "full_development_probe_map_frozen_for_external_scoring"
    ),
    gene_identity_resolution = list(
      scope = "TCGA_COAD_signature_rows_only",
      rule = "exact_SYMBOL_then_unique_reverse_unique_ENTREZ_alias",
      exact_symbol_priority = TRUE,
      alias_requirements = c("one_expression_row_candidate", "one_reverse_current_SYMBOL"),
      annotation = "org.Hs.eg.db", annotation_version = annotation_version,
      outcome_blind = TRUE, expression_values_used = FALSE,
      development_candidate_space_filtered_by_TCGA = FALSE
    ),
    scoring = list(scope = "cohort_adaptive_standardization",
                   risk_group_cutpoint = "cohort_median"),
    time_dependent_auc = list(years = c(1, 3, 5), min_events = 10L,
                              min_at_risk = 20L, weighting = "marginal"),
    survival = list(
      geo_month_to_day = 30.4375, year_days = 365.25, cox_ties = "efron",
      spline_df = 3L, clinical_min_n = 50L, clinical_min_events = 20L,
      os_meta_method = "REML", os_meta_primary_test = "knha"
    ),
    execution = list(force_nested = TRUE, entrypoint = "run_v3_local.ps1")
  )
  write_dput(params, file.path(results_dir, "run_parameters.txt"))

  genes <- c("GENE1", "GENE2")
  primary_coefficients <- data.frame(
    gene = genes, coef = c(0.2, -0.1), HR = exp(c(0.2, -0.1)),
    stringsAsFactors = FALSE
  )
  write_csv(primary_coefficients, file.path(results_dir, "final_model_coefficients.csv"))
  write_csv(
    make_mapping_audit(genes, annotation_version),
    file.path(results_dir, "TCGA_COAD_signature_gene_mapping.csv")
  )

  meta_cohorts <- sort(c("TCGA-COAD", "GSE17536", "GSE17537"))
  meta_inputs <- data.frame(
    cohort = meta_cohorts, endpoint = "OS", mapping_method =
      c("unique_mean", "unique_mean", "gene_level_input"),
    n = c(100L, 110L, 120L), events = c(30L, 35L, 40L),
    HR_per_SD = c(1.2, 1.3, 1.25), lower95 = c(1.01, 1.05, 1.02),
    upper95 = c(1.43, 1.61, 1.53), cox_p = c(0.04, 0.02, 0.03),
    yi = log(c(1.2, 1.3, 1.25)), sei = c(0.09, 0.11, 0.10),
    effect_scale = "log hazard ratio per 1-SD risk score",
    sei_method = "derived from cohort-level 95% Wald CI",
    stringsAsFactors = FALSE
  )
  write_csv(meta_inputs, file.path(results_dir, "os_validation_meta_inputs.csv"))
  meta <- make_meta_rows(paste(meta_cohorts, collapse = ";"), 3L, 2)
  write_csv(meta, file.path(results_dir, "os_validation_meta_analysis.csv"))
  loo <- do.call(rbind, lapply(meta_cohorts, function(drop) {
    cbind(
      dropped_cohort = drop,
      make_meta_rows(paste(setdiff(meta_cohorts, drop), collapse = ";"), 2L, 1),
      stringsAsFactors = FALSE
    )
  }))
  write_csv(loo, file.path(results_dir, "os_validation_meta_leave_one_out.csv"))

  clinical_cohorts <- c(
    "GSE39582 (training)", "TCGA-COAD", "GSE14333", "GSE17536", "GSE17537"
  )
  source_n <- c(573L, 279L, 226L, 177L, 55L)
  source_events <- c(190L, 69L, 50L, 73L, 20L)
  complete_n <- c(500L, 200L, 40L, 150L, 30L)
  complete_events <- c(170L, 55L, 15L, 60L, 10L)
  eligible <- c(TRUE, TRUE, FALSE, TRUE, FALSE)
  eligibility <- data.frame(
    cohort = clinical_cohorts, endpoint = c("OS", "OS", "DFS", "OS", "OS"),
    mapping_method = c("unique_mean", "gene_level_input", rep("unique_mean", 3L)),
    source_n = source_n, source_events = source_events,
    complete_case_n = complete_n, complete_case_events = complete_events,
    excluded_from_adjustment_n = source_n - complete_n,
    excluded_events = source_events - complete_events,
    minimum_n_required = 50L, minimum_events_required = 20L,
    eligible = eligible,
    reason = ifelse(eligible, "eligible", "complete_case_threshold_not_met"),
    stringsAsFactors = FALSE
  )
  write_csv(eligibility, file.path(results_dir, "clinical_adjustment_eligibility.csv"))
  write_csv(
    data.frame(term = genes, chisq = c(0.2, 0.3), df = 1L, p = c(0.65, 0.58)),
    file.path(results_dir, "proportional_hazards_test.csv")
  )
  write_csv(
    data.frame(
      cohort = clinical_cohorts, endpoint = c("OS", "OS", "DFS", "OS", "OS"),
      mapping_method = c("unique_mean", "gene_level_input", rep("unique_mean", 3L)),
      n = source_n, events = source_events, chisq = 0.2, p = 0.65
    ),
    file.path(results_dir, "risk_score_proportional_hazards_test.csv")
  )
  write_csv(
    data.frame(cohort = clinical_cohorts[eligible], endpoint = "OS", term = "risk_sd",
               chisq = 0.2, p = 0.65),
    file.path(results_dir, "clinical_model_proportional_hazards_test.csv")
  )
  write_csv(
    data.frame(
      cohort = clinical_cohorts, endpoint = c("OS", "OS", "DFS", "OS", "OS"),
      n = source_n, events = source_events, linear_loglik = -100,
      spline_loglik = -99.5, p_nonlinearity = 0.5,
      method = "LRT: linear Cox vs pspline(df=3)"
    ),
    file.path(results_dir, "risk_score_linearity_test.csv")
  )

  base_map <- data.frame(
    SYMBOL = genes, PROBEID = c("P1", "P2"),
    training_probe_mean = c(5.5, 6.5),
    training_probe_mean_hex = paste0("hex:", sprintf("%a", c(5.5, 6.5))),
    n_probes_gene = c(2L, 1L),
    stringsAsFactors = FALSE
  )
  outer_maps <- do.call(rbind, lapply(seq_len(5L), function(fold) {
    z <- base_map
    z$training_probe_mean <- z$training_probe_mean + fold / 100
    z$training_probe_mean_hex <- paste0(
      "hex:", sprintf("%a", z$training_probe_mean)
    )
    train <- samples[strict_outer_fold != fold]
    transform(
      z, repeat_id = 1L, outer_fold = fold, map_hash = qa_hash_probe_map(z),
      training_sample_hash = qa_hash_sample_ids(train), training_n = length(train),
      analysis_key = strict_analysis_key
    )
  }))
  write_csv(outer_maps, file.path(results_dir, "highest_mean_outer_probe_maps.csv"))

  inner_hash_rows <- list()
  counter <- 0L
  for (fold in seq_len(5L)) {
    current <- strict_inner[strict_inner$outer_fold == fold, , drop = FALSE]
    for (inner_fold in seq_len(5L)) {
      counter <- counter + 1L
      train <- current$sample[current$inner_fold != inner_fold]
      inner_hash_rows[[counter]] <- data.frame(
        repeat_id = 1L, outer_fold = fold, inner_fold = inner_fold,
        probe_map_hash = qa_hash_probe_map(base_map),
        training_sample_hash = qa_hash_sample_ids(train), training_n = length(train),
        analysis_key = strict_analysis_key, stringsAsFactors = FALSE
      )
    }
  }
  write_csv(
    do.call(rbind, inner_hash_rows),
    file.path(results_dir, "highest_mean_inner_probe_map_hashes.csv")
  )
  outer_frequency <- transform(
    base_map[, c("SYMBOL", "PROBEID")], selected_splits = 5L,
    denominator = 5L, frequency = 1, analysis_key = strict_analysis_key
  )
  inner_frequency <- transform(
    base_map[, c("SYMBOL", "PROBEID")], selected_splits = 25L,
    denominator = 25L, frequency = 1, analysis_key = strict_analysis_key
  )
  write_csv(
    outer_frequency, file.path(results_dir, "highest_mean_outer_probe_frequency.csv")
  )
  write_csv(
    inner_frequency, file.path(results_dir, "highest_mean_inner_probe_frequency.csv")
  )

  full_hash <- qa_hash_probe_map(base_map)
  full_sample_hash <- qa_hash_sample_ids(samples)
  full_map <- transform(
    base_map, map_hash = full_hash, training_sample_hash = full_sample_hash,
    training_n = length(samples), analysis_key = strict_analysis_key
  )
  write_csv(
    full_map, file.path(results_dir, "highest_mean_full_development_probe_map.csv")
  )
  model <- list(
    genes = genes,
    coefs = setNames(c(0.2, -0.1), genes),
    center = setNames(c(5.5, 6.5), genes),
    scale = setNames(c(1.1, 1.2), genes),
    train_lp_mean = 0,
    train_lp_sd = 1,
    lambda_ratio = 0.1,
    lambda = 0.01,
    probe_map = list(
      mapping = base_map, map_hash = full_hash,
      training_sample_hash = full_sample_hash, training_n = length(samples)
    ),
    probe_map_hash = full_hash
  )
  model$gene_model_fingerprint <- qa_hash_gene_model(model)
  model$model_fingerprint <- qa_combine_model_fingerprint(
    model$gene_model_fingerprint, model$probe_map_hash
  )
  saveRDS(
    model, file.path(results_dir, "highest_mean_full_development_model.rds"),
    version = 3
  )
  sensitivity_coefficients <- data.frame(
    gene = genes, coef = unname(model$coefs[genes]), probe_map_hash = full_hash,
    model_fingerprint = model$model_fingerprint, analysis_key = strict_analysis_key,
    stringsAsFactors = FALSE
  )
  write_csv(
    sensitivity_coefficients,
    file.path(results_dir, "highest_mean_final_model_coefficients.csv")
  )
  strict_mapping_audit <- make_mapping_audit(genes, annotation_version)
  strict_mapping_audit$probe_map_hash <- full_hash
  strict_mapping_audit$model_fingerprint <- model$model_fingerprint
  strict_mapping_audit$analysis_key <- strict_analysis_key
  write_csv(
    strict_mapping_audit,
    file.path(results_dir, "highest_mean_TCGA_COAD_signature_gene_mapping.csv")
  )
  external_cohorts <- c("TCGA-COAD", "GSE14333", "GSE17536", "GSE17537")
  scores <- do.call(rbind, lapply(seq_along(external_cohorts), function(i) {
    cohort <- external_cohorts[[i]]
    raw_score <- c(-1.7, -1.1, -0.6, -0.2, 0.15, 0.55, 1.05, 1.65) + i / 100
    data.frame(
      cohort = cohort,
      sample = paste0(cohort, "-", seq_along(raw_score)),
      time = c(90, 250, 170, 330, 120, 290, 210, 380) + i,
      status = c(1L, 0L, 1L, 0L, 0L, 1L, 1L, 1L),
      endpoint = if (cohort == "GSE14333") "DFS" else "OS",
      lp_raw = raw_score, risk_score = raw_score,
      risk_sd = as.numeric((raw_score - mean(raw_score)) / stats::sd(raw_score)),
      group = ifelse(raw_score > stats::median(raw_score), "High", "Low"),
      mapping_method = if (cohort == "TCGA-COAD") {
        "gene_level_input"
      } else {
        "full_development_frozen_unique_highest_mean"
      },
      probe_map_hash = full_hash, model_fingerprint = model$model_fingerprint,
      analysis_key = strict_analysis_key, stringsAsFactors = FALSE
    )
  }))
  write_csv(scores, file.path(results_dir, "highest_mean_external_scores.csv"))

  strict_performance <- make_external_performance(scores)
  strict_performance$probe_map_hash <- full_hash
  strict_performance$model_fingerprint <- model$model_fingerprint
  strict_performance$analysis_key <- strict_analysis_key
  write_csv(
    strict_performance,
    file.path(results_dir, "highest_mean_external_cohort_performance.csv")
  )
  strict_meta_cohorts <- sort(c("TCGA-COAD", "GSE17536", "GSE17537"))
  strict_os <- strict_performance[
    match(strict_meta_cohorts, strict_performance$cohort), , drop = FALSE
  ]
  strict_meta_inputs <- strict_os[, c(
    "cohort", "endpoint", "mapping_method", "n", "events", "HR_per_SD",
    "lower95", "upper95", "cox_p", "probe_map_hash", "model_fingerprint",
    "analysis_key"
  )]
  strict_meta_inputs$yi <- log(strict_meta_inputs$HR_per_SD)
  strict_meta_inputs$sei <- (
    log(strict_meta_inputs$upper95) - log(strict_meta_inputs$lower95)
  ) / (2 * stats::qnorm(0.975))
  strict_meta_inputs$effect_scale <-
    "log hazard ratio per 1-SD cohort-standardized risk score"
  strict_meta_inputs$sei_method <- "derived from cohort-level 95% Wald CI"
  write_csv(
    strict_meta_inputs,
    file.path(results_dir, "highest_mean_os_validation_meta_inputs.csv")
  )
  strict_meta <- make_meta_rows(paste(strict_meta_cohorts, collapse = ";"), 3L, 2)
  strict_meta$probe_map_hash <- full_hash
  strict_meta$model_fingerprint <- model$model_fingerprint
  strict_meta$analysis_key <- strict_analysis_key
  write_csv(
    strict_meta,
    file.path(results_dir, "highest_mean_os_validation_meta_analysis.csv")
  )
  strict_loo <- do.call(rbind, lapply(strict_meta_cohorts, function(drop) {
    z <- cbind(
      dropped_cohort = drop,
      make_meta_rows(
        paste(setdiff(strict_meta_cohorts, drop), collapse = ";"), 2L, 1
      ),
      stringsAsFactors = FALSE
    )
    z$probe_map_hash <- full_hash
    z$model_fingerprint <- model$model_fingerprint
    z$analysis_key <- strict_analysis_key
    z
  }))
  write_csv(
    strict_loo,
    file.path(results_dir, "highest_mean_os_validation_meta_leave_one_out.csv")
  )

  generated_abs <- list.files(
    run_dir, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE
  )
  generated_abs <- generated_abs[!file.info(generated_abs)$isdir]
  root_prefix <- paste0(normalizePath(run_dir, winslash = "/"), "/")
  generated_rel <- substring(
    normalizePath(generated_abs, winslash = "/"), nchar(root_prefix) + 1L
  )
  manifest <- data.frame(
    artifact_id = generated_rel,
    relative_path = generated_rel,
    required = TRUE,
    status = "generated",
    reason = "",
    bytes = as.numeric(file.info(generated_abs)$size),
    sha256 = unname(vapply(generated_abs, qa_sha256_file, character(1))),
    stringsAsFactors = FALSE
  )
  manifest <- manifest[order(manifest$relative_path), , drop = FALSE]
  write_csv(manifest, file.path(results_dir, "output_manifest.csv"))

  status <- data.frame(
    status = "complete", run_key = run_key, nested_analysis_key = analysis_key,
    input_manifest_sha256 = qa_sha256_file(file.path(results_dir, "input_manifest.csv")),
    parameter_manifest_sha256 = qa_sha256_file(file.path(results_dir, "run_parameters.txt")),
    code_manifest_sha256 = qa_sha256_file(file.path(results_dir, "code_manifest.csv")),
    output_manifest_sha256 = qa_sha256_file(file.path(results_dir, "output_manifest.csv")),
    started_at_utc = "2026-09-17T00:00:00Z",
    finished_at_utc = "2026-09-17T00:01:00Z",
    detail = "synthetic completed-run fixture",
    stringsAsFactors = FALSE
  )
  write_csv(status, file.path(results_dir, "run_status.csv"))
  run_dir
}

reseal_output_manifest <- function(run_dir) {
  manifest_path <- file.path(run_dir, "results", "output_manifest.csv")
  manifest <- utils::read.csv(
    manifest_path, stringsAsFactors = FALSE, check.names = FALSE,
    na.strings = "NA"
  )
  for (i in seq_len(nrow(manifest))) {
    if (manifest$status[[i]] == "generated") {
      path <- file.path(run_dir, manifest$relative_path[[i]])
      manifest$bytes[[i]] <- as.numeric(file.info(path)$size)
      manifest$sha256[[i]] <- qa_sha256_file(path)
    }
  }
  write_csv(manifest, manifest_path)
  status_path <- file.path(run_dir, "results", "run_status.csv")
  status <- utils::read.csv(status_path, stringsAsFactors = FALSE, check.names = FALSE)
  status$output_manifest_sha256 <- qa_sha256_file(manifest_path)
  write_csv(status, status_path)
  invisible(run_dir)
}

expect_error <- function(expr, pattern = NULL) {
  error <- NULL
  tryCatch(force(expr), error = function(e) error <<- conditionMessage(e))
  if (is.null(error)) stop("Expected an error, but expression returned normally")
  if (!is.null(pattern) && !grepl(pattern, error, ignore.case = TRUE)) {
    stop("Error did not match /", pattern, "/: ", error)
  }
  invisible(error)
}

results <- list()
run_test <- function(name, code) {
  error <- NULL
  tryCatch(force(code), error = function(e) error <<- conditionMessage(e))
  if (is.null(error)) {
    cat("PASS:", name, "\n")
    results[[length(results) + 1L]] <<- data.frame(
      test = name, status = "PASS", detail = "", stringsAsFactors = FALSE
    )
  } else {
    cat("FAIL:", name, "-", error, "\n")
    results[[length(results) + 1L]] <<- data.frame(
      test = name, status = "FAIL", detail = error, stringsAsFactors = FALSE
    )
  }
  invisible(NULL)
}

run_test("valid completed-run fixture passes function and CLI", {
  parent <- tempfile("crc9-completed-valid-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE, force = TRUE), add = TRUE)
  run_dir <- build_fixture(parent)
  result <- validate_completed_run(run_dir, quiet = TRUE)
  stopifnot(
    result$development_samples == 573L,
    identical(result$strict_analysis_key, paste(rep("e", 64L), collapse = ""))
  )
  cli <- system2(
    file.path(R.home("bin"), "Rscript.exe"),
    c(shQuote(validator_file), shQuote(run_dir)),
    stdout = TRUE, stderr = TRUE
  )
  exit_status <- attr(cli, "status")
  if (is.null(exit_status)) exit_status <- 0L
  if (exit_status != 0L) stop(paste(cli, collapse = "\n"))
  stopifnot(any(grepl("COMPLETED-RUN QA: PASS", cli, fixed = TRUE)))
})

run_test("unsealed artifact tamper breaks manifest verification", {
  parent <- tempfile("crc9-completed-byte-tamper-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE, force = TRUE), add = TRUE)
  run_dir <- build_fixture(parent)
  cat("\n# tampered\n", file = file.path(
    run_dir, "code_snapshot", "01_core_reanalysis_v3.R"
  ), append = TRUE)
  expect_error(validate_completed_run(run_dir, quiet = TRUE), "Manifest (byte|SHA-256) mismatch")
})

run_test("resealed frozen-map tamper still fails semantic QA", {
  parent <- tempfile("crc9-completed-semantic-tamper-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE, force = TRUE), add = TRUE)
  run_dir <- build_fixture(parent)
  score_path <- file.path(run_dir, "results", "highest_mean_external_scores.csv")
  scores <- utils::read.csv(score_path, stringsAsFactors = FALSE, check.names = FALSE)
  scores$probe_map_hash[[1L]] <- paste(rep("d", 64L), collapse = "")
  write_csv(scores, score_path)
  reseal_output_manifest(run_dir)
  expect_error(
    validate_completed_run(run_dir, quiet = TRUE),
    "External scores are not bound"
  )
})

run_test("strict nested analysis identity cannot reuse the primary key", {
  parent <- tempfile("crc9-completed-strict-key-tamper-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE, force = TRUE), add = TRUE)
  run_dir <- build_fixture(parent)
  primary_path <- file.path(run_dir, "results", "nested_cv_summary.csv")
  strict_path <- file.path(run_dir, "results", "highest_mean_nested_cv_summary.csv")
  primary <- utils::read.csv(primary_path, stringsAsFactors = FALSE, check.names = FALSE)
  strict <- utils::read.csv(strict_path, stringsAsFactors = FALSE, check.names = FALSE)
  strict$analysis_key <- primary$analysis_key
  write_csv(strict, strict_path)
  reseal_output_manifest(run_dir)
  expect_error(
    validate_completed_run(run_dir, quiet = TRUE),
    "Strict nested analysis_key must be independent"
  )
})

run_test("strict nested pipeline version cannot reuse the primary version", {
  parent <- tempfile("crc9-completed-strict-version-tamper-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE, force = TRUE), add = TRUE)
  run_dir <- build_fixture(parent)
  primary_path <- file.path(run_dir, "results", "nested_cv_summary.csv")
  strict_path <- file.path(run_dir, "results", "highest_mean_nested_cv_summary.csv")
  primary <- utils::read.csv(primary_path, stringsAsFactors = FALSE, check.names = FALSE)
  strict <- utils::read.csv(strict_path, stringsAsFactors = FALSE, check.names = FALSE)
  strict$pipeline_version <- primary$pipeline_version
  write_csv(strict, strict_path)
  reseal_output_manifest(run_dir)
  expect_error(
    validate_completed_run(run_dir, quiet = TRUE),
    "Strict highest-mean nested CV summary pipeline_version differs"
  )
})

run_test("strict external time must be positive", {
  parent <- tempfile("crc9-completed-strict-outcome-tamper-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE, force = TRUE), add = TRUE)
  run_dir <- build_fixture(parent)
  score_path <- file.path(run_dir, "results", "highest_mean_external_scores.csv")
  scores <- utils::read.csv(score_path, stringsAsFactors = FALSE, check.names = FALSE)
  scores$time[[1L]] <- 0
  write_csv(scores, score_path)
  reseal_output_manifest(run_dir)
  expect_error(validate_completed_run(run_dir, quiet = TRUE), "time must be positive")
})

run_test("strict external status must be binary", {
  parent <- tempfile("crc9-completed-strict-status-tamper-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE, force = TRUE), add = TRUE)
  run_dir <- build_fixture(parent)
  score_path <- file.path(run_dir, "results", "highest_mean_external_scores.csv")
  scores <- utils::read.csv(score_path, stringsAsFactors = FALSE, check.names = FALSE)
  scores$status[[1L]] <- 2L
  write_csv(scores, score_path)
  reseal_output_manifest(run_dir)
  expect_error(validate_completed_run(run_dir, quiet = TRUE), "status must be 0 or 1")
})

run_test("strict external lp_raw and risk_score aliases must agree", {
  parent <- tempfile("crc9-completed-strict-score-alias-tamper-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE, force = TRUE), add = TRUE)
  run_dir <- build_fixture(parent)
  score_path <- file.path(run_dir, "results", "highest_mean_external_scores.csv")
  scores <- utils::read.csv(score_path, stringsAsFactors = FALSE, check.names = FALSE)
  scores$risk_score[[1L]] <- scores$risk_score[[1L]] + 0.25
  write_csv(scores, score_path)
  reseal_output_manifest(run_dir)
  expect_error(validate_completed_run(run_dir, quiet = TRUE), "lp_raw and risk_score differ")
})

run_test("strict external risk_sd remains cohort-standardized", {
  parent <- tempfile("crc9-completed-strict-scaling-tamper-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE, force = TRUE), add = TRUE)
  run_dir <- build_fixture(parent)
  score_path <- file.path(run_dir, "results", "highest_mean_external_scores.csv")
  scores <- utils::read.csv(score_path, stringsAsFactors = FALSE, check.names = FALSE)
  idx <- scores$cohort == "GSE14333"
  scores$risk_sd[idx] <- scores$risk_sd[idx] * 2
  write_csv(scores, score_path)
  reseal_output_manifest(run_dir)
  expect_error(validate_completed_run(run_dir, quiet = TRUE), "score scaling is invalid")
})

run_test("strict TCGA mapping audit is bound to the certified model", {
  parent <- tempfile("crc9-completed-strict-mapping-identity-tamper-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE, force = TRUE), add = TRUE)
  run_dir <- build_fixture(parent)
  audit_path <- file.path(
    run_dir, "results", "highest_mean_TCGA_COAD_signature_gene_mapping.csv"
  )
  audit <- utils::read.csv(audit_path, stringsAsFactors = FALSE, check.names = FALSE)
  audit$model_fingerprint[[1L]] <- paste(rep("f", 64L), collapse = "")
  write_csv(audit, audit_path)
  reseal_output_manifest(run_dir)
  expect_error(validate_completed_run(run_dir, quiet = TRUE), "wrong model_fingerprint")
})

run_test("strict cohort performance must reconstruct from patient scores", {
  parent <- tempfile("crc9-completed-strict-performance-tamper-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE, force = TRUE), add = TRUE)
  run_dir <- build_fixture(parent)
  performance_path <- file.path(
    run_dir, "results", "highest_mean_external_cohort_performance.csv"
  )
  performance <- utils::read.csv(
    performance_path, stringsAsFactors = FALSE, check.names = FALSE
  )
  performance$HR_per_SD[[1L]] <- performance$HR_per_SD[[1L]] * 1.25
  write_csv(performance, performance_path)
  reseal_output_manifest(run_dir)
  expect_error(validate_completed_run(run_dir, quiet = TRUE), "HR_per_SD.*reconstruct")
})

run_test("strict meta-analysis artifacts retain certified model identity", {
  parent <- tempfile("crc9-completed-strict-meta-identity-tamper-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE, force = TRUE), add = TRUE)
  run_dir <- build_fixture(parent)
  meta_path <- file.path(
    run_dir, "results", "highest_mean_os_validation_meta_analysis.csv"
  )
  meta <- utils::read.csv(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
  meta$analysis_key[[1L]] <- paste(rep("f", 64L), collapse = "")
  write_csv(meta, meta_path)
  reseal_output_manifest(run_dir)
  expect_error(
    validate_completed_run(run_dir, quiet = TRUE),
    "Strict highest-mean OS meta-analysis analysis_key differs"
  )
})

run_test("resealed strict OOF duplicate still fails coverage QA", {
  parent <- tempfile("crc9-completed-oof-tamper-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE, force = TRUE), add = TRUE)
  run_dir <- build_fixture(parent)
  oof_path <- file.path(
    run_dir, "results", "highest_mean_nested_cv_oof_predictions.csv"
  )
  oof <- utils::read.csv(oof_path, stringsAsFactors = FALSE, check.names = FALSE)
  oof$sample[[1L]] <- oof$sample[[2L]]
  write_csv(oof, oof_path)
  reseal_output_manifest(run_dir)
  expect_error(validate_completed_run(run_dir, quiet = TRUE), "OOF.*unique")
})

report <- do.call(rbind, results)
cat("\nSummary:\n")
print(table(report$status))
if (any(report$status != "PASS")) quit(save = "no", status = 1L)
quit(save = "no", status = 0L)
