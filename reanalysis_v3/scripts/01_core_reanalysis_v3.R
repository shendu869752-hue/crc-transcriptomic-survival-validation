locale_result <- try(Sys.setlocale("LC_ALL", "Chinese_China.utf8"), silent = TRUE)
if (inherits(locale_result, "try-error") || is.na(locale_result) || !nzchar(locale_result)) {
  warning("Could not activate Chinese_China.utf8; retaining the current R locale")
}
work_root_boot <- Sys.getenv("CRC_WORK_ROOT", unset = getwd())
source(file.path(work_root_boot, "reanalysis_v3", "scripts", "utils.R"))
source(file.path(
  work_root_boot, "reanalysis_v3", "scripts", "highest_mean_sensitivity.R"
))

main <- function() {
set.seed(2024)
options(stringsAsFactors = FALSE)
run_started_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
PRIMARY_GPL570_METHOD <- "unique_mean"
SENSITIVITY_GPL570_METHOD <-
  "unique_highest_mean_training_partition_frozen"
SENSITIVITY_EXTERNAL_GPL570_METHOD <-
  "full_development_frozen_unique_highest_mean"
gse_expression_path <- file.path(
  PREPARED_INPUT_DIR, "GSE39582_expression_unique_mean.csv.gz"
)
gse_probe_expression_path <- file.path(
  V3_ROOT, "raw_inputs", "GSE39582_series_matrix.txt.gz"
)
gse_clinical_path <- file.path(PREPARED_INPUT_DIR, "GSE39582_clin.csv")
gse_probe_annotation_path <- file.path(
  PREPARED_INPUT_DIR, "GSE39582_probe_annotation.csv"
)
tcga_expression_path <- file.path(PREPARED_INPUT_DIR, "TCGA_COAD_expression.csv.gz")
tcga_survival_path <- file.path(PREPARED_INPUT_DIR, "COAD_survival.txt")
tcga_annotation_path <- file.path(
  PREPARED_INPUT_DIR, "TCGA_COAD_clinical_annotations.csv"
)
geo_paths <- setNames(
  file.path(
    EXTERNAL_GEO_DIR,
    paste0(c("GSE14333", "GSE17536", "GSE17537"), "_series_matrix.txt.gz")
  ),
  c("GSE14333", "GSE17536", "GSE17537")
)
required_inputs <- c(
  GSE39582_expression_unique_mean = gse_expression_path,
  GSE39582_probe_expression = gse_probe_expression_path,
  GSE39582_clinical = gse_clinical_path,
  GSE39582_frozen_probe_annotation = gse_probe_annotation_path,
  TCGA_COAD_expression = tcga_expression_path,
  TCGA_COAD_survival = tcga_survival_path,
  TCGA_COAD_clinical_annotations = tcga_annotation_path,
  GSE14333_series_matrix = unname(geo_paths["GSE14333"]),
  GSE17536_series_matrix = unname(geo_paths["GSE17536"]),
  GSE17537_series_matrix = unname(geo_paths["GSE17537"])
)
code_files <- c(
  core_analysis = file.path(SCRIPT_DIR, "01_core_reanalysis_v3.R"),
  analysis_utils = file.path(SCRIPT_DIR, "utils.R"),
  highest_mean_sensitivity = file.path(
    SCRIPT_DIR, "highest_mean_sensitivity.R"
  ),
  input_acquisition = file.path(SCRIPT_DIR, "00_acquire_prepare_inputs.R"),
  local_launcher = file.path(WORK_ROOT, "run_v3_local.ps1")
)
if (any(!file.exists(c(required_inputs, code_files)))) {
  missing_paths <- c(required_inputs, code_files)[
    !file.exists(c(required_inputs, code_files))
  ]
  stop("Missing run identity files: ", paste(missing_paths, collapse = ", "))
}

PIPELINE_VERSION <- paste0(
  "crc9-strict-nested-v3.5-unique-mean-primary-",
  "highest-mean-training-partition-learned-full-development-frozen-external"
)
requested_repeats <- suppressWarnings(as.integer(Sys.getenv("CRC_NESTED_REPEATS", "10")))
if (!is.finite(requested_repeats) || requested_repeats < 1L) {
  stop("CRC_NESTED_REPEATS must be a positive integer")
}
bootstrap_reps <- suppressWarnings(as.integer(Sys.getenv("CRC_BOOTSTRAP_REPS", "2000")))
if (!is.finite(bootstrap_reps) || bootstrap_reps < 100L) {
  stop("CRC_BOOTSTRAP_REPS must be an integer of at least 100")
}
force_nested <- identical(Sys.getenv("CRC_FORCE_NESTED", "0"), "1")
PIPELINE_PARAMETERS <- list(
  pipeline_version = PIPELINE_VERSION,
  outer_folds = 5L,
  outer_repeats = requested_repeats,
  inner_folds = 5L,
  base_seed = 20240916L,
  sd_cutoff = 0.2,
  univ_p_cutoff = 1e-4,
  min_candidates = 5L,
  candidate_fallback_n = 50L,
  coef_cutoff = 0.05,
  min_selected = 5L,
  selected_fallback_n = 8L,
  max_selected = 30L,
  # Explicit decimal literals are the frozen grid. Besides making the exact
  # search set auditable, these values survive the text dput/dget manifest
  # round-trip bit-for-bit on R 4.6.1.
  lambda_ratio_grid = c(
    0.562341325190349, 0.431933279405154, 0.331767112784286,
    0.254829674797935, 0.195734178148766, 0.150343041978733,
    0.115478198468946, 0.0886985799018192, 0.0681292069057962,
    0.0523299114681495, 0.0401945033361513, 0.0308733199257026,
    0.0237137370566166, 0.0182144753639595, 0.0139905031413729,
    0.0107460782832132, 0.00825404185268019, 0.00633991351172485,
    0.00486967525165863, 0.00374038810036779, 0.00287298483335367,
    0.00220673406908459, 0.00169498815139035, 0.00130191710619008,
    0.001
  ),
  lambda_rule = "one_se_larger_penalty",
  ties = "efron",
  development_probe_aggregation = PRIMARY_GPL570_METHOD
)
validate_pipeline_params(PIPELINE_PARAMETERS)
HIGHEST_MEAN_PARAMETERS <- PIPELINE_PARAMETERS
HIGHEST_MEAN_PARAMETERS$pipeline_version <- paste0(
  PIPELINE_VERSION, "-strict-highest-mean-full-pipeline"
)
HIGHEST_MEAN_PARAMETERS$development_probe_aggregation <-
  "unique_highest_mean_training_partition_learned"
HIGHEST_MEAN_PARAMETERS$external_probe_application <-
  "full_development_frozen"
validate_pipeline_params(HIGHEST_MEAN_PARAMETERS)
run_parameters <- list(
  output_contract_version = "crc9-artifact-contract-v3",
  pipeline_version = PIPELINE_VERSION,
  nested_pipeline = PIPELINE_PARAMETERS,
  highest_mean_nested_pipeline = HIGHEST_MEAN_PARAMETERS,
  bootstrap = list(
    replicates = bootstrap_reps,
    percentile_probabilities = c(0.025, 0.975),
    cohort_seed_base = 20240917L,
    ensemble_seed = 20240930L,
    minimum_valid_fraction = 0.95,
    estimand = "fixed-score patient percentile bootstrap"
  ),
  probe_mapping = list(
    primary = PRIMARY_GPL570_METHOD,
    sensitivity = "unique_highest_mean",
    sensitivity_result_label = SENSITIVITY_GPL570_METHOD,
    sensitivity_external_result_label =
      SENSITIVITY_EXTERNAL_GPL570_METHOD,
    eligibility = "exactly_one_nonempty_SYMBOL",
    primary_multi_probe_rule = "arithmetic_mean_on_log2_scale",
    sensitivity_tie_break = "PROBEID_ascending",
    annotation = "hgu133plus2.db",
    external_sensitivity_scope =
      "training_partition_learned_full_development_frozen_external",
    sensitivity_analysis_key_rule = paste(
      "SHA-256 over raw GSE39582 probe matrix, clinical table, frozen",
      "probe annotation, all code_files, HIGHEST_MEAN_PARAMETERS, R, and",
      "analysis package versions"
    ),
    sensitivity_analysis_key_inputs = c(
      "GSE39582_series_matrix.txt.gz",
      "GSE39582_clin.csv",
      "GSE39582_probe_annotation.csv"
    ),
    sensitivity_analysis_key_code = names(code_files)
  ),
  gene_identity_resolution = list(
    scope = "TCGA_COAD_signature_rows_only",
    rule = "exact_SYMBOL_then_unique_reverse_unique_ENTREZ_alias",
    exact_symbol_priority = TRUE,
    alias_requirements = c(
      "one_expression_row_candidate",
      "one_reverse_current_SYMBOL",
      "one_reverse_ENTREZID",
      "required_SYMBOL_and_ENTREZID_match"
    ),
    annotation = "org.Hs.eg.db",
    annotation_version = as.character(
      utils::packageVersion("org.Hs.eg.db")
    ),
    outcome_blind = TRUE,
    expression_values_used = FALSE,
    development_candidate_space_filtered_by_TCGA = FALSE
  ),
  scoring = list(
    scope = "cohort_adaptive_standardization",
    risk_group_cutpoint = "cohort_median"
  ),
  time_dependent_auc = list(
    years = c(1, 3, 5), min_events = 10L, min_at_risk = 20L,
    weighting = "marginal"
  ),
  survival = list(
    geo_month_to_day = 30.4375,
    year_days = 365.25,
    cox_ties = "efron",
    spline_df = 3L,
    clinical_min_n = 50L,
    clinical_min_events = 20L,
    os_meta_method = "REML",
    os_meta_primary_test = "knha"
  ),
  execution = list(
    force_nested = force_nested,
    entrypoint = "run_v3_local.ps1"
  )
)
run_identity <- make_run_key(
  input_paths = required_inputs,
  code_paths = code_files,
  params = run_parameters,
  return_manifest = TRUE
)
current_run_key <- run_identity$key
current_nested_key <- NA_character_

dir.create(RUNS_ROOT, recursive = TRUE, showWarnings = FALSE)
RUN_TMP_DIR <- file.path(RUNS_ROOT, paste0(current_run_key, ".tmp"))
RUN_FINAL_DIR <- file.path(RUNS_ROOT, current_run_key)
if (dir.exists(RUN_FINAL_DIR) || file.exists(RUN_FINAL_DIR)) {
  stop("Final run directory already exists; refusing to overwrite: ", RUN_FINAL_DIR)
}
if (dir.exists(RUN_TMP_DIR) || file.exists(RUN_TMP_DIR)) {
  stop("Temporary run directory already exists; preserving it for audit: ", RUN_TMP_DIR)
}
if (!dir.create(RUN_TMP_DIR, recursive = FALSE, showWarnings = FALSE)) {
  stop("Could not create temporary run directory: ", RUN_TMP_DIR)
}
RESULT_DIR <- file.path(RUN_TMP_DIR, "results")
FIGURE_DIR <- file.path(RUN_TMP_DIR, "figures")
LOG_DIR <- file.path(RUN_TMP_DIR, "logs")
for (path in c(RESULT_DIR, FIGURE_DIR, LOG_DIR)) {
  if (!dir.create(path, recursive = FALSE, showWarnings = FALSE)) {
    stop("Could not create run subdirectory: ", path)
  }
}

required_result_files <- c(
  "input_manifest.csv", "code_manifest.csv", "run_parameters.txt",
  "primary_cohort_qc.csv", "genomewide_univariate_cox.csv",
  "full_pipeline_audit.csv", "lasso_tuning_curve.csv",
  "full_model_inner_fold_performance.csv",
  "full_model_inner_fold_assignments.csv",
  "training_preprocessing_parameters.csv", "final_model_coefficients.csv",
  "proportional_hazards_test.csv", "final_model.rds", "lasso_path.rds",
  "full_pipeline.rds", "nested_cv_outer_fold_assignments.csv",
  "nested_cv_inner_fold_assignments.csv", "nested_cv_fold_performance.csv",
  "nested_cv_oof_predictions.csv", "nested_cv_repeat_performance.csv",
  "nested_cv_selected_coefficients.csv", "nested_cv_selection_frequency.csv",
  "nested_cv_selection_jaccard.csv", "nested_cv_tuning_performance.csv",
  "nested_cv_model_fingerprints.csv", "nested_cv_failure_log.csv",
  "nested_cv_summary.csv", "nested_cv_cache_metadata.csv",
  "TCGA_COAD_signature_gene_mapping.csv",
  "probe_mapping_audit.csv", "probe_gene_mapping.csv", "all_risk_scores.csv",
  "km_risk_table.csv",
  "cohort_performance.csv", "time_dependent_auc.csv",
  "cindex_bootstrap_intervals.csv", "nested_cv_oof_ensemble.csv",
  "nested_cv_oof_ensemble_cindex_bootstrap.csv",
  "internal_performance_summary.csv", "risk_score_linearity_test.csv",
  "risk_score_proportional_hazards_test.csv", "multivariable_clinical_cox.csv",
  "clinical_adjustment_eligibility.csv",
  "clinical_incremental_performance.csv",
  "clinical_model_proportional_hazards_test.csv",
  "os_validation_meta_inputs.csv", "os_validation_meta_analysis.csv",
  "os_validation_meta_leave_one_out.csv",
  "highest_mean_nested_cv_outer_fold_assignments.csv",
  "highest_mean_nested_cv_inner_fold_assignments.csv",
  "highest_mean_nested_cv_fold_performance.csv",
  "highest_mean_nested_cv_oof_predictions.csv",
  "highest_mean_nested_cv_repeat_performance.csv",
  "highest_mean_nested_cv_selected_coefficients.csv",
  "highest_mean_nested_cv_selection_frequency.csv",
  "highest_mean_nested_cv_selection_jaccard.csv",
  "highest_mean_nested_cv_tuning_performance.csv",
  "highest_mean_nested_cv_inner_tuning_performance.csv",
  "highest_mean_nested_cv_model_fingerprints.csv",
  "highest_mean_nested_cv_failure_log.csv",
  "highest_mean_nested_cv_summary.csv",
  "highest_mean_nested_cv_cache_metadata.csv",
  "highest_mean_outer_probe_maps.csv",
  "highest_mean_inner_probe_map_hashes.csv",
  "highest_mean_outer_probe_frequency.csv",
  "highest_mean_inner_probe_frequency_by_outer.csv",
  "highest_mean_inner_probe_frequency.csv",
  "highest_mean_full_development_probe_map.csv",
  "highest_mean_full_development_model.rds",
  "highest_mean_final_model_coefficients.csv",
  "highest_mean_external_scores.csv",
  "highest_mean_external_cohort_performance.csv",
  "highest_mean_TCGA_COAD_signature_gene_mapping.csv",
  "highest_mean_os_validation_meta_inputs.csv",
  "highest_mean_os_validation_meta_analysis.csv",
  "highest_mean_os_validation_meta_leave_one_out.csv",
  "probe_mapping_sensitivity.csv", "sessionInfo.txt"
)
required_figure_files <- c(
  "Fig1_nested_validation_stability.pdf",
  "Fig1_nested_validation_stability.png",
  "Fig2_cohort_forest.pdf", "Fig2_cohort_forest.png",
  "Fig3_multicohort_KM.pdf", "Fig3_multicohort_KM.png",
  "Fig4_probe_mapping_sensitivity.pdf", "Fig4_probe_mapping_sensitivity.png",
  "Supplementary_FigS1_LASSO_diagnostic.pdf",
  "Supplementary_FigS1_LASSO_diagnostic.png"
)
required_artifacts <- c(
  file.path("results", required_result_files),
  file.path("figures", required_figure_files),
  file.path("code_snapshot", basename(code_files)),
  file.path("logs", "01_core_reanalysis_v3.log")
)
optional_artifacts <- file.path(
  "results",
  c("os_adjusted_meta_inputs.csv", "os_adjusted_meta_analysis.csv")
)
artifact_spec <- data.frame(
  artifact_id = c(required_artifacts, optional_artifacts),
  relative_path = c(required_artifacts, optional_artifacts),
  required = c(
    rep(TRUE, length(required_artifacts)),
    rep(FALSE, length(optional_artifacts))
  ),
  stringsAsFactors = FALSE
)
artifact_registry <- new_artifact_registry(RUN_TMP_DIR, artifact_spec)
run_status_file <- file.path(RESULT_DIR, "run_status.csv")
output_manifest_path <- file.path(RESULT_DIR, "output_manifest.csv")
input_manifest_sha256 <- NA_character_
parameter_manifest_sha256 <- NA_character_
code_manifest_sha256 <- NA_character_
output_manifest_sha256 <- NA_character_
run_complete <- FALSE
log_con <- NULL
output_sink_active <- FALSE
message_sink_active <- FALSE

close_run_log <- function() {
  if (isTRUE(message_sink_active)) {
    sink(type = "message")
    message_sink_active <<- FALSE
  }
  if (isTRUE(output_sink_active)) {
    sink()
    output_sink_active <<- FALSE
  }
  if (!is.null(log_con)) {
    close(log_con)
    log_con <<- NULL
  }
  invisible(TRUE)
}
write_run_status <- function(status, detail = "") {
  write_run_status_file(
    path = run_status_file,
    status = status,
    run_key = current_run_key,
    nested_analysis_key = current_nested_key,
    input_manifest_sha256 = input_manifest_sha256,
    parameter_manifest_sha256 = parameter_manifest_sha256,
    code_manifest_sha256 = code_manifest_sha256,
    output_manifest_sha256 = output_manifest_sha256,
    started_at_utc = run_started_utc,
    detail = detail
  )
}
on.exit({
  close_run_log()
  if (!isTRUE(run_complete) && dir.exists(RUN_TMP_DIR)) {
    try(write_run_status(
      "failed_or_interrupted",
      "The script did not reach the atomic commit gate; inspect the run log"
    ), silent = TRUE)
  }
}, add = TRUE)

input_norm <- normalizePath(unname(required_inputs), winslash = "/")
root_prefix <- paste0(normalizePath(PKG_ROOT, winslash = "/"), "/")
relative_input_paths <- ifelse(
  startsWith(input_norm, root_prefix),
  substring(input_norm, nchar(root_prefix) + 1L),
  basename(input_norm)
)
input_manifest <- data.frame(
  source_id = names(required_inputs),
  relative_path = relative_input_paths,
  bytes = file.info(unname(required_inputs))$size,
  sha256 = unname(run_identity$manifest$inputs[names(required_inputs)]),
  stringsAsFactors = FALSE
)
input_manifest_path <- file.path(RESULT_DIR, "input_manifest.csv")
fwrite(input_manifest, input_manifest_path)
input_manifest_sha256 <- sha256_file(input_manifest_path)

code_snapshot <- snapshot_code_files(code_files, RUN_TMP_DIR)
expected_snapshot_hashes <- unname(
  run_identity$manifest$code[code_snapshot$name]
)
if (!identical(code_snapshot$sha256, expected_snapshot_hashes)) {
  stop("Code changed after run-key construction; refusing an inconsistent snapshot")
}
code_manifest <- rbind(
  data.frame(
    category = "code_sha256",
    name = names(run_identity$manifest$code),
    value = unname(run_identity$manifest$code),
    stringsAsFactors = FALSE
  ),
  data.frame(
    category = "code_snapshot_relative_path",
    name = code_snapshot$name,
    value = code_snapshot$relative_path,
    stringsAsFactors = FALSE
  ),
  data.frame(
    category = "code_snapshot_bytes",
    name = code_snapshot$name,
    value = as.character(code_snapshot$bytes),
    stringsAsFactors = FALSE
  ),
  data.frame(
    category = "runtime",
    name = c("R.version.string", "R.platform"),
    value = c(run_identity$manifest$r_version, run_identity$manifest$r_platform),
    stringsAsFactors = FALSE
  ),
  data.frame(
    category = "package_version",
    name = names(run_identity$manifest$packages),
    value = unname(run_identity$manifest$packages),
    stringsAsFactors = FALSE
  )
)
code_manifest_path <- file.path(RESULT_DIR, "code_manifest.csv")
fwrite(code_manifest, code_manifest_path)
code_manifest_sha256 <- sha256_file(code_manifest_path)

parameter_manifest_path <- file.path(RESULT_DIR, "run_parameters.txt")
write_dput_manifest(run_parameters, parameter_manifest_path)
parameter_manifest_sha256 <- sha256_file(parameter_manifest_path)
write_run_status("running", "Run identity established; outputs are not valid until status=complete")

log_file <- file.path(LOG_DIR, "01_core_reanalysis_v3.log")
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
output_sink_active <- TRUE
sink(log_con, type = "message")
message_sink_active <- TRUE

cat("CRC transcriptomic survival-score A-level reanalysis v3\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Run key:", current_run_key, "\n")
cat("Temporary run directory:", RUN_TMP_DIR, "\n")
cat("Source root:", PKG_ROOT, "\n")
cat("Writable work root:", WORK_ROOT, "\n\n")

cat("Loading gene-level training and TCGA matrices...\n")
gse <- read_expr_csv(gse_expression_path)
tcga <- read_expr_csv(tcga_expression_path)
gse_clin <- fread(gse_clinical_path, data.table = FALSE)
gse_clin$os_event <- suppressWarnings(as.numeric(gse_clin$os_event))
gse_clin$os_delay_months <- suppressWarnings(as.numeric(gse_clin$os_delay_months))
gse_clin$time <- gse_clin$os_delay_months * 30.4375
gse_clin$status <- gse_clin$os_event
gse_clin <- gse_clin[
  gse_clin$sample %in% colnames(gse) & is.finite(gse_clin$time) &
    gse_clin$time > 0 & gse_clin$status %in% 0:1,
  , drop = FALSE
]
gse <- gse[, gse_clin$sample, drop = FALSE]

tcga_surv <- fread(tcga_survival_path, data.table = FALSE)
tcga_surv$OS <- suppressWarnings(as.numeric(tcga_surv$OS))
tcga_surv$OS.time <- suppressWarnings(as.numeric(tcga_surv$OS.time))
tcga_surv <- tcga_surv[
  tcga_surv$sample %in% colnames(tcga) & is.finite(tcga_surv$OS.time) &
    tcga_surv$OS.time > 0 & tcga_surv$OS %in% 0:1,
  c("sample", "OS.time", "OS"), drop = FALSE
]
tcga <- tcga[, tcga_surv$sample, drop = FALSE]

reverse_km_followup <- function(time, status) {
  fit <- survfit(Surv(time, 1 - status) ~ 1)
  tab <- summary(fit)$table
  values <- unname(tab[c("median", "0.95LCL", "0.95UCL")])
  names(values) <- c("median", "lower95", "upper95")
  values
}
gse_followup <- reverse_km_followup(gse_clin$time, gse_clin$status)
tcga_followup <- reverse_km_followup(tcga_surv$OS.time, tcga_surv$OS)

qc <- data.frame(
  cohort = c("GSE39582", "TCGA-COAD"),
  n = c(ncol(gse), ncol(tcga)),
  events = c(sum(gse_clin$status), sum(tcga_surv$OS)),
  genes = c(nrow(gse), nrow(tcga)),
  median_observed_time_days = c(median(gse_clin$time), median(tcga_surv$OS.time)),
  reverse_km_median_followup_days = c(gse_followup["median"], tcga_followup["median"]),
  reverse_km_followup_lower95_days = c(gse_followup["lower95"], tcga_followup["lower95"]),
  reverse_km_followup_upper95_days = c(gse_followup["upper95"], tcga_followup["upper95"]),
  stringsAsFactors = FALSE
)
fwrite(qc, file.path(RESULT_DIR, "primary_cohort_qc.csv"))
print(qc)

# -------------------------------------------------------------------------
# Reproduce the disclosed discovery pipeline on GSE39582.
# -------------------------------------------------------------------------
cat("\nFitting the disclosed development algorithm with training-only preprocessing...\n")

is_noncoding <- function(x) {
  grepl("^LINC|^MIR|^SNORA|^SNORD|^RNU|^RN7S|^SCARNA|^RNY|^RPPH|^MIRLET", x) |
    grepl("-AS[0-9]*$|-IT[0-9]*$", x) | grepl("^AC[0-9]{6}", x)
}

# The feature universe is defined from the development matrix alone. No
# validation-cohort expression values, feature availability, or outcomes are
# used to screen the development features.
feature_universe <- rownames(gse)[!is_noncoding(rownames(gse))]
discovery_expr <- gse[feature_universe, , drop = FALSE]
full_pipeline <- fit_signature_pipeline(
  discovery_expr,
  gse_clin$time,
  gse_clin$status,
  seed = PIPELINE_PARAMETERS$base_seed,
  params = PIPELINE_PARAMETERS
)
final_fit <- full_pipeline$fit
final_coefs <- full_pipeline$coefs
lasso_coef <- full_pipeline$lasso_coefficients
univ_df <- full_pipeline$univ_df
fwrite(univ_df, file.path(RESULT_DIR, "genomewide_univariate_cox.csv"))
fwrite(full_pipeline$audit, file.path(RESULT_DIR, "full_pipeline_audit.csv"))
fwrite(full_pipeline$tuning_curve, file.path(RESULT_DIR, "lasso_tuning_curve.csv"))
fwrite(
  full_pipeline$inner_fold_performance,
  file.path(RESULT_DIR, "full_model_inner_fold_performance.csv")
)
fwrite(
  full_pipeline$inner_fold_assignments,
  file.path(RESULT_DIR, "full_model_inner_fold_assignments.csv")
)
fwrite(data.frame(
  gene = names(final_coefs),
  training_center = unname(full_pipeline$center),
  training_scale = unname(full_pipeline$scale),
  stringsAsFactors = FALSE
), file.path(RESULT_DIR, "training_preprocessing_parameters.csv"))

model_table <- data.frame(
  gene = names(final_coefs),
  coef = unname(final_coefs),
  HR = unname(exp(final_coefs)),
  lower95 = unname(exp(confint(final_fit)[, 1])),
  upper95 = unname(exp(confint(final_fit)[, 2])),
  p = summary(final_fit)$coefficients[, "Pr(>|z|)"],
  lasso_coef = unname(lasso_coef[names(final_coefs)]),
  stringsAsFactors = FALSE
)
fwrite(model_table, file.path(RESULT_DIR, "final_model_coefficients.csv"))
ph <- cox.zph(final_fit)
ph_df <- data.frame(term = rownames(ph$table), ph$table, check.names = FALSE)
fwrite(ph_df, file.path(RESULT_DIR, "proportional_hazards_test.csv"))
saveRDS(final_fit, file.path(RESULT_DIR, "final_model.rds"))
saveRDS(full_pipeline$lasso_path, file.path(RESULT_DIR, "lasso_path.rds"))
saveRDS(full_pipeline, file.path(RESULT_DIR, "full_pipeline.rds"))

cat("Selected genes:", paste(names(final_coefs), collapse = ", "), "\n")
cat("Matches legacy nine-gene set:", setequal(names(final_coefs), SIGNATURE_GENES), "\n")
cat("Apparent training C-index:", round(concordance(final_fit)$concordance, 3), "\n")
print(full_pipeline$audit)

# Base graphics preserve the glmnet path while the second panel reports the
# custom relative-lambda, full-pipeline inner validation used for tuning.
plot_lasso_diagnostic <- function() {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 4.3, 1), cex.main = 1.05)
  plot(
    full_pipeline$lasso_path,
    xvar = "lambda",
    label = FALSE,
    main = "A  LASSO coefficient paths"
  )
  abline(v = log(full_pipeline$lambda), lty = 2, col = "#B2182B")
  tuning_plot <- full_pipeline$tuning_curve[
    order(full_pipeline$tuning_curve$lambda_ratio), , drop = FALSE
  ]
  x <- log10(tuning_plot$lambda_ratio)
  y <- tuning_plot$mean_c_index
  ylim <- range(y - tuning_plot$se_c_index, y + tuning_plot$se_c_index)
  plot(
    x, y, type = "b", pch = 16, ylim = ylim,
    xlab = expression(log[10](lambda/lambda[max])),
    ylab = "Mean inner-fold C-index",
    main = "B  Five-fold full-pipeline tuning"
  )
  arrows(
    x, y - tuning_plot$se_c_index, x, y + tuning_plot$se_c_index,
    angle = 90, code = 3, length = 0.035, col = "grey45"
  )
  abline(v = log10(full_pipeline$lambda_ratio), lty = 2, col = "#B2182B")
  abline(h = unique(tuning_plot$one_se_threshold), lty = 3, col = "grey45")
}
grDevices::cairo_pdf(
  file.path(FIGURE_DIR, "Supplementary_FigS1_LASSO_diagnostic.pdf"),
  width = 7.2, height = 3.8
)
plot_lasso_diagnostic()
dev.off()
png(
  file.path(FIGURE_DIR, "Supplementary_FigS1_LASSO_diagnostic.png"),
  width = 7.2, height = 3.8, units = "in", res = 300,
  type = "cairo-png"
)
plot_lasso_diagnostic()
dev.off()

# -------------------------------------------------------------------------
# Strict repeated nested cross-validation: every data-dependent discovery step
# occurs within each outer training fold. Test predictions are put on a common
# training-linear-predictor SD scale before repeat-level OOF concordance.
# -------------------------------------------------------------------------
training_input_files <- c(
  gse_expression_path,
  gse_clinical_path,
  gse_probe_annotation_path
)
analysis_identity <- make_analysis_key(
  input_paths = training_input_files,
  code_paths = code_files,
  params = PIPELINE_PARAMETERS,
  return_manifest = TRUE
)
analysis_key <- analysis_identity$key
current_nested_key <- analysis_key
write_run_status(
  "running",
  paste0("run_key=", current_run_key, "; nested_analysis_key=", current_nested_key)
)
cache_dir <- file.path(V3_ROOT, "cache", "nested_cv")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
cache_file <- file.path(cache_dir, paste0("nested_cv_", analysis_key, ".rds"))

if (file.exists(cache_file) && !force_nested) {
  cat("\nLoading strict nested-CV cache:", basename(cache_file), "\n")
  nested_result <- readRDS(cache_file)
  validate_nested_result(
    nested_result,
    discovery_expr,
    gse_clin$time,
    gse_clin$status,
    gse_clin$sample,
    PIPELINE_PARAMETERS,
    analysis_key
  )
} else {
  cat(
    "\nStrict nested CV:", PIPELINE_PARAMETERS$outer_repeats, "repeats x",
    PIPELINE_PARAMETERS$outer_folds, "outer folds x",
    PIPELINE_PARAMETERS$inner_folds, "inner folds\n"
  )
  nested_result <- tryCatch(
    run_nested_cv(
      discovery_expr,
      gse_clin$time,
      gse_clin$status,
      gse_clin$sample,
      PIPELINE_PARAMETERS,
      analysis_key
    ),
    error = function(e) {
      failure <- data.frame(
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
        analysis_key = analysis_key,
        error = conditionMessage(e),
        stringsAsFactors = FALSE
      )
      fwrite(failure, file.path(LOG_DIR, "nested_cv_failure_log.csv"))
      stop(conditionMessage(e), call. = FALSE)
    }
  )
  nested_result$hash_manifest <- analysis_identity$manifest
  if (!identical(nested_result$pipeline_version, PIPELINE_VERSION)) {
    stop("Fresh nested-CV result pipeline version mismatch")
  }
  validate_nested_result(
    nested_result,
    discovery_expr,
    gse_clin$time,
    gse_clin$status,
    gse_clin$sample,
    PIPELINE_PARAMETERS,
    analysis_key
  )
  saveRDS(nested_result, cache_file, version = 3)
}
if (is.null(nested_result$hash_manifest) ||
    !identical(nested_result$hash_manifest, analysis_identity$manifest)) {
  stop("Nested-CV cache hash manifest does not match the current analysis identity")
}
if (!identical(nested_result$pipeline_version, PIPELINE_VERSION)) {
  stop("Nested-CV cache pipeline version mismatch")
}

nested_df <- nested_result$fold_performance
nested_oof <- nested_result$oof
nested_repeat <- nested_result$repeat_performance
selection_frequency <- nested_result$selection_frequency
nested_summary <- nested_result$summary
fwrite(nested_result$outer_fold_assignments, file.path(RESULT_DIR, "nested_cv_outer_fold_assignments.csv"))
fwrite(nested_result$inner_fold_assignments, file.path(RESULT_DIR, "nested_cv_inner_fold_assignments.csv"))
fwrite(nested_df, file.path(RESULT_DIR, "nested_cv_fold_performance.csv"))
fwrite(nested_oof, file.path(RESULT_DIR, "nested_cv_oof_predictions.csv"))
fwrite(nested_repeat, file.path(RESULT_DIR, "nested_cv_repeat_performance.csv"))
fwrite(nested_result$selection, file.path(RESULT_DIR, "nested_cv_selected_coefficients.csv"))
fwrite(selection_frequency, file.path(RESULT_DIR, "nested_cv_selection_frequency.csv"))
fwrite(nested_result$jaccard, file.path(RESULT_DIR, "nested_cv_selection_jaccard.csv"))
fwrite(nested_result$tuning_performance, file.path(RESULT_DIR, "nested_cv_tuning_performance.csv"))
fwrite(nested_result$model_fingerprints, file.path(RESULT_DIR, "nested_cv_model_fingerprints.csv"))
fwrite(nested_result$failure_log, file.path(RESULT_DIR, "nested_cv_failure_log.csv"))
fwrite(nested_summary, file.path(RESULT_DIR, "nested_cv_summary.csv"))

cache_metadata <- rbind(
  data.frame(
    category = "analysis", name = c("analysis_key", "pipeline_version"),
    value = c(analysis_key, PIPELINE_VERSION), stringsAsFactors = FALSE
  ),
  data.frame(
    category = "input_sha256", name = names(analysis_identity$manifest$inputs),
    value = unname(analysis_identity$manifest$inputs), stringsAsFactors = FALSE
  ),
  data.frame(
    category = "code_sha256", name = names(analysis_identity$manifest$code),
    value = unname(analysis_identity$manifest$code), stringsAsFactors = FALSE
  ),
  data.frame(
    category = "package_version", name = names(analysis_identity$manifest$packages),
    value = unname(analysis_identity$manifest$packages), stringsAsFactors = FALSE
  ),
  data.frame(
    category = "runtime", name = "R.version.string",
    value = analysis_identity$manifest$r_version, stringsAsFactors = FALSE
  ),
  data.frame(
    category = "parameters", name = "dput",
    value = paste(capture.output(dput(PIPELINE_PARAMETERS)), collapse = " "),
    stringsAsFactors = FALSE
  )
)
fwrite(cache_metadata, file.path(RESULT_DIR, "nested_cv_cache_metadata.csv"))
print(nested_summary)

# Main Fig1 is generated only after the primary repeated nested-CV object has
# passed validation. It separates apparent fit, patient-level OOF performance,
# outer-model feature-selection frequency, and within-repeat model overlap.
flow_nodes <- data.frame(
  x = c(0.5, 1.5, 1.5, 0.5, 0.5, 1.5),
  y = c(3, 3, 2, 2, 1, 1),
  label = c(
    sprintf("GSE39582 development\nn=%d; OS events=%d", nrow(gse_clin), sum(gse_clin$status)),
    sprintf(
      "Repeated nested CV\n%d repeats x %d outer\n%d inner folds",
      PIPELINE_PARAMETERS$outer_repeats,
      PIPELINE_PARAMETERS$outer_folds,
      PIPELINE_PARAMETERS$inner_folds
    ),
    "Patient-level OOF\nC-index + stability",
    "Full-development refit\nsame training function",
    "Four external cohorts\ncontinuous score; OS/DFS",
    "REML-KH synthesis\nthree external OS cohorts"
  ),
  node_type = c(
    "development", "internal", "internal",
    "development", "external", "external"
  ),
  stringsAsFactors = FALSE
)
flow_edges <- data.frame(
  x = c(0.82, 1.5, 1.18, 0.5, 0.82),
  y = c(3, 2.72, 2, 1.72, 1),
  xend = c(1.18, 1.5, 0.82, 0.5, 1.18),
  yend = c(3, 2.28, 2, 1.28, 1)
)
p_flow <- ggplot() +
  geom_segment(
    data = flow_edges,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.55, color = "#555555",
    arrow = grid::arrow(length = grid::unit(0.055, "inches"), type = "closed")
  ) +
  geom_label(
    data = flow_nodes,
    aes(x = x, y = y, label = label, fill = node_type),
    family = "Arial", size = 2.35, lineheight = 0.92,
    linewidth = 0.25, label.padding = grid::unit(0.10, "lines")
  ) +
  scale_fill_manual(values = c(
    development = "#DCEAF5", internal = "#EEF3F7", external = "#F7E6E4"
  )) +
  coord_cartesian(xlim = c(0.05, 1.95), ylim = c(0.55, 3.45), clip = "off") +
  labs(title = "Analysis design", tag = "A") +
  theme_void(base_size = 8, base_family = "Arial") +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 8, hjust = 0),
    plot.tag = element_text(face = "bold", size = 9)
  )

apparent_c_index <- unname(concordance(final_fit)$concordance)
n_oof_repeats <- nrow(nested_repeat)
oof_offsets <- if (n_oof_repeats == 1L) {
  0
} else {
  seq(-0.18, 0.18, length.out = n_oof_repeats)
}
performance_points <- rbind(
  data.frame(
    x = 1, c_index = apparent_c_index,
    estimate = "Apparent full-data", stringsAsFactors = FALSE
  ),
  data.frame(
    x = 2 + oof_offsets, c_index = nested_repeat$c_index,
    estimate = "Patient-level OOF by repeat", stringsAsFactors = FALSE
  )
)
performance_summary <- data.frame(
  x = c(1, 2),
  c_index = c(apparent_c_index, median(nested_repeat$c_index)),
  stringsAsFactors = FALSE
)
performance_range <- range(c(0.5, performance_points$c_index))
performance_pad <- max(0.03, diff(performance_range) * 0.12)
performance_limits <- c(
  max(0, performance_range[1] - performance_pad),
  min(1, performance_range[2] + performance_pad)
)
p_performance <- ggplot(
  performance_points,
  aes(x = x, y = c_index, fill = estimate, shape = estimate)
) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "grey65", linewidth = 0.45) +
  geom_point(size = 2.5, color = "#222222", stroke = 0.4) +
  geom_segment(
    data = performance_summary,
    aes(x = x - 0.22, xend = x + 0.22, y = c_index, yend = c_index),
    inherit.aes = FALSE, linewidth = 0.8, color = "#222222"
  ) +
  scale_fill_manual(values = c(
    "Apparent full-data" = "#999999",
    "Patient-level OOF by repeat" = "#2166AC"
  )) +
  scale_shape_manual(values = c(
    "Apparent full-data" = 23,
    "Patient-level OOF by repeat" = 21
  )) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Apparent\nfull-data", "Patient-level OOF\nby repeat"),
    limits = c(0.62, 2.38)
  ) +
  scale_y_continuous(limits = performance_limits) +
  labs(
    title = "Apparent versus internal-validation discrimination",
    subtitle = "Dots are repeat-level estimates; bar denotes the median",
    x = NULL, y = "Harrell C-index", tag = "B"
  ) +
  theme_pub(8) +
  theme(
    text = element_text(family = "Arial"), legend.position = "none",
    plot.title = element_text(face = "bold", size = 8),
    plot.subtitle = element_text(size = 6.5),
    plot.tag = element_text(face = "bold", size = 9)
  )

if (!nrow(selection_frequency)) {
  stop("Primary nested-CV selection-frequency table is empty")
}
selection_top_n <- min(12L, nrow(selection_frequency))
selection_plot_data <- selection_frequency[
  order(-selection_frequency$frequency, -selection_frequency$sign_consistency,
        selection_frequency$gene),
  , drop = FALSE
][seq_len(selection_top_n), , drop = FALSE]
selection_plot_data$gene <- factor(
  selection_plot_data$gene,
  levels = rev(selection_plot_data$gene)
)
outer_model_n <- nrow(nested_df)
selection_plot_data$count_label <- sprintf(
  "%d/%d", selection_plot_data$selected_folds, outer_model_n
)
p_selection <- ggplot(selection_plot_data, aes(frequency, gene)) +
  geom_segment(
    aes(x = 0, xend = frequency, yend = gene),
    linewidth = 0.55, color = "#9EBCD4"
  ) +
  geom_point(shape = 21, size = 2.2, fill = "#2166AC", color = "#222222", stroke = 0.3) +
  geom_text(
    aes(label = count_label), nudge_x = 0.025, hjust = 0,
    size = 2.15, family = "Arial"
  ) +
  scale_x_continuous(
    breaks = seq(0, 1, 0.25),
    labels = paste0(seq(0, 100, 25), "%"),
    limits = c(0, 1.14), expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Outer-model gene-selection frequency",
    subtitle = sprintf(
      "Top %d of %d genes selected across %d outer models",
      selection_top_n, nrow(selection_frequency), outer_model_n
    ),
    x = "Selection frequency", y = NULL, tag = "C"
  ) +
  theme_pub(8) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(face = "bold", size = 8),
    plot.subtitle = element_text(size = 6.5),
    plot.tag = element_text(face = "bold", size = 9)
  )

jaccard_source <- nested_result$jaccard
repeat_levels <- sort(unique(jaccard_source$repeat_id))
jaccard_plot_data <- do.call(rbind, lapply(seq_along(repeat_levels), function(i) {
  z <- jaccard_source[jaccard_source$repeat_id == repeat_levels[[i]], , drop = FALSE]
  offsets <- if (nrow(z) == 1L) 0 else seq(-0.17, 0.17, length.out = nrow(z))
  z$x <- i
  z$x_jitter <- i + offsets
  z
}))
p_jaccard <- ggplot(jaccard_plot_data, aes(x = x, y = jaccard)) +
  geom_boxplot(
    aes(group = x), width = 0.48, outlier.shape = NA,
    fill = "#DCEAF5", color = "#555555", linewidth = 0.45
  ) +
  geom_point(
    aes(x = x_jitter), shape = 21, size = 1.45,
    fill = "#2166AC", color = "#222222", stroke = 0.25, alpha = 0.85
  ) +
  scale_x_continuous(breaks = seq_along(repeat_levels), labels = repeat_levels) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    title = "Within-repeat feature-set overlap",
    subtitle = "Pairwise Jaccard across five outer-fold models",
    x = "Nested-CV repeat", y = "Jaccard index", tag = "D"
  ) +
  theme_pub(8) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(face = "bold", size = 8),
    plot.subtitle = element_text(size = 6.5),
    plot.tag = element_text(face = "bold", size = 9)
  )

draw_nested_validation_figure <- function() {
  grid::grid.newpage()
  layout <- grid::grid.layout(nrow = 2L, ncol = 2L)
  grid::pushViewport(grid::viewport(layout = layout))
  figure_panels <- list(p_flow, p_performance, p_selection, p_jaccard)
  for (i in seq_along(figure_panels)) {
    print(
      figure_panels[[i]],
      vp = grid::viewport(
        layout.pos.row = ((i - 1L) %/% 2L) + 1L,
        layout.pos.col = ((i - 1L) %% 2L) + 1L
      )
    )
  }
  grid::popViewport()
}
grDevices::cairo_pdf(
  file.path(FIGURE_DIR, "Fig1_nested_validation_stability.pdf"),
  width = 7.2, height = 6.6
)
draw_nested_validation_figure()
dev.off()
png(
  file.path(FIGURE_DIR, "Fig1_nested_validation_stability.png"),
  width = 7.2, height = 6.6, units = "in", res = 300,
  type = "cairo-png"
)
draw_nested_validation_figure()
dev.off()

# -------------------------------------------------------------------------
# Pre-specified full-pipeline GPL570 highest-mean sensitivity. Probe choice is
# learned inside every inner/outer training partition. The full-development
# map is then frozen unchanged for all GPL570 external cohorts.
# -------------------------------------------------------------------------
cat("\nLoading probe-level GSE39582 data for strict highest-mean sensitivity...\n")
gse_probe <- read_geo_matrix(gse_probe_expression_path)
missing_development_samples <- setdiff(gse_clin$sample, colnames(gse_probe))
if (length(missing_development_samples)) {
  stop(
    "Probe-level GSE39582 matrix is missing eligible development samples: ",
    paste(missing_development_samples, collapse = ", ")
  )
}
gse_probe <- gse_probe[, gse_clin$sample, drop = FALSE]
if (!identical(colnames(gse_probe), as.character(gse_clin$sample))) {
  stop("Probe-level GSE39582 and clinical sample orders differ")
}
highest_probe_annotation <- fread(
  gse_probe_annotation_path, data.table = FALSE
)
if (!all(c("PROBEID", "SYMBOL") %in% names(highest_probe_annotation))) {
  stop("Frozen GSE39582 probe annotation lacks PROBEID or SYMBOL")
}
highest_probe_annotation$PROBEID <- as.character(
  highest_probe_annotation$PROBEID
)
highest_probe_annotation$SYMBOL <- as.character(
  highest_probe_annotation$SYMBOL
)
coding_symbol <- !is.na(highest_probe_annotation$SYMBOL) &
  nzchar(highest_probe_annotation$SYMBOL) &
  !is_noncoding(highest_probe_annotation$SYMBOL)
highest_probe_annotation <- highest_probe_annotation[
  coding_symbol, c("PROBEID", "SYMBOL"), drop = FALSE
]
if (!nrow(highest_probe_annotation)) {
  stop("No coding one-to-one probe annotations remain for sensitivity analysis")
}

highest_training_input_files <- c(
  gse_probe_expression_path,
  gse_clinical_path,
  gse_probe_annotation_path
)
highest_analysis_identity <- make_analysis_key(
  input_paths = highest_training_input_files,
  code_paths = code_files,
  params = HIGHEST_MEAN_PARAMETERS,
  return_manifest = TRUE
)
highest_mean_analysis_key <- highest_analysis_identity$key
write_run_status(
  "running",
  paste0(
    "run_key=", current_run_key,
    "; primary_nested_analysis_key=", analysis_key,
    "; highest_mean_nested_analysis_key=", highest_mean_analysis_key
  )
)
highest_cache_file <- file.path(
  cache_dir,
  paste0("highest_mean_nested_cv_", highest_mean_analysis_key, ".rds")
)

if (file.exists(highest_cache_file) && !force_nested) {
  cat("Loading strict highest-mean nested-CV cache:",
      basename(highest_cache_file), "\n")
  highest_nested_result <- readRDS(highest_cache_file)
  validate_probe_nested_result(
    highest_nested_result,
    gse_probe,
    gse_clin$time,
    gse_clin$status,
    gse_clin$sample,
    HIGHEST_MEAN_PARAMETERS,
    highest_mean_analysis_key,
    highest_probe_annotation
  )
} else {
  cat(
    "Strict highest-mean nested CV:",
    HIGHEST_MEAN_PARAMETERS$outer_repeats, "repeats x",
    HIGHEST_MEAN_PARAMETERS$outer_folds, "outer folds x",
    HIGHEST_MEAN_PARAMETERS$inner_folds, "inner folds\n"
  )
  highest_nested_result <- tryCatch(
    run_probe_nested_cv(
      gse_probe,
      gse_clin$time,
      gse_clin$status,
      gse_clin$sample,
      HIGHEST_MEAN_PARAMETERS,
      highest_mean_analysis_key,
      highest_probe_annotation
    ),
    error = function(e) {
      failure <- data.frame(
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
        analysis_key = highest_mean_analysis_key,
        error = conditionMessage(e),
        stringsAsFactors = FALSE
      )
      fwrite(
        failure,
        file.path(LOG_DIR, "highest_mean_nested_cv_failure_log.csv")
      )
      stop(conditionMessage(e), call. = FALSE)
    }
  )
  highest_nested_result$hash_manifest <- highest_analysis_identity$manifest
  if (!identical(
      highest_nested_result$pipeline_version,
      HIGHEST_MEAN_PARAMETERS$pipeline_version)) {
    stop("Fresh highest-mean nested-CV pipeline version mismatch")
  }
  validate_probe_nested_result(
    highest_nested_result,
    gse_probe,
    gse_clin$time,
    gse_clin$status,
    gse_clin$sample,
    HIGHEST_MEAN_PARAMETERS,
    highest_mean_analysis_key,
    highest_probe_annotation
  )
  saveRDS(highest_nested_result, highest_cache_file, version = 3)
}
if (is.null(highest_nested_result$hash_manifest) ||
    !identical(
      highest_nested_result$hash_manifest,
      highest_analysis_identity$manifest
    )) {
  stop(
    "Highest-mean nested-CV cache hash manifest does not match current identity"
  )
}
if (!identical(
    highest_nested_result$pipeline_version,
    HIGHEST_MEAN_PARAMETERS$pipeline_version)) {
  stop("Highest-mean nested-CV cache pipeline version mismatch")
}

highest_nested_files <- list(
  highest_mean_nested_cv_outer_fold_assignments.csv =
    highest_nested_result$outer_fold_assignments,
  highest_mean_nested_cv_inner_fold_assignments.csv =
    highest_nested_result$inner_fold_assignments,
  highest_mean_nested_cv_fold_performance.csv =
    highest_nested_result$fold_performance,
  highest_mean_nested_cv_oof_predictions.csv = highest_nested_result$oof,
  highest_mean_nested_cv_repeat_performance.csv =
    highest_nested_result$repeat_performance,
  highest_mean_nested_cv_selected_coefficients.csv =
    highest_nested_result$selection,
  highest_mean_nested_cv_selection_frequency.csv =
    highest_nested_result$selection_frequency,
  highest_mean_nested_cv_selection_jaccard.csv = highest_nested_result$jaccard,
  highest_mean_nested_cv_tuning_performance.csv =
    highest_nested_result$tuning_performance,
  highest_mean_nested_cv_inner_tuning_performance.csv =
    highest_nested_result$inner_tuning_performance,
  highest_mean_nested_cv_model_fingerprints.csv =
    highest_nested_result$model_fingerprints,
  highest_mean_nested_cv_failure_log.csv = highest_nested_result$failure_log,
  highest_mean_nested_cv_summary.csv = highest_nested_result$summary,
  highest_mean_outer_probe_maps.csv = highest_nested_result$outer_probe_maps,
  highest_mean_inner_probe_map_hashes.csv =
    highest_nested_result$inner_probe_map_hashes,
  highest_mean_outer_probe_frequency.csv =
    highest_nested_result$outer_probe_frequency,
  highest_mean_inner_probe_frequency_by_outer.csv =
    highest_nested_result$inner_probe_frequency_by_outer,
  highest_mean_inner_probe_frequency.csv =
    highest_nested_result$inner_probe_frequency
)
for (file_name in names(highest_nested_files)) {
  fwrite(
    highest_nested_files[[file_name]],
    file.path(RESULT_DIR, file_name)
  )
}

highest_cache_metadata <- rbind(
  data.frame(
    category = "analysis",
    name = c("analysis_key", "pipeline_version"),
    value = c(
      highest_mean_analysis_key,
      HIGHEST_MEAN_PARAMETERS$pipeline_version
    ),
    stringsAsFactors = FALSE
  ),
  data.frame(
    category = "input_sha256",
    name = names(highest_analysis_identity$manifest$inputs),
    value = unname(highest_analysis_identity$manifest$inputs),
    stringsAsFactors = FALSE
  ),
  data.frame(
    category = "code_sha256",
    name = names(highest_analysis_identity$manifest$code),
    value = unname(highest_analysis_identity$manifest$code),
    stringsAsFactors = FALSE
  ),
  data.frame(
    category = "package_version",
    name = names(highest_analysis_identity$manifest$packages),
    value = unname(highest_analysis_identity$manifest$packages),
    stringsAsFactors = FALSE
  ),
  data.frame(
    category = "runtime",
    name = "R.version.string",
    value = highest_analysis_identity$manifest$r_version,
    stringsAsFactors = FALSE
  ),
  data.frame(
    category = "parameters",
    name = "dput",
    value = paste(
      capture.output(dput(HIGHEST_MEAN_PARAMETERS)), collapse = " "
    ),
    stringsAsFactors = FALSE
  )
)
fwrite(
  highest_cache_metadata,
  file.path(RESULT_DIR, "highest_mean_nested_cv_cache_metadata.csv")
)

highest_full_model <- highest_nested_result$full_development_model
highest_full_map <- highest_nested_result$full_development_probe_map
highest_full_map$map_hash <- highest_nested_result$full_development_probe_map_hash
highest_full_map$training_sample_hash <-
  highest_nested_result$full_development_training_sample_hash
highest_full_map$training_n <- highest_full_model$probe_map$training_n
highest_full_map$analysis_key <- highest_mean_analysis_key
fwrite(
  highest_full_map,
  file.path(RESULT_DIR, "highest_mean_full_development_probe_map.csv")
)
saveRDS(
  highest_full_model,
  file.path(RESULT_DIR, "highest_mean_full_development_model.rds"),
  version = 3
)
highest_final_coefs <- highest_full_model$coefs[highest_full_model$genes]
highest_coef_table <- data.frame(
  gene = names(highest_final_coefs),
  coef = unname(highest_final_coefs),
  probe_map_hash = highest_full_model$probe_map_hash,
  model_fingerprint = highest_full_model$model_fingerprint,
  analysis_key = highest_mean_analysis_key,
  stringsAsFactors = FALSE
)
fwrite(
  highest_coef_table,
  file.path(RESULT_DIR, "highest_mean_final_model_coefficients.csv")
)
print(highest_nested_result$summary)
rm(gse_probe)
invisible(gc())

# -------------------------------------------------------------------------
# Apply fixed coefficients with cohort-adaptive preprocessing. These analyses
# estimate within-cohort standardized prognostic associations, not a fully
# locked patient-level absolute-risk model.
# -------------------------------------------------------------------------
train_risk <- full_pipeline$train_lp
train_d <- data.frame(
  sample = gse_clin$sample, time = gse_clin$time, status = gse_clin$status,
  stage = gse_clin$tnm_stage, age = suppressWarnings(as.numeric(gse_clin$age)),
  endpoint = "OS", cohort = "GSE39582 (training)",
  mapping_method = PRIMARY_GPL570_METHOD,
  preprocessing_scope = "development_unique_probe_mean_then_training_standardization",
  estimand = "apparent_development_association",
  risk = train_risk, risk_sd = as.numeric(scale(train_risk)),
  group = factor(ifelse(train_risk > median(train_risk), "High", "Low"), levels = c("Low", "High"))
)

tcga_clin <- data.frame(
  sample = tcga_surv$sample,
  time = tcga_surv$OS.time,
  status = tcga_surv$OS,
  stage = NA_character_, age = NA_real_, endpoint = "OS",
  stringsAsFactors = FALSE
)
tcga_gene_resolution <- resolve_expression_gene_rows(
  tcga,
  required_genes = names(final_coefs),
  cohort = "TCGA-COAD"
)
fwrite(
  tcga_gene_resolution$audit,
  file.path(RESULT_DIR, "TCGA_COAD_signature_gene_mapping.csv")
)
tcga_d <- apply_signature(
  tcga_gene_resolution$expr,
  tcga_clin,
  final_coefs,
  "TCGA-COAD",
  "gene_level_input"
)
tcga_annot <- fread(tcga_annotation_path, data.table = FALSE)
tcga_d$stage <- tcga_annot$pathologic_stage[match(tcga_d$sample, tcga_annot$sampleID)]
tcga_d$age <- suppressWarnings(as.numeric(tcga_annot$age_at_initial_pathologic_diagnosis[match(tcga_d$sample, tcga_annot$sampleID)]))

# The sensitivity signature uses its own selected genes. TCGA is gene-level,
# so probe identities are not applicable; resolve SYMBOL/alias rows under the
# same fail-closed rules and retain a separate audit bound to the sensitivity
# model and its originating GPL570 map.
highest_tcga_gene_resolution <- resolve_expression_gene_rows(
  tcga,
  required_genes = names(highest_final_coefs),
  cohort = "TCGA-COAD strict highest-mean sensitivity"
)
highest_tcga_mapping_audit <- highest_tcga_gene_resolution$audit
highest_tcga_mapping_audit$probe_map_hash <- highest_full_model$probe_map_hash
highest_tcga_mapping_audit$model_fingerprint <-
  highest_full_model$model_fingerprint
highest_tcga_mapping_audit$analysis_key <- highest_mean_analysis_key
fwrite(
  highest_tcga_mapping_audit,
  file.path(
    RESULT_DIR, "highest_mean_TCGA_COAD_signature_gene_mapping.csv"
  )
)
highest_tcga_d <- apply_signature(
  highest_tcga_gene_resolution$expr,
  tcga_clin,
  highest_final_coefs,
  "TCGA-COAD",
  "gene_level_input"
)
highest_tcga_d$stage <- tcga_annot$pathologic_stage[
  match(highest_tcga_d$sample, tcga_annot$sampleID)
]
highest_tcga_d$age <- suppressWarnings(as.numeric(
  tcga_annot$age_at_initial_pathologic_diagnosis[
    match(highest_tcga_d$sample, tcga_annot$sampleID)
  ]
))
highest_tcga_d$source_mapping <-
  "TCGA_gene_level_exact_SYMBOL_or_unique_reverse_ENTREZ_alias"

geo_names <- c("GSE14333", "GSE17536", "GSE17537")
geo_data <- list()
highest_geo_data <- list()
probe_mappings <- list()
probe_gene_mappings <- list()
for (nm in geo_names) {
  cat("\nParsing", nm, "with per-sample clinical fields and one-to-one probe mapping...\n")
  path <- geo_paths[[nm]]
  probe_mat <- read_geo_matrix(path)
  meta <- read_geo_metadata(path)
  clin <- if (nm == "GSE14333") parse_gse14333_clin(meta) else parse_gse17x_clin(meta)
  cat("Clinical usable:", sum(is.finite(clin$time) & clin$time > 0 & clin$status %in% 0:1), "/", nrow(clin), "\n")
  primary_mapped <- map_probes(
    probe_mat,
    method = PRIMARY_GPL570_METHOD,
    required_genes = names(final_coefs)
  )
  probe_mappings[[paste(nm, PRIMARY_GPL570_METHOD, sep = "_")]] <-
    transform(primary_mapped$probe_audit, cohort = nm)
  probe_gene_mappings[[paste(nm, PRIMARY_GPL570_METHOD, sep = "_")]] <-
    transform(primary_mapped$mapping, cohort = nm)
  geo_data[[paste(nm, PRIMARY_GPL570_METHOD, sep = "_")]] <-
    apply_signature(
      primary_mapped$expr,
      clin,
      final_coefs,
      nm,
      PRIMARY_GPL570_METHOD
    )

  frozen_highest_expr <- apply_frozen_probe_map(
    probe_mat, highest_full_model$probe_map
  )
  frozen_hash <- attr(frozen_highest_expr, "probe_map_hash")
  if (!identical(frozen_hash, highest_full_model$probe_map_hash)) {
    stop(nm, " frozen highest-mean probe-map hash mismatch")
  }
  highest_geo_data[[nm]] <- apply_signature(
    frozen_highest_expr,
    clin,
    highest_final_coefs,
    nm,
    SENSITIVITY_EXTERNAL_GPL570_METHOD
  )
  highest_geo_data[[nm]]$source_mapping <-
    "GPL570_full_development_frozen_PROBEID"
  frozen_mapping_audit <- highest_full_model$probe_map$mapping
  frozen_mapping_audit$method <- SENSITIVITY_EXTERNAL_GPL570_METHOD
  frozen_mapping_audit$cohort <- nm
  frozen_mapping_audit$map_hash <- highest_full_model$probe_map_hash
  probe_gene_mappings[[paste(nm, SENSITIVITY_GPL570_METHOD, sep = "_")]] <-
    frozen_mapping_audit
  rm(probe_mat)
  invisible(gc())
}
fwrite(rbindlist(probe_mappings, fill = TRUE), file.path(RESULT_DIR, "probe_mapping_audit.csv"))
fwrite(rbindlist(probe_gene_mappings, fill = TRUE), file.path(RESULT_DIR, "probe_gene_mapping.csv"))

highest_external_d <- rbindlist(
  c(list(highest_tcga_d), highest_geo_data), fill = TRUE
)
if (!setequal(
    unique(highest_external_d$cohort),
    c("TCGA-COAD", "GSE14333", "GSE17536", "GSE17537"))) {
  stop("Strict highest-mean external scores do not cover all four cohorts")
}
highest_external_d$probe_map_hash <- highest_full_model$probe_map_hash
highest_external_d$model_fingerprint <-
  highest_full_model$model_fingerprint
highest_external_d$analysis_key <- highest_mean_analysis_key
# Stable schema aliases used by the completed-run validator and downstream
# reporting. `risk` is the raw cohort-adaptive linear score returned by
# apply_signature(); no training-cohort centering or scaling is reintroduced.
highest_external_d$lp_raw <- highest_external_d$risk
highest_external_d$risk_score <- highest_external_d$risk
if (any(!is.finite(highest_external_d$risk)) ||
    any(!is.finite(highest_external_d$risk_sd)) ||
    anyNA(highest_external_d$probe_map_hash) ||
    anyNA(highest_external_d$model_fingerprint) ||
    anyNA(highest_external_d$analysis_key)) {
  stop("Strict highest-mean external scores contain invalid values or identity")
}
highest_score_scaling <- highest_external_d[, .(
  risk_mean = mean(risk_sd),
  risk_sd_observed = sd(risk_sd),
  raw_risk_sd = sd(risk),
  probe_map_hashes = uniqueN(probe_map_hash),
  model_fingerprints = uniqueN(model_fingerprint),
  analysis_keys = uniqueN(analysis_key)
), by = cohort]
if (any(abs(highest_score_scaling$risk_mean) > 1e-10) ||
    any(abs(highest_score_scaling$risk_sd_observed - 1) > 1e-10) ||
    any(!is.finite(highest_score_scaling$raw_risk_sd)) ||
    any(highest_score_scaling$raw_risk_sd <= 0) ||
    any(highest_score_scaling$probe_map_hashes != 1L) ||
    any(highest_score_scaling$model_fingerprints != 1L) ||
    any(highest_score_scaling$analysis_keys != 1L)) {
  stop(
    "Strict highest-mean scores violate cohort SD scaling or identity binding"
  )
}
fwrite(
  highest_external_d,
  file.path(RESULT_DIR, "highest_mean_external_scores.csv")
)
highest_stats <- rbindlist(lapply(
  split(highest_external_d, highest_external_d$cohort),
  cohort_stats
), fill = TRUE)
highest_stats$probe_map_hash <- highest_full_model$probe_map_hash
highest_stats$model_fingerprint <- highest_full_model$model_fingerprint
highest_stats$analysis_key <- highest_mean_analysis_key
fwrite(
  highest_stats,
  file.path(RESULT_DIR, "highest_mean_external_cohort_performance.csv")
)

all_d <- rbindlist(c(list(train_d, tcga_d), geo_data), fill = TRUE)
fwrite(all_d, file.path(RESULT_DIR, "all_risk_scores.csv"))

stats_list <- lapply(split(all_d, interaction(all_d$cohort, all_d$mapping_method, drop = TRUE)), cohort_stats)
stats <- rbindlist(stats_list, fill = TRUE)
auc_list <- lapply(split(all_d, interaction(all_d$cohort, all_d$mapping_method, drop = TRUE)), function(d) {
  cbind(
    cohort = unique(d$cohort), endpoint = unique(d$endpoint),
    mapping_method = unique(d$mapping_method), eligible_auc(d)
  )
})
auc <- rbindlist(auc_list, fill = TRUE)
fwrite(stats, file.path(RESULT_DIR, "cohort_performance.csv"))
fwrite(auc, file.path(RESULT_DIR, "time_dependent_auc.csv"))
print(stats)

# Patient-level bootstrap intervals for Harrell's C-index. Scores are held
# fixed; resampling quantifies sampling uncertainty within each evaluated
# cohort rather than refitting the development pipeline.
bootstrap_cindex <- function(
    d,
    b = bootstrap_reps,
    seed = 1L,
    evaluation_type = "external_fixed_score") {
  d <- as.data.frame(d)
  required <- c("sample", "time", "status", "risk_sd")
  if (!all(required %in% names(d)) || any(!complete.cases(d[, required, drop = FALSE])) ||
      anyDuplicated(as.character(d$sample)) || any(!is.finite(d$time)) ||
      any(d$time <= 0) || any(!d$status %in% 0:1) ||
      any(!is.finite(d$risk_sd)) || sd(d$risk_sd) <= 0 || sum(d$status) < 1L) {
    stop("C-index bootstrap requires unique patients and complete finite survival scores")
  }
  set.seed(seed)
  n <- nrow(d)
  observed <- concordance(
    Surv(time, status) ~ risk_sd,
    data = d,
    reverse = TRUE,
    timewt = "n"
  )$concordance
  if (!is.finite(observed)) stop("Observed Harrell C-index is non-finite")
  boot <- replicate(b, {
    idx <- sample.int(n, n, replace = TRUE)
    z <- d[idx, , drop = FALSE]
    if (sum(z$status) == 0L) return(NA_real_)
    tryCatch(
      concordance(
        Surv(time, status) ~ risk_sd,
        data = z,
        reverse = TRUE,
        timewt = "n"
      )$concordance,
      error = function(e) NA_real_
    )
  })
  boot <- boot[is.finite(boot)]
  min_valid <- ceiling(0.95 * b)
  if (length(boot) < min_valid) {
    stop(
      "C-index bootstrap produced only ", length(boot), "/", b,
      " finite replicates; required at least ", min_valid
    )
  }
  data.frame(
    n = n, events = sum(d$status), c_index = observed,
    lower95 = unname(quantile(boot, 0.025, na.rm = TRUE)),
    upper95 = unname(quantile(boot, 0.975, na.rm = TRUE)),
    interval_method = "patient-level percentile bootstrap conditional on precomputed scores",
    evaluation_type = evaluation_type,
    harrell_timewt = "n",
    seed = as.integer(seed),
    bootstrap_replicates_requested = b,
    bootstrap_replicates_valid = length(boot),
    bootstrap_valid_fraction = length(boot) / b,
    minimum_valid_fraction = 0.95,
    stringsAsFactors = FALSE
  )
}

primary_data <- all_d[
  (cohort == "GSE39582 (training)" & mapping_method == PRIMARY_GPL570_METHOD) |
    (cohort == "TCGA-COAD" & mapping_method == "gene_level_input") |
    (cohort %in% geo_names & mapping_method == PRIMARY_GPL570_METHOD)
]
primary_split <- split(primary_data, primary_data$cohort)
cindex_ci <- rbindlist(lapply(seq_along(primary_split), function(i) {
  d <- primary_split[[i]]
  cbind(
    cohort = unique(d$cohort), endpoint = unique(d$endpoint),
    mapping_method = unique(d$mapping_method),
    preprocessing_scope = unique(d$preprocessing_scope),
    estimand = unique(d$estimand),
    bootstrap_cindex(
      d,
      b = bootstrap_reps,
      seed = 20240917L + i,
      evaluation_type = if (unique(d$cohort) == "GSE39582 (training)") {
        "apparent_training_fixed_score"
      } else {
        "external_cohort_fixed_score"
      }
    )
  )
}), fill = TRUE)
fwrite(cindex_ci, file.path(RESULT_DIR, "cindex_bootstrap_intervals.csv"))

# A single patient-level summary of repeated OOF predictions. The percentile
# interval conditions on the already generated OOF scores and does not include
# uncertainty from rerunning the full feature-selection algorithm.
nested_oof_dt <- as.data.table(nested_oof)
nested_oof_ensemble <- summarize_repeated_oof(
  nested_oof_dt,
  expected_repeats = PIPELINE_PARAMETERS$outer_repeats
)
nested_oof_ensemble$time <- as.numeric(nested_oof_ensemble$time)
nested_oof_ensemble$status <- as.numeric(nested_oof_ensemble$status)
nested_oof_ensemble$risk_sd <- nested_oof_ensemble$lp_oof_mean
ensemble_cindex <- bootstrap_cindex(
  nested_oof_ensemble,
  b = bootstrap_reps,
  seed = 20240930L,
  evaluation_type = "repeated_nested_cv_oof_ensemble_fixed_predictions"
)
ensemble_cindex$method <- paste0(
  "Patient bootstrap conditional on mean of ",
  PIPELINE_PARAMETERS$outer_repeats,
  " precomputed OOF predictions"
)
fwrite(nested_oof_ensemble, file.path(RESULT_DIR, "nested_cv_oof_ensemble.csv"))
fwrite(
  ensemble_cindex,
  file.path(RESULT_DIR, "nested_cv_oof_ensemble_cindex_bootstrap.csv")
)
internal_performance <- data.frame(
  apparent_c_index = unname(concordance(final_fit)$concordance),
  median_repeat_oof_c_index = nested_summary$median_repeat_oof_c_index,
  repeat_oof_q1 = nested_summary$repeat_oof_c_index_q1,
  repeat_oof_q3 = nested_summary$repeat_oof_c_index_q3,
  repeat_oof_min = nested_summary$min_repeat_oof_c_index,
  repeat_oof_max = nested_summary$max_repeat_oof_c_index,
  ensemble_oof_c_index = ensemble_cindex$c_index,
  ensemble_oof_lower95_conditional = ensemble_cindex$lower95,
  ensemble_oof_upper95_conditional = ensemble_cindex$upper95,
  apparent_to_median_repeat_gap = unname(concordance(final_fit)$concordance) -
    nested_summary$median_repeat_oof_c_index,
  stringsAsFactors = FALSE
)
fwrite(internal_performance, file.path(RESULT_DIR, "internal_performance_summary.csv"))

# Diagnostic check of the linear-score assumption on the log-hazard scale.
# The spline comparison is exploratory, especially in cohorts with few events.
linearity_rows <- lapply(primary_split, function(d) {
  linear_fit <- coxph(Surv(time, status) ~ risk_sd, data = d, ties = "efron")
  spline_fit <- coxph(Surv(time, status) ~ pspline(risk_sd, df = 3), data = d, ties = "efron")
  cmp <- anova(linear_fit, spline_fit, test = "LRT")
  p_col <- grep("Pr|P\\(", colnames(cmp), value = TRUE)
  p_nonlinearity <- if (length(p_col)) as.numeric(cmp[2, p_col[1]]) else NA_real_
  data.frame(
    cohort = unique(d$cohort), endpoint = unique(d$endpoint),
    n = nrow(d), events = sum(d$status),
    linear_loglik = as.numeric(logLik(linear_fit)),
    spline_loglik = as.numeric(logLik(spline_fit)),
    p_nonlinearity = p_nonlinearity,
    method = "LRT: linear Cox vs pspline(df=3)",
    stringsAsFactors = FALSE
  )
})
linearity_df <- rbindlist(linearity_rows, fill = TRUE)
fwrite(linearity_df, file.path(RESULT_DIR, "risk_score_linearity_test.csv"))

# Test the proportional-hazards assumption for the transported risk score in
# every main-analysis cohort. This is distinct from the gene-level test above.
risk_ph_rows <- lapply(
  split(all_d, interaction(all_d$cohort, all_d$mapping_method, drop = TRUE)),
  function(d) {
    fit <- coxph(Surv(time, status) ~ risk_sd, data = d, x = TRUE)
    zph <- cox.zph(fit)
    data.frame(
      cohort = unique(d$cohort), endpoint = unique(d$endpoint),
      mapping_method = unique(d$mapping_method), n = nrow(d),
      events = sum(d$status), chisq = unname(zph$table["risk_sd", "chisq"]),
      p = unname(zph$table["risk_sd", "p"]), stringsAsFactors = FALSE
    )
  }
)
risk_ph <- rbindlist(risk_ph_rows, fill = TRUE)
fwrite(risk_ph, file.path(RESULT_DIR, "risk_score_proportional_hazards_test.csv"))

# Multivariable clinical adjustment where age and stage are usable.
stage_binary <- function(x) {
  out <- rep(NA_character_, length(x))
  z <- trimws(toupper(as.character(x)))
  # TNM stage 0/I/II versus III/IV.
  out[z %in% c("0", "1", "2")] <- "early"
  out[z %in% c("3", "4")] <- "advanced"
  out[grepl("^STAGE\\s*(0|I([^I]|$)|II([^I]|$))", z)] <- "early"
  out[grepl("^STAGE\\s*(III|IV)", z)] <- "advanced"
  # GSE14333 reports Dukes stage rather than TNM stage.
  out[z %in% c("A", "B")] <- "early"
  out[z %in% c("C", "D")] <- "advanced"
  factor(out, levels = c("early", "advanced"))
}
adjusted_rows <- list()
incremental_rows <- list()
clinical_ph_rows <- list()
clinical_eligibility_rows <- list()
for (cohort_name in c("GSE39582 (training)", "TCGA-COAD", "GSE14333", "GSE17536", "GSE17537")) {
  d_source <- all_d[
    all_d$cohort == cohort_name &
      ((cohort_name == "TCGA-COAD" & all_d$mapping_method == "gene_level_input") |
       (cohort_name != "TCGA-COAD" & all_d$mapping_method == PRIMARY_GPL570_METHOD)),
  ]
  d_source$stage_group <- stage_binary(d_source$stage)
  d_source$age10 <- suppressWarnings(as.numeric(d_source$age)) / 10
  complete_idx <- complete.cases(
    d_source[, c("time", "status", "risk_sd", "age10", "stage_group")]
  )
  d <- d_source[complete_idx, ]
  ineligibility_reasons <- character()
  if (nrow(d) < run_parameters$survival$clinical_min_n) {
    ineligibility_reasons <- c(
      ineligibility_reasons,
      paste0("complete_case_n_below_", run_parameters$survival$clinical_min_n)
    )
  }
  if (sum(d$status) < run_parameters$survival$clinical_min_events) {
    ineligibility_reasons <- c(
      ineligibility_reasons,
      paste0("complete_case_events_below_", run_parameters$survival$clinical_min_events)
    )
  }
  if (nrow(d) && length(unique(d$status)) != 2L) {
    ineligibility_reasons <- c(ineligibility_reasons, "outcome_lacks_both_states")
  }
  if (nrow(d) && length(unique(d$stage_group)) != 2L) {
    ineligibility_reasons <- c(ineligibility_reasons, "stage_strata_lack_both_groups")
  }
  if (nrow(d) && (!is.finite(sd(d$age10)) || sd(d$age10) <= 0)) {
    ineligibility_reasons <- c(ineligibility_reasons, "age_has_zero_or_invalid_variance")
  }
  if (nrow(d) && (!is.finite(sd(d$risk_sd)) || sd(d$risk_sd) <= 0)) {
    ineligibility_reasons <- c(ineligibility_reasons, "risk_has_zero_or_invalid_variance")
  }
  clinically_eligible <- !length(ineligibility_reasons)
  clinical_eligibility_rows[[cohort_name]] <- data.frame(
    cohort = cohort_name,
    endpoint = if (nrow(d_source)) unique(d_source$endpoint)[[1L]] else NA_character_,
    mapping_method = if (nrow(d_source)) unique(d_source$mapping_method)[[1L]] else NA_character_,
    source_n = nrow(d_source),
    source_events = sum(d_source$status),
    missing_age_n = sum(is.na(d_source$age10)),
    missing_stage_n = sum(is.na(d_source$stage_group)),
    complete_case_n = nrow(d),
    complete_case_events = sum(d$status),
    excluded_from_adjustment_n = nrow(d_source) - nrow(d),
    excluded_events = sum(d_source$status) - sum(d$status),
    complete_case_fraction = if (nrow(d_source)) nrow(d) / nrow(d_source) else NA_real_,
    stage_groups_present = length(unique(d$stage_group)),
    minimum_n_required = run_parameters$survival$clinical_min_n,
    minimum_events_required = run_parameters$survival$clinical_min_events,
    eligible = clinically_eligible,
    reason = if (clinically_eligible) "eligible" else paste(ineligibility_reasons, collapse = ";"),
    adjustment = "age_per_10y + stage_strata(early/advanced)",
    stringsAsFactors = FALSE
  )
  if (!clinically_eligible) next
  # Stage violated proportional hazards in TCGA. Stratification adjusts for
  # early/advanced stage without imposing a constant stage hazard ratio and is
  # applied consistently across cohorts.
  f0 <- coxph(Surv(time, status) ~ age10 + strata(stage_group), data = d, x = TRUE)
  f <- coxph(Surv(time, status) ~ risk_sd + age10 + strata(stage_group), data = d, x = TRUE)
  sf <- summary(f)
  term_labels <- rownames(sf$coefficients)
  term_labels[term_labels == "age10"] <- "age_per_10y"
  adjusted_rows[[cohort_name]] <- data.frame(
    cohort = cohort_name, endpoint = unique(d$endpoint), n = nrow(d), events = sum(d$status),
    adjustment = "age_per_10y + stage_strata(early/advanced)",
    term = term_labels, coef = sf$coefficients[, "coef"],
    HR = sf$conf.int[, "exp(coef)"], lower95 = sf$conf.int[, "lower .95"],
    upper95 = sf$conf.int[, "upper .95"], p = sf$coefficients[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
  lrt <- anova(f0, f, test = "LRT")
  incremental_rows[[cohort_name]] <- data.frame(
    cohort = cohort_name, endpoint = unique(d$endpoint), n = nrow(d), events = sum(d$status),
    clinical_c_index = concordance(f0)$concordance,
    clinical_plus_risk_c_index = concordance(f)$concordance,
    apparent_delta_c_index = concordance(f)$concordance - concordance(f0)$concordance,
    likelihood_ratio_chisq = lrt[2, "Chisq"], likelihood_ratio_p = lrt[2, "Pr(>|Chi|)"],
    stringsAsFactors = FALSE
  )
  zph <- cox.zph(f)
  clinical_ph_rows[[cohort_name]] <- data.frame(
    cohort = cohort_name, endpoint = unique(d$endpoint),
    term = rownames(zph$table), chisq = zph$table[, "chisq"],
    p = zph$table[, "p"], stringsAsFactors = FALSE
  )
}
clinical_eligibility <- rbindlist(clinical_eligibility_rows, fill = TRUE)
fwrite(
  clinical_eligibility,
  file.path(RESULT_DIR, "clinical_adjustment_eligibility.csv")
)
adjusted_df <- rbindlist(adjusted_rows, fill = TRUE)
fwrite(adjusted_df, file.path(RESULT_DIR, "multivariable_clinical_cox.csv"))
fwrite(rbindlist(incremental_rows, fill = TRUE), file.path(RESULT_DIR, "clinical_incremental_performance.csv"))
fwrite(rbindlist(clinical_ph_rows, fill = TRUE), file.path(RESULT_DIR, "clinical_model_proportional_hazards_test.csv"))

# Small-k random-effects synthesis of OS validation cohorts. Knapp-Hartung/t
# inference is primary; conventional Wald inference is retained transparently.
os_val <- as.data.frame(stats[
  ((stats$cohort == "TCGA-COAD" & stats$mapping_method == "gene_level_input") |
   (stats$cohort != "TCGA-COAD" & stats$mapping_method == PRIMARY_GPL570_METHOD)) &
    stats$endpoint == "OS" & stats$cohort %in% c("TCGA-COAD", "GSE17536", "GSE17537"),
])
os_val$yi <- log(os_val$HR_per_SD)
os_val$sei <- (log(os_val$upper95) - log(os_val$lower95)) / (2 * qnorm(0.975))
os_meta_inputs <- os_val[, c(
  "cohort", "endpoint", "mapping_method", "n", "events", "HR_per_SD",
  "lower95", "upper95", "cox_p", "yi", "sei"
)]
os_meta_inputs$effect_scale <- "log hazard ratio per 1-SD risk score"
os_meta_inputs$sei_method <- "derived from cohort-level 95% Wald CI"
fwrite(
  os_meta_inputs,
  file.path(RESULT_DIR, "os_validation_meta_inputs.csv")
)
meta_summary <- fit_reml_meta_dual(
  os_val,
  analysis = "unadjusted external OS association",
  endpoint = "OS"
)
fwrite(meta_summary, file.path(RESULT_DIR, "os_validation_meta_analysis.csv"))
print(meta_summary)

leave_one_out <- leave_one_out_reml_meta_dual(
  os_val,
  analysis = "leave-one-out unadjusted external OS association",
  endpoint = "OS"
)
fwrite(leave_one_out, file.path(RESULT_DIR, "os_validation_meta_leave_one_out.csv"))

# Apply the identical REML Knapp-Hartung/t primary and Wald sensitivity
# meta-analysis contract to the strict full-pipeline highest-mean sensitivity.
highest_os_val <- as.data.frame(highest_stats[
  highest_stats$endpoint == "OS" &
    highest_stats$cohort %in% c("TCGA-COAD", "GSE17536", "GSE17537"),
])
if (nrow(highest_os_val) != 3L ||
    anyDuplicated(highest_os_val$cohort) ||
    !setequal(
      highest_os_val$cohort,
      c("TCGA-COAD", "GSE17536", "GSE17537")
    )) {
  stop("Strict highest-mean OS meta-analysis lacks three external cohorts")
}
highest_os_val$yi <- log(highest_os_val$HR_per_SD)
highest_os_val$sei <- (
  log(highest_os_val$upper95) - log(highest_os_val$lower95)
) / (2 * qnorm(0.975))
highest_os_meta_inputs <- highest_os_val[, c(
  "cohort", "endpoint", "mapping_method", "n", "events", "HR_per_SD",
  "lower95", "upper95", "cox_p", "yi", "sei", "probe_map_hash",
  "model_fingerprint", "analysis_key"
)]
highest_os_meta_inputs$effect_scale <-
  "log hazard ratio per 1-SD cohort-standardized risk score"
highest_os_meta_inputs$sei_method <-
  "derived from cohort-level 95% Wald CI"
fwrite(
  highest_os_meta_inputs,
  file.path(RESULT_DIR, "highest_mean_os_validation_meta_inputs.csv")
)
highest_meta_summary <- fit_reml_meta_dual(
  highest_os_val,
  analysis = "strict highest-mean full-pipeline external OS association",
  endpoint = "OS"
)
highest_meta_summary$probe_map_hash <- highest_full_model$probe_map_hash
highest_meta_summary$model_fingerprint <-
  highest_full_model$model_fingerprint
highest_meta_summary$analysis_key <- highest_mean_analysis_key
fwrite(
  highest_meta_summary,
  file.path(RESULT_DIR, "highest_mean_os_validation_meta_analysis.csv")
)
highest_leave_one_out <- leave_one_out_reml_meta_dual(
  highest_os_val,
  analysis = paste(
    "leave-one-out strict highest-mean full-pipeline external OS association"
  ),
  endpoint = "OS"
)
highest_leave_one_out$probe_map_hash <- highest_full_model$probe_map_hash
highest_leave_one_out$model_fingerprint <-
  highest_full_model$model_fingerprint
highest_leave_one_out$analysis_key <- highest_mean_analysis_key
fwrite(
  highest_leave_one_out,
  file.path(
    RESULT_DIR, "highest_mean_os_validation_meta_leave_one_out.csv"
  )
)

adjusted_os <- as.data.frame(adjusted_df[
  adjusted_df$endpoint == "OS" & adjusted_df$term == "risk_sd" &
    adjusted_df$cohort %in% c("TCGA-COAD", "GSE17536", "GSE17537"),
])
if (nrow(adjusted_os) == 3L) {
  adjusted_os$yi <- log(adjusted_os$HR)
  adjusted_os$sei <- (log(adjusted_os$upper95) - log(adjusted_os$lower95)) / (2 * qnorm(0.975))
  adjusted_meta_inputs <- adjusted_os[, c(
    "cohort", "endpoint", "n", "events", "adjustment", "HR",
    "lower95", "upper95", "p", "yi", "sei"
  )]
  adjusted_meta_inputs$effect_scale <- "log hazard ratio per 1-SD risk score"
  adjusted_meta_inputs$sei_method <- "derived from cohort-level 95% Wald CI"
  fwrite(
    adjusted_meta_inputs,
    file.path(RESULT_DIR, "os_adjusted_meta_inputs.csv")
  )
  adjusted_meta <- fit_reml_meta_dual(
    adjusted_os,
    analysis = "age-adjusted, stage-stratified cohort estimates",
    endpoint = "OS"
  )
  fwrite(adjusted_meta, file.path(RESULT_DIR, "os_adjusted_meta_analysis.csv"))
} else {
  for (artifact_id in optional_artifacts) {
    register_not_generated(
      artifact_registry,
      artifact_id,
      paste0(
        "fewer_than_3_eligible_cohorts: eligible=", nrow(adjusted_os),
        "; required=3"
      )
    )
  }
}

# Full-pipeline probe-mapping sensitivity comparison. The sensitivity rows
# come from a separately refitted strict nested pipeline and a single frozen
# full-development probe map, never from cohort-adaptive probe reselection.
sensitivity_cohort_order <- data.frame(
  cohort = c("TCGA-COAD", "GSE14333", "GSE17536", "GSE17537"),
  endpoint = c("OS", "DFS", "OS", "OS"),
  figure_order = 1:4,
  stringsAsFactors = FALSE
)
for (nm in sensitivity_cohort_order$cohort) {
  primary_identity <- as.data.frame(all_d[all_d$cohort == nm, ])
  strict_identity <- as.data.frame(
    highest_external_d[highest_external_d$cohort == nm, ]
  )
  primary_identity <- primary_identity[
    order(primary_identity$sample), c("sample", "time", "status"), drop = FALSE
  ]
  strict_identity <- strict_identity[
    order(strict_identity$sample), c("sample", "time", "status"), drop = FALSE
  ]
  if (anyDuplicated(primary_identity$sample) ||
      anyDuplicated(strict_identity$sample) ||
      !identical(primary_identity$sample, strict_identity$sample) ||
      !identical(
        as.numeric(primary_identity$time),
        as.numeric(strict_identity$time)
      ) ||
      !identical(
        as.numeric(primary_identity$status),
        as.numeric(strict_identity$status)
      )) {
    stop(
      nm,
      " primary and strict sensitivity pipelines differ in sample/time/status"
    )
  }
}
primary_sens <- as.data.frame(stats[
  (stats$cohort == "TCGA-COAD" &
     stats$mapping_method == "gene_level_input") |
    (stats$cohort %in% geo_names &
       stats$mapping_method == PRIMARY_GPL570_METHOD),
])
primary_sens$analysis_pipeline <- "primary_unique_mean_full_pipeline"
primary_mapping_hashes <- setNames(vapply(
  sensitivity_cohort_order$cohort,
  function(nm) {
    mapping_object <- if (nm == "TCGA-COAD") {
      tcga_gene_resolution$audit
    } else {
      key <- paste(nm, PRIMARY_GPL570_METHOD, sep = "_")
      if (is.null(probe_gene_mappings[[key]])) {
        stop("Missing primary mapping identity for ", nm)
      }
      probe_gene_mappings[[key]]
    }
    digest(
      canonicalize_for_hash(mapping_object),
      algo = "sha256", serialize = TRUE
    )
  },
  character(1)
), sensitivity_cohort_order$cohort)
primary_sens$probe_map_hash <- unname(
  primary_mapping_hashes[primary_sens$cohort]
)
primary_sens$model_fingerprint <- full_pipeline$model_fingerprint
primary_sens$analysis_key <- analysis_key
strict_sens <- as.data.frame(highest_stats)
strict_sens$analysis_pipeline <- "strict_highest_mean_full_pipeline"
sens <- rbindlist(list(primary_sens, strict_sens), fill = TRUE)
sens_counts <- table(sens$cohort, sens$analysis_pipeline)
if (!setequal(
    rownames(sens_counts),
    c("TCGA-COAD", "GSE14333", "GSE17536", "GSE17537")) ||
    any(sens_counts != 1L)) {
  stop("Primary-versus-strict sensitivity table is incomplete")
}
identity_fields <- c("analysis_key", "probe_map_hash", "model_fingerprint")
if (anyNA(sens[, ..identity_fields]) ||
    any(!nzchar(unlist(sens[, ..identity_fields], use.names = FALSE)))) {
  stop("Primary-versus-strict sensitivity identity fields are incomplete")
}
sens_effect_wide <- dcast(
  sens, cohort + endpoint + n + events ~ analysis_pipeline,
  value.var = c("HR_per_SD", "lower95", "upper95", "cox_p", "c_index")
)
sens_identity_wide <- dcast(
  sens, cohort + endpoint + n + events ~ analysis_pipeline,
  value.var = identity_fields
)
sens_wide <- merge(
  sens_effect_wide, sens_identity_wide,
  by = c("cohort", "endpoint", "n", "events"), sort = FALSE
)
sens_wide$figure_order <- sensitivity_cohort_order$figure_order[
  match(sens_wide$cohort, sensitivity_cohort_order$cohort)
]
if (anyNA(sens_wide$figure_order)) {
  stop("Sensitivity output contains an unexpected cohort")
}
sens_wide <- sens_wide[order(sens_wide$figure_order), ]
sens_wide$figure_order <- NULL
fwrite(sens_wide, file.path(RESULT_DIR, "probe_mapping_sensitivity.csv"))

# -------------------------------------------------------------------------
# Publication-ready core figures.
# -------------------------------------------------------------------------
plot_stats <- stats[
  (stats$cohort == "GSE39582 (training)" & stats$mapping_method == PRIMARY_GPL570_METHOD) |
    (stats$cohort == "TCGA-COAD" & stats$mapping_method == "gene_level_input") |
    (stats$cohort %in% geo_names & stats$mapping_method == PRIMARY_GPL570_METHOD),
]
forest_cohort_order <- data.frame(
  cohort = c(
    "GSE39582 (training)", "GSE14333", "TCGA-COAD",
    "GSE17536", "GSE17537"
  ),
  endpoint = c("OS", "DFS", "OS", "OS", "OS"),
  analysis_role = c(
    "Development (apparent)", "External DFS (not pooled)",
    rep("External OS", 3)
  ),
  row_order = 1:5,
  stringsAsFactors = FALSE
)
plot_stats <- as.data.frame(plot_stats)
forest_match <- match(forest_cohort_order$cohort, plot_stats$cohort)
if (anyNA(forest_match) || anyDuplicated(plot_stats$cohort) ||
    nrow(plot_stats) != nrow(forest_cohort_order) ||
    !identical(
      as.character(plot_stats$endpoint[forest_match]),
      forest_cohort_order$endpoint
    )) {
  stop("Forest plot cohort/endpoint roles are incomplete or misassigned")
}
plot_stats <- plot_stats[forest_match, , drop = FALSE]
plot_stats$analysis_role <- forest_cohort_order$analysis_role
plot_stats$row_order <- forest_cohort_order$row_order
plot_stats$label <- sprintf(
  "%s [%s]\nn=%d; events=%d",
  plot_stats$cohort, plot_stats$endpoint, plot_stats$n, plot_stats$events
)
pooled_primary <- meta_summary[meta_summary$primary, , drop = FALSE]
if (nrow(pooled_primary) != 1L ||
    any(!is.finite(c(
      pooled_primary$prediction_lower95,
      pooled_primary$prediction_upper95
    ))) || pooled_primary$prediction_lower95 <= 0) {
  stop("Primary KH meta-analysis lacks a finite positive prediction interval")
}
pooled_plot <- data.frame(
  cohort = "Pooled external OS",
  endpoint = "OS",
  mapping_method = PRIMARY_GPL570_METHOD,
  n = sum(os_val$n),
  events = sum(os_val$events),
  HR_per_SD = pooled_primary$pooled_HR,
  lower95 = pooled_primary$lower95,
  upper95 = pooled_primary$upper95,
  prediction_lower95 = pooled_primary$prediction_lower95,
  prediction_upper95 = pooled_primary$prediction_upper95,
  analysis_role = "Pooled external OS (REML-KH)",
  row_order = 6L,
  label = sprintf(
    "Pooled external OS [OS]\nk=%d; n=%d; events=%d",
    nrow(os_val), sum(os_val$n), sum(os_val$events)
  ),
  stringsAsFactors = FALSE
)
forest_data <- rbindlist(list(
  plot_stats[, c(
    "cohort", "endpoint", "mapping_method", "n", "events", "HR_per_SD",
    "lower95", "upper95", "analysis_role", "row_order", "label"
  )],
  pooled_plot
), fill = TRUE)
forest_data$analysis_role <- factor(
  forest_data$analysis_role,
  levels = c(
    "Development (apparent)", "External DFS (not pooled)",
    "External OS", "Pooled external OS (REML-KH)"
  )
)
forest_data$label <- factor(
  forest_data$label,
  levels = rev(as.character(forest_data$label[order(forest_data$row_order)]))
)
forest_nonpooled <- forest_data[forest_data$cohort != "Pooled external OS", ]
forest_pooled <- forest_data[forest_data$cohort == "Pooled external OS", ]
p_forest <- ggplot(
  forest_data, aes(HR_per_SD, label)
) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
  geom_errorbar(
    data = forest_pooled,
    aes(xmin = prediction_lower95, xmax = prediction_upper95),
    orientation = "y", width = 0.22, linewidth = 1.05,
    linetype = "dashed", color = "grey45"
  ) +
  geom_errorbar(
    aes(xmin = lower95, xmax = upper95),
    orientation = "y", width = 0.16, linewidth = 0.7
  ) +
  geom_point(
    data = forest_nonpooled,
    aes(color = analysis_role, shape = analysis_role), size = 2.5
  ) +
  geom_point(
    data = forest_pooled,
    shape = 23, size = 4.0, fill = "#222222", color = "#222222"
  ) +
  facet_grid(
    rows = vars(analysis_role), scales = "free_y", space = "free_y",
    switch = "y"
  ) +
  scale_x_log10() +
  scale_color_manual(values = c(
    "Development (apparent)" = "#666666",
    "External DFS (not pooled)" = "#2166AC",
    "External OS" = "#B2182B"
  )) +
  scale_shape_manual(values = c(
    "Development (apparent)" = 15,
    "External DFS (not pooled)" = 17,
    "External OS" = 16
  )) +
  labs(
    x = "Hazard ratio per 1-SD higher risk score (95% CI)",
    y = NULL, color = NULL, shape = NULL,
    title = "Cohort associations by prespecified analysis role",
    caption = sprintf(
      paste(
        "Development is apparent; GSE14333 reports DFS and is not pooled; external OS is pooled by REML-KH.",
        "\nKH 95%% prediction interval %.2f-%.2f."
      ),
      pooled_primary$prediction_lower95,
      pooled_primary$prediction_upper95
    )
  ) +
  theme_pub(8) + theme(
    text = element_text(family = "Arial"),
    legend.position = "none",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7, face = "bold"),
    panel.spacing.y = grid::unit(0.12, "lines"),
    plot.title = element_text(size = 9, face = "bold"),
    plot.caption = element_text(hjust = 0, size = 6.5)
  )
save_pub_plot(p_forest, file.path(FIGURE_DIR, "Fig2_cohort_forest"), width = 7.2, height = 5.4)

km_external_order <- data.frame(
  cohort = c("TCGA-COAD", "GSE14333", "GSE17536", "GSE17537"),
  endpoint = c("OS", "DFS", "OS", "OS"),
  stringsAsFactors = FALSE
)
km_data <- all_d[
  (all_d$cohort == "TCGA-COAD" & all_d$mapping_method == "gene_level_input") |
    (all_d$cohort %in% geo_names & all_d$mapping_method == PRIMARY_GPL570_METHOD),
]
if (!setequal(unique(km_data$cohort), km_external_order$cohort)) {
  stop("Main KM figure does not contain exactly the four external cohorts")
}
km_endpoints <- vapply(km_external_order$cohort, function(nm) {
  z <- unique(as.character(km_data$endpoint[km_data$cohort == nm]))
  if (length(z) != 1L) stop("Non-unique KM endpoint for ", nm)
  z
}, character(1))
if (!identical(unname(km_endpoints), km_external_order$endpoint)) {
  stop("Main KM figure cohort/endpoint order is incorrect")
}
km_max_year <- max(km_data$time / 365.25)
km_break_by <- if (km_max_year <= 6) 1 else if (km_max_year <= 12) 2 else 5
km_x_limit <- ceiling(km_max_year / km_break_by) * km_break_by
km_x_breaks <- seq(0, km_x_limit, by = km_break_by)
km_curve_rows <- list()
risk_table_rows <- list()
for (nm in km_external_order$cohort) {
  d <- km_data[km_data$cohort == nm, ]
  fit <- survfit(Surv(time / 365.25, status) ~ group, data = d)
  # `censored = TRUE` is required here: the default summary reports event
  # times only and would place/omit censor marks at the wrong x-coordinates.
  sm <- summary(fit, censored = TRUE)
  curve_rows <- data.frame(
    time_years = sm$time, survival = sm$surv, lower = sm$lower, upper = sm$upper,
    n_censor = sm$n.censor,
    group = sub("group=", "", sm$strata), cohort = nm,
    endpoint = unique(d$endpoint)
  )
  risk_groups <- levels(droplevels(d$group))
  if (!identical(risk_groups, c("Low", "High"))) {
    stop("KM risk groups are not the prespecified Low/High pair for ", nm)
  }
  # survfit summaries begin at the first observed event/censor time. Add the
  # mathematical KM origin explicitly so every displayed curve starts at
  # S(0) = 1 instead of appearing to begin at its first observed time.
  baseline_rows <- data.frame(
    time_years = 0, survival = 1, lower = 1, upper = 1,
    n_censor = 0L, group = risk_groups, cohort = nm,
    endpoint = unique(d$endpoint), stringsAsFactors = FALSE
  )
  km_curve_rows[[nm]] <- rbind(baseline_rows, curve_rows)
  for (risk_group in risk_groups) {
    group_time <- d$time[d$group == risk_group] / 365.25
    risk_table_rows[[paste(nm, risk_group, sep = "::")]] <- data.frame(
      cohort = nm,
      endpoint = unique(d$endpoint),
      mapping_method = unique(d$mapping_method),
      time_years = km_x_breaks,
      group = risk_group,
      n_at_risk = vapply(
        km_x_breaks,
        function(tt) sum(group_time >= tt),
        integer(1)
      ),
      stringsAsFactors = FALSE
    )
  }
}
km_curve <- rbindlist(km_curve_rows)
km_curve$group <- factor(km_curve$group, levels = c("Low", "High"))
km_risk_table <- rbindlist(risk_table_rows, fill = TRUE)
km_risk_table$group <- factor(km_risk_table$group, levels = c("Low", "High"))
fwrite(km_risk_table, file.path(RESULT_DIR, "km_risk_table.csv"))
ann <- plot_stats[
  match(km_external_order$cohort, plot_stats$cohort),
  c("cohort", "endpoint", "n", "events", "HR_per_SD", "lower95", "upper95", "logrank_p")
]
ann$text <- sprintf(
  "n=%d; events=%d\nHR/SD %.2f (%.2f-%.2f)\nlog-rank p=%s",
  ann$n, ann$events, ann$HR_per_SD, ann$lower95, ann$upper95, p_format(ann$logrank_p)
)
km_curve_plots <- list()
km_risk_plots <- list()
for (i in seq_len(nrow(km_external_order))) {
  nm <- km_external_order$cohort[[i]]
  endpoint_i <- km_external_order$endpoint[[i]]
  curve_i <- km_curve[km_curve$cohort == nm, ]
  risk_i <- km_risk_table[km_risk_table$cohort == nm, ]
  ann_i <- ann[ann$cohort == nm, , drop = FALSE]
  km_curve_plots[[nm]] <- ggplot(
    curve_i,
    aes(time_years, survival, color = group, linetype = group)
  ) +
    geom_step(linewidth = 0.72) +
    geom_point(
      data = curve_i[curve_i$n_censor > 0, ],
      aes(shape = group), size = 0.8, stroke = 0.4,
      show.legend = FALSE
    ) +
    geom_text(
      data = ann_i,
      aes(x = 0.02 * km_x_limit, y = 0.05, label = text),
      inherit.aes = FALSE, hjust = 0, vjust = 0,
      size = 2.05, family = "Arial", color = "#222222"
    ) +
    scale_color_manual(values = c(Low = "#2166AC", High = "#B2182B")) +
    scale_linetype_manual(values = c(Low = "solid", High = "22")) +
    scale_shape_manual(values = c(Low = 3, High = 4)) +
    scale_x_continuous(
      limits = c(0, km_x_limit), breaks = km_x_breaks,
      expand = expansion(mult = c(0, 0.01))
    ) +
    scale_y_continuous(
      limits = c(0, 1), breaks = seq(0, 1, 0.25),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = paste0(nm, " - ", endpoint_i),
      x = NULL, y = "Survival probability", tag = LETTERS[[i]]
    ) +
    theme_pub(8) +
    theme(
      text = element_text(family = "Arial"), legend.position = "none",
      plot.title = element_text(face = "bold", size = 8),
      plot.tag = element_text(face = "bold", size = 9),
      plot.margin = margin(4, 5, 0, 5)
    )
  km_risk_plots[[nm]] <- ggplot(
    risk_i,
    aes(time_years, group, label = n_at_risk, color = group)
  ) +
    geom_text(
      data = risk_i[risk_i$time_years > 0, , drop = FALSE],
      size = 2.15, family = "Arial", show.legend = FALSE
    ) +
    geom_text(
      data = risk_i[risk_i$time_years == 0, , drop = FALSE],
      aes(x = time_years + 0.12), hjust = 0,
      size = 2.15, family = "Arial", show.legend = FALSE
    ) +
    scale_color_manual(values = c(Low = "#2166AC", High = "#B2182B")) +
    scale_x_continuous(
      limits = c(0, km_x_limit), breaks = km_x_breaks,
      expand = expansion(mult = c(0, 0.01))
    ) +
    scale_y_discrete(drop = FALSE) +
    labs(x = "Time (years)", y = "At risk") +
    theme_pub(7) +
    theme(
      text = element_text(family = "Arial"), legend.position = "none",
      axis.line.y = element_blank(), axis.ticks.y = element_blank(),
      plot.margin = margin(0, 5, 4, 5)
    )
}

draw_km_with_risk_table <- function() {
  grid::grid.newpage()
  outer_layout <- grid::grid.layout(
    nrow = 3L, ncol = 2L,
    heights = grid::unit.c(
      grid::unit(0.22, "inches"),
      grid::unit(1, "null"), grid::unit(1, "null")
    )
  )
  grid::pushViewport(grid::viewport(layout = outer_layout))
  grid::pushViewport(grid::viewport(layout.pos.row = 1L, layout.pos.col = 1:2))
  grid::grid.text(
    "Risk groups", x = 0.26, y = 0.52,
    gp = grid::gpar(fontfamily = "Arial", fontsize = 7, fontface = "bold")
  )
  grid::grid.segments(
    x0 = 0.36, x1 = 0.43, y0 = 0.52, y1 = 0.52,
    gp = grid::gpar(col = "#2166AC", lwd = 1.4, lty = 1)
  )
  grid::grid.text(
    "Low", x = 0.47, y = 0.52,
    gp = grid::gpar(fontfamily = "Arial", fontsize = 7)
  )
  grid::grid.segments(
    x0 = 0.55, x1 = 0.62, y0 = 0.52, y1 = 0.52,
    gp = grid::gpar(col = "#B2182B", lwd = 1.4, lty = 2)
  )
  grid::grid.text(
    "High", x = 0.67, y = 0.52,
    gp = grid::gpar(fontfamily = "Arial", fontsize = 7)
  )
  grid::popViewport()
  for (i in seq_len(nrow(km_external_order))) {
    nm <- km_external_order$cohort[[i]]
    grid::pushViewport(grid::viewport(
      layout.pos.row = ((i - 1L) %/% 2L) + 2L,
      layout.pos.col = ((i - 1L) %% 2L) + 1L
    ))
    pair_layout <- grid::grid.layout(
      nrow = 2L, ncol = 1L,
      heights = grid::unit(c(3.2, 1.0), "null")
    )
    grid::pushViewport(grid::viewport(layout = pair_layout))
    print(
      km_curve_plots[[nm]],
      vp = grid::viewport(layout.pos.row = 1L, layout.pos.col = 1L)
    )
    print(
      km_risk_plots[[nm]],
      vp = grid::viewport(layout.pos.row = 2L, layout.pos.col = 1L)
    )
    grid::popViewport(2L)
  }
  grid::popViewport()
}
grDevices::cairo_pdf(
  file.path(FIGURE_DIR, "Fig3_multicohort_KM.pdf"),
  width = 7.2, height = 7.4
)
draw_km_with_risk_table()
dev.off()
png(
  file.path(FIGURE_DIR, "Fig3_multicohort_KM.png"),
  width = 7.2, height = 7.4, units = "in", res = 300,
  type = "cairo-png"
)
draw_km_with_risk_table()
dev.off()

sens_long <- as.data.frame(sens)[, c(
  "cohort", "endpoint", "analysis_pipeline", "HR_per_SD", "lower95",
  "upper95", "analysis_key", "probe_map_hash", "model_fingerprint"
), drop = FALSE]
sens_long$figure_order <- sensitivity_cohort_order$figure_order[
  match(sens_long$cohort, sensitivity_cohort_order$cohort)
]
expected_endpoint <- sensitivity_cohort_order$endpoint[
  match(sens_long$cohort, sensitivity_cohort_order$cohort)
]
if (anyNA(sens_long$figure_order) ||
    any(as.character(sens_long$endpoint) != expected_endpoint)) {
  stop("Fig4 cohort/endpoint order is incomplete or incorrect")
}
sens_long$cohort_endpoint <- factor(
  paste0(sens_long$cohort, " [", sens_long$endpoint, "]"),
  levels = rev(paste0(
    sensitivity_cohort_order$cohort, " [",
    sensitivity_cohort_order$endpoint, "]"
  ))
)
probe_method_labels <- c(
  primary_unique_mean_full_pipeline = "Primary unique-mean development pipeline",
  strict_highest_mean_full_pipeline =
    "Strict highest-mean development pipeline"
)
p_sens <- ggplot(
  sens_long,
  aes(
    HR_per_SD, cohort_endpoint, color = analysis_pipeline,
    shape = analysis_pipeline
  )
) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
  geom_errorbar(
    aes(xmin = lower95, xmax = upper95),
    orientation = "y", width = 0.12,
    position = position_dodge(width = 0.35)
  ) +
  geom_point(size = 2.5, position = position_dodge(width = 0.35)) +
  scale_x_log10() +
  scale_color_manual(
    values = c(
      primary_unique_mean_full_pipeline = "#2166AC",
      strict_highest_mean_full_pipeline = "#B2182B"
    ),
    breaks = names(probe_method_labels), labels = unname(probe_method_labels)
  ) +
  scale_shape_manual(
    values = c(
      primary_unique_mean_full_pipeline = 16,
      strict_highest_mean_full_pipeline = 17
    ),
    breaks = names(probe_method_labels), labels = unname(probe_method_labels)
  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  labs(
    title = "External-cohort sensitivity to the development mapping pipeline",
    x = "Hazard ratio per 1-SD risk score (95% CI)",
    y = NULL, color = NULL, shape = NULL,
    caption = paste(
      "Each paired estimate uses identical samples, survival times, and event",
      "indicators; identity hashes remain in probe_mapping_sensitivity.csv."
    )
  ) +
  theme_pub(8) + theme(
    text = element_text(family = "Arial"),
    legend.position = "top",
    legend.box = "vertical",
    plot.title = element_text(size = 9, face = "bold"),
    plot.caption = element_text(size = 6.5, hjust = 0)
  )
save_pub_plot(p_sens, file.path(FIGURE_DIR, "Fig4_probe_mapping_sensitivity"), width = 7.2, height = 4.0)

writeLines(capture.output(sessionInfo()), file.path(RESULT_DIR, "sessionInfo.txt"))
cat("\nCompleted analysis calculations:",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
close_run_log()

output_manifest <- finalize_artifact_registry(
  artifact_registry,
  output_manifest_path,
  control_relative_paths = c(
    file.path("results", "run_status.csv"),
    file.path("results", "output_manifest.csv")
  )
)
output_manifest_sha256 <- sha256_file(output_manifest_path)
write_run_status(
  "complete",
  "All registered artifacts passed hash and allow-list validation"
)
commit_run_directory(RUN_TMP_DIR, RUN_FINAL_DIR, current_run_key)
run_complete <- TRUE
message("Committed run directory: ", RUN_FINAL_DIR)
}

main()
