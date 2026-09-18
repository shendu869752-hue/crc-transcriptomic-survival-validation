# Strict training-partition-only GPL570 highest-mean probe sensitivity.
#
# This file is intentionally separate from the primary unique_mean pipeline.
# It expects utils.R to have been sourced first and never learns a probe choice
# from a validation or external partition.

.probe_map_schema <- "unique_highest_mean_probe_map_v4_radix_textsha256"
.probe_model_schema <- "probe_signature_model_v1"
.probe_sample_schema <- "probe_training_sample_set_v1"

.validate_probe_expression <- function(probe_expr, context = "probe expression",
                                       min_samples = 1L) {
  x <- as.matrix(probe_expr)
  storage.mode(x) <- "double"
  if (nrow(x) < 1L || ncol(x) < as.integer(min_samples)) {
    stop(context, ": insufficient probes or samples")
  }
  if (is.null(rownames(x)) || anyNA(rownames(x)) ||
      any(!nzchar(rownames(x))) || anyDuplicated(rownames(x))) {
    stop(context, ": probe row names must be nonempty and unique")
  }
  if (is.null(colnames(x)) || anyNA(colnames(x)) ||
      any(!nzchar(colnames(x))) || anyDuplicated(colnames(x))) {
    stop(context, ": sample column names must be nonempty and unique")
  }
  if (any(!is.finite(x))) stop(context, ": non-finite expression value")
  x
}

.validate_probe_survival_inputs <- function(probe_expr, time, status,
                                             context = "probe training") {
  x <- .validate_probe_expression(probe_expr, context, min_samples = 2L)
  time <- as.numeric(time)
  status <- as.numeric(status)
  if (ncol(x) != length(time) || length(time) != length(status)) {
    stop(context, ": expression columns, time, and status have different lengths")
  }
  if (any(!is.finite(time)) || any(time <= 0)) {
    stop(context, ": survival times must be finite and positive")
  }
  if (any(!status %in% 0:1) || length(unique(status)) != 2L) {
    stop(context, ": status must contain both 0 and 1 only")
  }
  list(probe_expr = x, time = time, status = status)
}

.validate_probe_annotation <- function(probe_annotation) {
  ann <- as.data.frame(probe_annotation, stringsAsFactors = FALSE)
  required <- c("PROBEID", "SYMBOL")
  missing <- setdiff(required, names(ann))
  if (length(missing)) {
    stop("Probe annotation missing columns: ", paste(missing, collapse = ", "))
  }
  ann <- ann[, required, drop = FALSE]
  ann$PROBEID <- as.character(ann$PROBEID)
  ann$SYMBOL <- as.character(ann$SYMBOL)
  ann$PROBEID <- trimws(ann$PROBEID)
  ann$SYMBOL <- trimws(ann$SYMBOL)
  ann
}

.validate_sample_set <- function(sample_ids, expected = NULL,
                                 context = "sample IDs") {
  sample_ids <- as.character(sample_ids)
  if (!length(sample_ids) || anyNA(sample_ids) || any(!nzchar(sample_ids)) ||
      anyDuplicated(sample_ids)) {
    stop(context, " must be nonempty and unique")
  }
  if (!is.null(expected) && !setequal(sample_ids, as.character(expected))) {
    stop(context, " do not match the expression columns")
  }
  sample_ids
}

.hash_probe_sample_ids <- function(sample_ids) {
  sample_ids <- .validate_sample_set(sample_ids)
  digest::digest(
    list(schema = .probe_sample_schema, sample_ids = sort(sample_ids)),
    algo = "sha256", serialize = TRUE
  )
}

.canonical_probe_mapping <- function(probe_map) {
  mapping <- if (is.list(probe_map) && !is.data.frame(probe_map)) {
    probe_map$mapping
  } else {
    probe_map
  }
  mapping <- as.data.frame(mapping, stringsAsFactors = FALSE)
  required <- c(
    "SYMBOL", "PROBEID", "training_probe_mean", "training_probe_mean_hex",
    "n_probes_gene"
  )
  missing <- setdiff(required, names(mapping))
  if (length(missing)) {
    stop("Probe map missing columns: ", paste(missing, collapse = ", "))
  }
  mapping <- mapping[, required, drop = FALSE]
  mapping$SYMBOL <- as.character(mapping$SYMBOL)
  mapping$PROBEID <- as.character(mapping$PROBEID)
  reported_mean <- as.numeric(mapping$training_probe_mean)
  mapping$training_probe_mean_hex <-
    as.character(mapping$training_probe_mean_hex)
  valid_hex <- grepl(
    "^hex:[+-]?0x[0-9a-f]+(?:\\.[0-9a-f]*)?p[+-]?[0-9]+$",
    mapping$training_probe_mean_hex,
    ignore.case = TRUE,
    perl = TRUE
  )
  exact_mean <- suppressWarnings(as.numeric(sub(
    "^hex:", "", mapping$training_probe_mean_hex
  )))
  if (any(!valid_hex) || any(!is.finite(exact_mean)) ||
      any(!is.finite(reported_mean)) ||
      any(abs(reported_mean - exact_mean) >
          1e-12 * pmax(1, abs(exact_mean)))) {
    stop("Probe map has an invalid or inconsistent exact training mean")
  }
  # Hash the exact value recovered from the persisted hexadecimal text. The
  # ordinary numeric column remains human-readable but is not trusted for a
  # bitwise hash after CSV serialization.
  mapping$training_probe_mean <- exact_mean
  mapping$training_probe_mean_hex <- paste0(
    "hex:", sprintf("%a", exact_mean)
  )
  mapping$n_probes_gene <- as.numeric(mapping$n_probes_gene)
  if (!nrow(mapping) || anyNA(mapping$SYMBOL) || anyNA(mapping$PROBEID) ||
      any(!nzchar(mapping$SYMBOL)) || any(!nzchar(mapping$PROBEID)) ||
      anyDuplicated(mapping$SYMBOL) || anyDuplicated(mapping$PROBEID)) {
    stop("Probe map identifiers must be nonempty and one-to-one")
  }
  if (any(!is.finite(mapping$training_probe_mean)) ||
      any(!is.finite(mapping$n_probes_gene)) ||
      any(mapping$n_probes_gene < 1L) ||
      any(mapping$n_probes_gene != as.integer(mapping$n_probes_gene))) {
    stop("Probe map contains invalid means or eligible-probe counts")
  }
  mapping$n_probes_gene <- as.integer(mapping$n_probes_gene)
  # `method = "radix"` is deliberately locale-independent. The production
  # script runs under Chinese_China.utf8 while the external validator may run
  # under C; default collation would give identical maps different hashes.
  mapping <- mapping[
    order(mapping$SYMBOL, mapping$PROBEID, method = "radix"),
    , drop = FALSE
  ]
  rownames(mapping) <- NULL
  mapping
}

hash_probe_map <- function(probe_map) {
  mapping <- .canonical_probe_mapping(probe_map)
  if (any(grepl("[[:cntrl:]]", mapping$SYMBOL)) ||
      any(grepl("[[:cntrl:]]", mapping$PROBEID))) {
    stop("Probe-map identifiers cannot contain tab or newline characters")
  }
  rows <- paste(
    enc2utf8(mapping$SYMBOL), enc2utf8(mapping$PROBEID),
    mapping$training_probe_mean_hex,
    sprintf("%d", mapping$n_probes_gene),
    sep = "\t"
  )
  payload <- paste(
    c(
      .probe_map_schema,
      "SYMBOL\tPROBEID\ttraining_probe_mean_hex\tn_probes_gene",
      rows
    ),
    collapse = "\n"
  )
  digest::digest(payload, algo = "sha256", serialize = FALSE)
}

learn_unique_highest_probe_map <- function(
    probe_expr, probe_annotation,
    training_sample_ids = colnames(probe_expr)) {
  x <- .validate_probe_expression(probe_expr, "probe-map training expression")
  training_sample_ids <- .validate_sample_set(
    training_sample_ids, colnames(x), "probe-map training sample IDs"
  )
  ann <- .validate_probe_annotation(probe_annotation)
  ann <- ann[
    !is.na(ann$PROBEID) & nzchar(ann$PROBEID) &
      ann$PROBEID %in% rownames(x) &
      !is.na(ann$SYMBOL) & nzchar(ann$SYMBOL),
    , drop = FALSE
  ]
  ann <- unique(ann)
  if (!nrow(ann)) {
    stop("No observed probe maps to a nonempty SYMBOL")
  }

  symbols_by_probe <- split(ann$SYMBOL, ann$PROBEID)
  symbols_by_probe <- lapply(symbols_by_probe, function(z) sort(unique(z)))
  eligible_ids <- rownames(x)[vapply(rownames(x), function(probe_id) {
    length(symbols_by_probe[[probe_id]]) == 1L
  }, logical(1))]
  if (!length(eligible_ids)) {
    stop("No observed probe maps to exactly one nonempty SYMBOL")
  }

  eligible <- data.frame(
    SYMBOL = vapply(eligible_ids, function(probe_id) {
      symbols_by_probe[[probe_id]][[1L]]
    }, character(1)),
    PROBEID = eligible_ids,
    training_probe_mean = rowMeans(x[eligible_ids, , drop = FALSE]),
    stringsAsFactors = FALSE
  )
  eligible$training_probe_mean_hex <- paste0(
    "hex:", sprintf("%a", eligible$training_probe_mean)
  )
  gene_counts <- table(eligible$SYMBOL)
  eligible$n_probes_gene <- as.integer(gene_counts[eligible$SYMBOL])
  eligible <- eligible[
    order(eligible$SYMBOL, -eligible$training_probe_mean, eligible$PROBEID),
    , drop = FALSE
  ]
  mapping <- eligible[!duplicated(eligible$SYMBOL), , drop = FALSE]
  mapping <- .canonical_probe_mapping(mapping)
  map_hash <- hash_probe_map(mapping)
  structure(
    list(
      mapping = mapping,
      map_hash = map_hash,
      training_sample_hash = .hash_probe_sample_ids(training_sample_ids),
      training_n = length(training_sample_ids),
      method = "unique_highest_mean",
      schema = .probe_map_schema
    ),
    class = c("crc_unique_highest_probe_map", "list")
  )
}

apply_frozen_probe_map <- function(probe_expr, probe_map) {
  x <- .validate_probe_expression(probe_expr, "frozen-map prediction expression")
  if (!is.list(probe_map) || is.null(probe_map$mapping) ||
      !identical(probe_map$method, "unique_highest_mean") ||
      !identical(probe_map$schema, .probe_map_schema) ||
      !is.character(probe_map$map_hash) || length(probe_map$map_hash) != 1L) {
    stop("Invalid frozen unique-highest-mean probe-map object")
  }
  mapping <- .canonical_probe_mapping(probe_map)
  recomputed_hash <- hash_probe_map(mapping)
  if (!identical(probe_map$map_hash, recomputed_hash)) {
    stop("Frozen probe-map hash mismatch")
  }
  missing <- setdiff(mapping$PROBEID, rownames(x))
  if (length(missing)) {
    stop("Prediction data missing frozen probes: ", paste(missing, collapse = ", "))
  }
  out <- x[match(mapping$PROBEID, rownames(x)), , drop = FALSE]
  rownames(out) <- mapping$SYMBOL
  attr(out, "probe_map_hash") <- recomputed_hash
  out
}

.combine_probe_model_fingerprint <- function(gene_model_fingerprint,
                                             probe_map_hash) {
  if (!is.character(gene_model_fingerprint) ||
      length(gene_model_fingerprint) != 1L ||
      !nzchar(gene_model_fingerprint) ||
      !is.character(probe_map_hash) || length(probe_map_hash) != 1L ||
      !nzchar(probe_map_hash)) {
    stop("Model and probe-map fingerprints must be nonempty strings")
  }
  digest::digest(
    list(
      schema = .probe_model_schema,
      gene_model_fingerprint = gene_model_fingerprint,
      probe_map_hash = probe_map_hash
    ),
    algo = "sha256", serialize = TRUE
  )
}

.hash_gene_signature_model <- function(model) {
  required <- c(
    "genes", "coefs", "center", "scale", "train_lp_mean", "train_lp_sd",
    "lambda_ratio", "lambda"
  )
  missing <- setdiff(required, names(model))
  if (length(missing)) {
    stop("Gene-signature model missing fields: ", paste(missing, collapse = ", "))
  }
  genes <- as.character(model$genes)
  if (!length(genes) || anyNA(genes) || any(!nzchar(genes)) ||
      anyDuplicated(genes) || any(!genes %in% names(model$coefs)) ||
      any(!genes %in% names(model$center)) || any(!genes %in% names(model$scale))) {
    stop("Gene-signature model has invalid feature metadata")
  }
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

.summarize_probe_frequency <- function(mapping_rows, denominator,
                                       analysis_key = NULL) {
  dt <- data.table::as.data.table(mapping_rows)
  if (!nrow(dt) || !all(c("SYMBOL", "PROBEID") %in% names(dt))) {
    stop("Cannot summarize an empty probe-mapping table")
  }
  out <- dt[, .(selected_splits = .N), by = .(SYMBOL, PROBEID)]
  out$denominator <- as.integer(denominator)
  out$frequency <- out$selected_splits / out$denominator
  if (!is.null(analysis_key)) out$analysis_key <- analysis_key
  data.table::setorder(out, SYMBOL, -frequency, PROBEID)
  as.data.frame(out)
}

fit_probe_signature_pipeline <- function(
    probe_expr, time, status, seed, params, probe_annotation) {
  validate_pipeline_params(params)
  checked <- .validate_probe_survival_inputs(
    probe_expr, time, status, "probe-signature training"
  )
  probe_expr <- checked$probe_expr
  time <- checked$time
  status <- checked$status
  if (any(table(factor(status, levels = 0:1)) < params$inner_folds)) {
    stop("Too few events or censored observations for the requested inner folds")
  }

  inner_foldid <- make_stratified_folds(status, params$inner_folds, seed)
  ratios <- sort(unique(params$lambda_ratio_grid), decreasing = TRUE)
  perf_rows <- vector("list", params$inner_folds * length(ratios))
  hash_rows <- vector("list", params$inner_folds)
  map_rows <- vector("list", params$inner_folds)
  counter <- 0L

  for (inner_fold in seq_len(params$inner_folds)) {
    inner_train <- which(inner_foldid != inner_fold)
    inner_valid <- which(inner_foldid == inner_fold)
    fold_map <- tryCatch(
      learn_unique_highest_probe_map(
        probe_expr[, inner_train, drop = FALSE], probe_annotation,
        training_sample_ids = colnames(probe_expr)[inner_train]
      ),
      error = function(e) {
        stop("Inner fold ", inner_fold, " probe-map learning failed: ",
             conditionMessage(e), call. = FALSE)
      }
    )
    train_expr <- apply_frozen_probe_map(
      probe_expr[, inner_train, drop = FALSE], fold_map
    )
    valid_expr <- apply_frozen_probe_map(
      probe_expr[, inner_valid, drop = FALSE], fold_map
    )
    prep <- tryCatch(
      prepare_training(train_expr, time[inner_train], status[inner_train], params),
      error = function(e) {
        stop("Inner fold ", inner_fold, " preprocessing failed: ",
             conditionMessage(e), call. = FALSE)
      }
    )
    path <- tryCatch(
      fit_lasso_path(prep, time[inner_train], status[inner_train], params),
      error = function(e) {
        stop("Inner fold ", inner_fold, " Cox LASSO path failed: ",
             conditionMessage(e), call. = FALSE)
      }
    )

    hash_rows[[inner_fold]] <- data.frame(
      inner_fold = inner_fold,
      probe_map_hash = fold_map$map_hash,
      training_sample_hash = fold_map$training_sample_hash,
      training_n = fold_map$training_n,
      stringsAsFactors = FALSE
    )
    map_rows[[inner_fold]] <- transform(
      fold_map$mapping, inner_fold = inner_fold
    )

    for (ratio in ratios) {
      counter <- counter + 1L
      candidate_model <- tryCatch(
        fit_at_ratio(
          prep, time[inner_train], status[inner_train], ratio, params,
          path = path
        ),
        error = function(e) {
          stop(
            "Inner fold ", inner_fold, " at lambda ratio ", signif(ratio, 6),
            " failed: ", conditionMessage(e), call. = FALSE
          )
        }
      )
      gene_fingerprint <- candidate_model$model_fingerprint
      combined_fingerprint <- .combine_probe_model_fingerprint(
        gene_fingerprint, fold_map$map_hash
      )
      pred <- predict_signature_pipeline(candidate_model, valid_expr)$lp_train_sd
      cidx <- unname(survival::concordance(
        survival::Surv(time[inner_valid], status[inner_valid]) ~ pred,
        reverse = TRUE
      )$concordance)
      if (!is.finite(cidx)) {
        stop("Non-finite inner-validation C-index in fold ", inner_fold,
             " at lambda ratio ", signif(ratio, 5))
      }
      perf_rows[[counter]] <- data.frame(
        inner_fold = inner_fold,
        lambda_ratio = ratio,
        c_index = cidx,
        n_train = length(inner_train),
        events_train = sum(status[inner_train]),
        n_valid = length(inner_valid),
        events_valid = sum(status[inner_valid]),
        n_genes = length(candidate_model$genes),
        probe_map_hash = fold_map$map_hash,
        gene_model_fingerprint = gene_fingerprint,
        model_fingerprint = combined_fingerprint,
        stringsAsFactors = FALSE
      )
    }
  }

  inner_perf <- data.table::rbindlist(perf_rows, fill = TRUE)
  tuning <- inner_perf[, .(
    mean_c_index = mean(c_index),
    sd_c_index = stats::sd(c_index),
    se_c_index = stats::sd(c_index) / sqrt(.N),
    min_c_index = min(c_index),
    max_c_index = max(c_index),
    inner_folds = .N
  ), by = lambda_ratio]
  data.table::setorder(tuning, -lambda_ratio)
  best_mean <- max(tuning$mean_c_index)
  best_candidates <- tuning[mean_c_index == best_mean]
  best_row <- best_candidates[which.max(lambda_ratio)]
  one_se_threshold <- best_row$mean_c_index - best_row$se_c_index
  eligible <- tuning[tuning$mean_c_index >= one_se_threshold]
  chosen_ratio <- max(eligible$lambda_ratio)
  tuning$best_mean <- tuning$mean_c_index == best_row$mean_c_index &
    tuning$lambda_ratio == best_row$lambda_ratio
  tuning$one_se_eligible <- tuning$mean_c_index >= one_se_threshold
  tuning$selected <- tuning$lambda_ratio == chosen_ratio
  tuning$one_se_threshold <- one_se_threshold

  final_map <- learn_unique_highest_probe_map(
    probe_expr, probe_annotation, training_sample_ids = colnames(probe_expr)
  )
  full_expr <- apply_frozen_probe_map(probe_expr, final_map)
  full_prep <- prepare_training(full_expr, time, status, params)
  full_path <- fit_lasso_path(full_prep, time, status, params)
  final <- fit_at_ratio(
    full_prep, time, status, chosen_ratio, params, path = full_path
  )
  final$gene_model_fingerprint <- final$model_fingerprint
  final$probe_map <- final_map
  final$probe_map_hash <- final_map$map_hash
  final$model_fingerprint <- .combine_probe_model_fingerprint(
    final$gene_model_fingerprint, final$probe_map_hash
  )
  final$lasso_path <- full_path
  final$candidate_genes <- full_prep$candidates
  final$variance_genes <- full_prep$variance_genes
  final$univ_df <- full_prep$univ_table
  final$tuning_curve <- as.data.frame(tuning)
  final$inner_fold_performance <- as.data.frame(inner_perf)
  final$inner_fold_assignments <- data.frame(
    sample = colnames(probe_expr),
    inner_fold = inner_foldid,
    stringsAsFactors = FALSE
  )
  final$inner_probe_map_hashes <- as.data.frame(
    data.table::rbindlist(hash_rows, fill = TRUE)
  )
  inner_mapping_rows <- data.table::rbindlist(map_rows, fill = TRUE)
  final$inner_probe_frequency <- .summarize_probe_frequency(
    inner_mapping_rows, params$inner_folds
  )
  final$params <- params
  final$seed <- as.integer(seed)
  final$probe_method <- "unique_highest_mean"
  final$audit <- data.frame(
    seed = as.integer(seed),
    input_probes = nrow(probe_expr),
    input_genes = nrow(full_expr),
    input_samples = ncol(probe_expr),
    input_events = sum(status),
    variance_retained = length(full_prep$variance_genes),
    univ_candidates = sum(
      is.finite(full_prep$univ_table$p) &
        full_prep$univ_table$p < params$univ_p_cutoff
    ),
    candidate_fallback = full_prep$candidate_fallback,
    candidates_entering_lasso = length(full_prep$candidates),
    selected_lambda_ratio = chosen_ratio,
    selected_lambda = final$lambda,
    lasso_nonzero = length(final$lasso_nonzero),
    selected_fallback = final$selected_fallback,
    final_genes = length(final$genes),
    inner_folds = params$inner_folds,
    probe_map_hash = final$probe_map_hash,
    stringsAsFactors = FALSE
  )
  final
}

predict_probe_signature_pipeline <- function(model, probe_expr_new) {
  required <- c(
    "probe_map", "probe_map_hash", "gene_model_fingerprint",
    "model_fingerprint"
  )
  missing <- setdiff(required, names(model))
  if (length(missing)) {
    stop("Probe-signature model missing fields: ", paste(missing, collapse = ", "))
  }
  map_hash <- hash_probe_map(model$probe_map)
  if (!identical(map_hash, model$probe_map_hash) ||
      !identical(map_hash, model$probe_map$map_hash)) {
    stop("Probe-signature model contains an inconsistent frozen map")
  }
  gene_hash <- .hash_gene_signature_model(model)
  if (!identical(gene_hash, model$gene_model_fingerprint)) {
    stop("Probe-signature gene-model fingerprint mismatch")
  }
  expected_model_hash <- .combine_probe_model_fingerprint(gene_hash, map_hash)
  if (!identical(expected_model_hash, model$model_fingerprint)) {
    stop("Probe-signature model fingerprint mismatch")
  }
  gene_expr <- apply_frozen_probe_map(probe_expr_new, model$probe_map)
  predict_signature_pipeline(model, gene_expr)
}

run_probe_nested_cv <- function(
    probe_expr, time, status, sample_ids, params, analysis_key,
    probe_annotation) {
  validate_pipeline_params(params)
  checked <- .validate_probe_survival_inputs(
    probe_expr, time, status, "probe nested CV"
  )
  probe_expr <- checked$probe_expr
  time <- checked$time
  status <- checked$status
  sample_ids <- .validate_sample_set(
    sample_ids, colnames(probe_expr), "nested-CV sample_ids"
  )
  if (!identical(sample_ids, colnames(probe_expr))) {
    stop("sample_ids must be in the same order as expression columns")
  }
  if (!is.character(analysis_key) || length(analysis_key) != 1L ||
      is.na(analysis_key) || !nzchar(analysis_key)) {
    stop("analysis_key must be one nonempty character value")
  }
  if (any(table(factor(status, levels = 0:1)) < params$outer_folds)) {
    stop("Too few events or censored observations for the requested outer folds")
  }

  fold_rows <- list()
  oof_rows <- list()
  selection_rows <- list()
  assignment_rows <- list()
  inner_assignment_rows <- list()
  tuning_rows <- list()
  inner_tuning_rows <- list()
  fingerprint_rows <- list()
  outer_map_rows <- list()
  inner_hash_rows <- list()
  inner_frequency_rows <- list()
  counter <- 0L

  for (repeat_id in seq_len(params$outer_repeats)) {
    outer_foldid <- make_stratified_folds(
      status, params$outer_folds, params$base_seed + repeat_id
    )
    assignment_rows[[repeat_id]] <- data.frame(
      sample = sample_ids,
      repeat_id = repeat_id,
      outer_fold = outer_foldid,
      analysis_key = analysis_key,
      stringsAsFactors = FALSE
    )
    for (outer_fold in seq_len(params$outer_folds)) {
      counter <- counter + 1L
      train_idx <- which(outer_foldid != outer_fold)
      test_idx <- which(outer_foldid == outer_fold)
      if (length(intersect(train_idx, test_idx)) ||
          !setequal(c(train_idx, test_idx), seq_len(ncol(probe_expr)))) {
        stop("Outer training/test partitions are not a disjoint exhaustive split")
      }
      if (!setequal(unique(status[test_idx]), 0:1) ||
          !setequal(unique(status[train_idx]), 0:1)) {
        stop("Outer fold ", outer_fold, " lacks events or censoring")
      }
      fold_seed <- params$base_seed + repeat_id * 1000L + outer_fold
      model <- tryCatch(
        fit_probe_signature_pipeline(
          probe_expr[, train_idx, drop = FALSE],
          time[train_idx], status[train_idx], seed = fold_seed,
          params = params, probe_annotation = probe_annotation
        ),
        error = function(e) {
          stop(
            "Probe nested CV failed at repeat ", repeat_id,
            ", outer fold ", outer_fold, ": ", conditionMessage(e),
            call. = FALSE
          )
        }
      )
      pred <- predict_probe_signature_pipeline(
        model, probe_expr[, test_idx, drop = FALSE]
      )
      cidx <- unname(survival::concordance(
        survival::Surv(time[test_idx], status[test_idx]) ~ pred$lp_train_sd,
        reverse = TRUE
      )$concordance)
      if (!is.finite(cidx)) {
        stop("Non-finite outer-fold C-index at repeat ", repeat_id,
             ", fold ", outer_fold)
      }

      fold_rows[[counter]] <- cbind(
        data.frame(
          repeat_id = repeat_id,
          outer_fold = outer_fold,
          n_train = length(train_idx),
          events_train = sum(status[train_idx]),
          n_test = length(test_idx),
          events_test = sum(status[test_idx]),
          c_index = cidx,
          analysis_key = analysis_key,
          model_fingerprint = model$model_fingerprint,
          stringsAsFactors = FALSE
        ),
        model$audit
      )
      oof_rows[[counter]] <- data.frame(
        sample = sample_ids[test_idx],
        repeat_id = repeat_id,
        outer_fold = outer_fold,
        time = time[test_idx],
        status = status[test_idx],
        lp_raw = pred$lp,
        lp_oof = pred$lp_train_sd,
        n_train = length(train_idx),
        events_train = sum(status[train_idx]),
        n_test = length(test_idx),
        events_test = sum(status[test_idx]),
        n_genes = length(model$genes),
        analysis_key = analysis_key,
        stringsAsFactors = FALSE
      )
      selection_rows[[counter]] <- data.frame(
        repeat_id = repeat_id,
        outer_fold = outer_fold,
        gene = model$genes,
        coefficient = unname(model$coefs[model$genes]),
        coefficient_sign = ifelse(
          model$coefs[model$genes] > 0, "positive",
          ifelse(model$coefs[model$genes] < 0, "negative", "zero")
        ),
        analysis_key = analysis_key,
        stringsAsFactors = FALSE
      )
      inner_assignment_rows[[counter]] <- transform(
        model$inner_fold_assignments,
        repeat_id = repeat_id,
        outer_fold = outer_fold,
        analysis_key = analysis_key
      )
      tuning_rows[[counter]] <- transform(
        model$tuning_curve,
        repeat_id = repeat_id,
        outer_fold = outer_fold,
        analysis_key = analysis_key
      )
      inner_tuning_rows[[counter]] <- transform(
        model$inner_fold_performance,
        repeat_id = repeat_id,
        outer_fold = outer_fold,
        analysis_key = analysis_key
      )
      fingerprint_rows[[counter]] <- data.frame(
        repeat_id = repeat_id,
        outer_fold = outer_fold,
        probe_map_hash = model$probe_map_hash,
        gene_model_fingerprint = model$gene_model_fingerprint,
        model_fingerprint = model$model_fingerprint,
        selected_lambda_ratio = model$lambda_ratio,
        analysis_key = analysis_key,
        stringsAsFactors = FALSE
      )
      outer_map_rows[[counter]] <- transform(
        model$probe_map$mapping,
        repeat_id = repeat_id,
        outer_fold = outer_fold,
        map_hash = model$probe_map_hash,
        training_sample_hash = model$probe_map$training_sample_hash,
        training_n = model$probe_map$training_n,
        analysis_key = analysis_key
      )
      inner_hash_rows[[counter]] <- transform(
        model$inner_probe_map_hashes,
        repeat_id = repeat_id,
        outer_fold = outer_fold,
        analysis_key = analysis_key
      )
      inner_frequency_rows[[counter]] <- transform(
        model$inner_probe_frequency,
        repeat_id = repeat_id,
        outer_fold = outer_fold,
        analysis_key = analysis_key
      )
    }
  }

  fold_performance <- data.table::rbindlist(fold_rows, fill = TRUE)
  oof <- data.table::rbindlist(oof_rows, fill = TRUE)
  selection <- data.table::rbindlist(selection_rows, fill = TRUE)
  outer_assignments <- data.table::rbindlist(assignment_rows, fill = TRUE)
  inner_assignments <- data.table::rbindlist(inner_assignment_rows, fill = TRUE)
  tuning_performance <- data.table::rbindlist(tuning_rows, fill = TRUE)
  inner_tuning_performance <- data.table::rbindlist(
    inner_tuning_rows, fill = TRUE
  )
  model_fingerprints <- data.table::rbindlist(fingerprint_rows, fill = TRUE)
  outer_probe_maps <- data.table::rbindlist(outer_map_rows, fill = TRUE)
  inner_probe_map_hashes <- data.table::rbindlist(inner_hash_rows, fill = TRUE)
  inner_probe_frequency_by_outer <- data.table::rbindlist(
    inner_frequency_rows, fill = TRUE
  )

  repeat_performance <- oof[, .(
    n = .N,
    events = sum(status),
    c_index = unname(survival::concordance(
      survival::Surv(time, status) ~ lp_oof, reverse = TRUE
    )$concordance)
  ), by = repeat_id]
  repeat_performance$analysis_key <- analysis_key
  total_models <- params$outer_repeats * params$outer_folds
  selection_frequency <- selection[, .(
    selected_folds = .N,
    frequency = .N / total_models,
    positive_folds = sum(coefficient > 0),
    negative_folds = sum(coefficient < 0),
    zero_folds = sum(coefficient == 0),
    modal_sign = if (sum(coefficient > 0) >= sum(coefficient < 0)) {
      "positive"
    } else {
      "negative"
    },
    sign_consistency = max(sum(coefficient > 0), sum(coefficient < 0)) / .N,
    coefficient_mean = mean(coefficient),
    coefficient_sd = if (.N > 1L) stats::sd(coefficient) else NA_real_,
    coefficient_median = stats::median(coefficient),
    coefficient_min = min(coefficient),
    coefficient_max = max(coefficient)
  ), by = gene]
  data.table::setorder(selection_frequency, -frequency, -sign_consistency, gene)
  selection_frequency$analysis_key <- analysis_key

  jaccard_rows <- list()
  j_counter <- 0L
  for (current_repeat in seq_len(params$outer_repeats)) {
    for (fold_a in seq_len(params$outer_folds - 1L)) {
      for (fold_b in seq.int(fold_a + 1L, params$outer_folds)) {
        j_counter <- j_counter + 1L
        genes_a <- selection[
          repeat_id == current_repeat & outer_fold == fold_a, gene
        ]
        genes_b <- selection[
          repeat_id == current_repeat & outer_fold == fold_b, gene
        ]
        jaccard_rows[[j_counter]] <- data.frame(
          repeat_id = current_repeat,
          fold_a = fold_a,
          fold_b = fold_b,
          intersection_n = length(intersect(genes_a, genes_b)),
          union_n = length(union(genes_a, genes_b)),
          jaccard = model_set_jaccard(genes_a, genes_b),
          analysis_key = analysis_key,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  jaccard <- data.table::rbindlist(jaccard_rows, fill = TRUE)

  outer_probe_frequency <- .summarize_probe_frequency(
    outer_probe_maps, total_models, analysis_key
  )
  inner_counts <- inner_probe_frequency_by_outer[, .(
    selected_splits = sum(selected_splits)
  ), by = .(SYMBOL, PROBEID)]
  inner_counts$denominator <- as.integer(total_models * params$inner_folds)
  inner_counts$frequency <- inner_counts$selected_splits / inner_counts$denominator
  inner_counts$analysis_key <- analysis_key
  data.table::setorder(inner_counts, SYMBOL, -frequency, PROBEID)
  inner_probe_frequency <- as.data.frame(inner_counts)

  full_seed <- as.integer(params$base_seed + 900000L)
  full_development_model <- fit_probe_signature_pipeline(
    probe_expr, time, status, seed = full_seed, params = params,
    probe_annotation = probe_annotation
  )

  summary <- data.frame(
    pipeline_version = as.character(params$pipeline_version),
    analysis_key = analysis_key,
    repeats = params$outer_repeats,
    outer_folds = params$outer_folds,
    inner_folds = params$inner_folds,
    completed_folds = nrow(fold_performance),
    median_repeat_oof_c_index = stats::median(repeat_performance$c_index),
    repeat_oof_c_index_q1 = unname(stats::quantile(
      repeat_performance$c_index, 0.25
    )),
    repeat_oof_c_index_q3 = unname(stats::quantile(
      repeat_performance$c_index, 0.75
    )),
    min_repeat_oof_c_index = min(repeat_performance$c_index),
    max_repeat_oof_c_index = max(repeat_performance$c_index),
    median_within_repeat_jaccard = stats::median(jaccard$jaccard),
    jaccard_q1 = unname(stats::quantile(jaccard$jaccard, 0.25)),
    jaccard_q3 = unname(stats::quantile(jaccard$jaccard, 0.75)),
    stringsAsFactors = FALSE
  )
  failure_log <- data.frame(
    repeat_id = integer(), outer_fold = integer(), stage = character(),
    message = character(), stringsAsFactors = FALSE
  )
  result <- list(
    pipeline_version = params$pipeline_version,
    analysis_key = analysis_key,
    params = params,
    probe_mapping_method = "unique_highest_mean",
    outer_fold_assignments = as.data.frame(outer_assignments),
    inner_fold_assignments = as.data.frame(inner_assignments),
    oof = as.data.frame(oof),
    fold_performance = as.data.frame(fold_performance),
    repeat_performance = as.data.frame(repeat_performance),
    selection = as.data.frame(selection),
    selection_frequency = as.data.frame(selection_frequency),
    jaccard = as.data.frame(jaccard),
    tuning_performance = as.data.frame(tuning_performance),
    inner_tuning_performance = as.data.frame(inner_tuning_performance),
    model_fingerprints = as.data.frame(model_fingerprints),
    outer_probe_maps = as.data.frame(outer_probe_maps),
    inner_probe_map_hashes = as.data.frame(inner_probe_map_hashes),
    outer_probe_frequency = as.data.frame(outer_probe_frequency),
    inner_probe_frequency_by_outer = as.data.frame(
      inner_probe_frequency_by_outer
    ),
    inner_probe_frequency = inner_probe_frequency,
    full_development_model = full_development_model,
    full_development_probe_map = full_development_model$probe_map$mapping,
    full_development_probe_map_hash = full_development_model$probe_map_hash,
    full_development_training_sample_hash =
      full_development_model$probe_map$training_sample_hash,
    full_development_seed = full_seed,
    summary = summary,
    failure_log = failure_log
  )
  validate_probe_nested_result(
    result, probe_expr, time, status, sample_ids, params, analysis_key,
    probe_annotation
  )
  result
}

.validate_probe_frequency <- function(tab, denominator, label,
                                      analysis_key = NULL) {
  tab <- as.data.frame(tab, stringsAsFactors = FALSE)
  required <- c(
    "SYMBOL", "PROBEID", "selected_splits", "denominator", "frequency"
  )
  missing <- setdiff(required, names(tab))
  if (length(missing) || !nrow(tab)) {
    stop(label, " is missing required frequency columns or rows")
  }
  if (anyNA(tab[, required]) || anyDuplicated(tab[c("SYMBOL", "PROBEID")]) ||
      any(tab$selected_splits < 1L) ||
      any(tab$selected_splits != as.integer(tab$selected_splits)) ||
      any(tab$denominator != as.integer(denominator)) ||
      any(abs(tab$frequency - tab$selected_splits / denominator) > 1e-12)) {
    stop(label, " contains invalid counts, denominators, or frequencies")
  }
  totals <- stats::aggregate(
    tab$selected_splits, list(SYMBOL = tab$SYMBOL), sum
  )
  if (any(totals$x != as.integer(denominator))) {
    stop(label, " does not conserve one probe choice per gene per split")
  }
  if (!is.null(analysis_key) &&
      (!"analysis_key" %in% names(tab) || any(tab$analysis_key != analysis_key))) {
    stop(label, " contains an incorrect analysis key")
  }
  invisible(TRUE)
}

.same_probe_frequency <- function(actual, expected, label) {
  cols <- c("SYMBOL", "PROBEID", "selected_splits", "denominator", "frequency")
  actual <- as.data.frame(actual, stringsAsFactors = FALSE)
  expected <- as.data.frame(expected, stringsAsFactors = FALSE)
  actual <- actual[, cols, drop = FALSE]
  expected <- expected[, cols, drop = FALSE]
  actual <- actual[order(actual$SYMBOL, actual$PROBEID), , drop = FALSE]
  expected <- expected[order(expected$SYMBOL, expected$PROBEID), , drop = FALSE]
  rownames(actual) <- NULL
  rownames(expected) <- NULL
  if (!isTRUE(all.equal(actual, expected, tolerance = 1e-12,
                        check.attributes = FALSE))) {
    stop(label, " does not match its auditable split-level counts")
  }
  invisible(TRUE)
}

validate_probe_nested_result <- function(
    result, probe_expr, time, status, sample_ids, params, analysis_key,
    probe_annotation) {
  checked <- .validate_probe_survival_inputs(
    probe_expr, time, status, "probe nested-result validation"
  )
  probe_expr <- checked$probe_expr
  time <- checked$time
  status <- checked$status
  sample_ids <- .validate_sample_set(
    sample_ids, colnames(probe_expr), "nested-result sample_ids"
  )
  if (!identical(sample_ids, colnames(probe_expr))) {
    stop("sample_ids must be in the same order as expression columns")
  }
  probe_annotation <- .validate_probe_annotation(probe_annotation)
  validate_nested_result(
    result, probe_expr, time, status, sample_ids, params, analysis_key
  )

  required <- c(
    "probe_mapping_method", "inner_tuning_performance", "outer_probe_maps",
    "inner_probe_map_hashes", "outer_probe_frequency",
    "inner_probe_frequency_by_outer", "inner_probe_frequency",
    "full_development_model", "full_development_probe_map",
    "full_development_probe_map_hash",
    "full_development_training_sample_hash", "full_development_seed"
  )
  missing <- setdiff(required, names(result))
  if (length(missing)) {
    stop("Probe nested-CV result missing fields: ", paste(missing, collapse = ", "))
  }
  if (!identical(result$probe_mapping_method, "unique_highest_mean")) {
    stop("Unexpected probe-mapping sensitivity method")
  }

  total_models <- params$outer_repeats * params$outer_folds
  outer_maps <- as.data.frame(result$outer_probe_maps, stringsAsFactors = FALSE)
  required_outer <- c(
    "repeat_id", "outer_fold", "SYMBOL", "PROBEID",
    "training_probe_mean", "training_probe_mean_hex", "n_probes_gene", "map_hash",
    "training_sample_hash", "training_n", "analysis_key"
  )
  if (!all(required_outer %in% names(outer_maps)) || !nrow(outer_maps) ||
      any(outer_maps$analysis_key != analysis_key)) {
    stop("Outer probe-map audit table is incomplete or incorrectly keyed")
  }
  group_key <- paste(outer_maps$repeat_id, outer_maps$outer_fold, sep = "\r")
  if (length(unique(group_key)) != total_models) {
    stop("Outer probe-map audit table has the wrong number of model groups")
  }

  assignments <- as.data.frame(result$outer_fold_assignments)
  fingerprints <- as.data.frame(result$model_fingerprints)
  fold_performance <- as.data.frame(result$fold_performance)
  for (repeat_id in seq_len(params$outer_repeats)) {
    repeat_assignment <- assignments[
      assignments$repeat_id == repeat_id, , drop = FALSE
    ]
    for (outer_fold in seq_len(params$outer_folds)) {
      rows <- outer_maps[
        outer_maps$repeat_id == repeat_id &
          outer_maps$outer_fold == outer_fold,
        , drop = FALSE
      ]
      if (!nrow(rows) || anyDuplicated(rows$SYMBOL) ||
          length(unique(rows$map_hash)) != 1L ||
          length(unique(rows$training_sample_hash)) != 1L ||
          length(unique(rows$training_n)) != 1L) {
        stop("Outer probe map is not unique within repeat/fold")
      }
      recomputed_map_hash <- hash_probe_map(rows)
      if (!identical(unique(rows$map_hash), recomputed_map_hash)) {
        stop("Outer probe-map hash mismatch")
      }
      test_samples <- repeat_assignment$sample[
        repeat_assignment$outer_fold == outer_fold
      ]
      train_samples <- setdiff(sample_ids, test_samples)
      train_idx <- match(train_samples, colnames(probe_expr))
      if (anyNA(train_idx)) {
        stop("Outer probe-map training samples are absent from source expression")
      }
      expected_outer_map <- learn_unique_highest_probe_map(
        probe_expr[, train_idx, drop = FALSE],
        probe_annotation,
        training_sample_ids = train_samples
      )
      if (!identical(
          .canonical_probe_mapping(rows),
          .canonical_probe_mapping(expected_outer_map$mapping))) {
        stop(
          "Outer probe map does not match source expression, annotation, ",
          "and training partition"
        )
      }
      if (unique(rows$training_n) != length(train_samples) ||
          !identical(
            unique(rows$training_sample_hash),
            .hash_probe_sample_ids(train_samples)
          ) || !identical(recomputed_map_hash, expected_outer_map$map_hash)) {
        stop("Outer probe map training-sample identity mismatch")
      }
      fp <- fingerprints[
        fingerprints$repeat_id == repeat_id &
          fingerprints$outer_fold == outer_fold,
        , drop = FALSE
      ]
      fold <- fold_performance[
        fold_performance$repeat_id == repeat_id &
          fold_performance$outer_fold == outer_fold,
        , drop = FALSE
      ]
      if (nrow(fp) != 1L || nrow(fold) != 1L ||
          fp$probe_map_hash != recomputed_map_hash ||
          fp$model_fingerprint != fold$model_fingerprint ||
          fp$model_fingerprint != .combine_probe_model_fingerprint(
            fp$gene_model_fingerprint, recomputed_map_hash
          )) {
        stop("Outer model/probe-map fingerprint linkage is inconsistent")
      }
    }
  }

  expected_outer_frequency <- .summarize_probe_frequency(
    outer_maps, total_models, analysis_key
  )
  .validate_probe_frequency(
    result$outer_probe_frequency, total_models,
    "Outer probe frequency", analysis_key
  )
  .same_probe_frequency(
    result$outer_probe_frequency, expected_outer_frequency,
    "Outer probe frequency"
  )

  inner_hashes <- as.data.frame(
    result$inner_probe_map_hashes, stringsAsFactors = FALSE
  )
  required_hashes <- c(
    "repeat_id", "outer_fold", "inner_fold", "probe_map_hash",
    "training_sample_hash", "training_n", "analysis_key"
  )
  expected_inner_splits <- total_models * params$inner_folds
  if (!all(required_hashes %in% names(inner_hashes)) ||
      nrow(inner_hashes) != expected_inner_splits ||
      anyDuplicated(inner_hashes[c("repeat_id", "outer_fold", "inner_fold")]) ||
      any(inner_hashes$analysis_key != analysis_key) ||
      any(!grepl("^[0-9a-f]{64}$", inner_hashes$probe_map_hash))) {
    stop("Inner probe-map hash table is incomplete or invalid")
  }
  inner_assignments <- as.data.frame(result$inner_fold_assignments)
  expected_inner_mapping_rows <- vector("list", nrow(inner_hashes))
  for (i in seq_len(nrow(inner_hashes))) {
    row <- inner_hashes[i, , drop = FALSE]
    inner <- inner_assignments[
      inner_assignments$repeat_id == row$repeat_id &
        inner_assignments$outer_fold == row$outer_fold,
      , drop = FALSE
    ]
    inner_train_samples <- inner$sample[inner$inner_fold != row$inner_fold]
    inner_train_idx <- match(inner_train_samples, colnames(probe_expr))
    if (anyNA(inner_train_idx)) {
      stop("Inner probe-map training samples are absent from source expression")
    }
    expected_inner_map <- learn_unique_highest_probe_map(
      probe_expr[, inner_train_idx, drop = FALSE],
      probe_annotation,
      training_sample_ids = inner_train_samples
    )
    if (row$training_n != length(inner_train_samples) ||
        row$training_sample_hash != .hash_probe_sample_ids(inner_train_samples) ||
        row$probe_map_hash != expected_inner_map$map_hash) {
      stop("Inner probe map training-sample identity mismatch")
    }
    expected_inner_mapping_rows[[i]] <- transform(
      expected_inner_map$mapping,
      repeat_id = row$repeat_id,
      outer_fold = row$outer_fold,
      inner_fold = row$inner_fold
    )
  }

  inner_tuning <- as.data.frame(
    result$inner_tuning_performance, stringsAsFactors = FALSE
  )
  required_tuning <- c(
    "repeat_id", "outer_fold", "inner_fold", "lambda_ratio",
    "probe_map_hash", "gene_model_fingerprint", "model_fingerprint",
    "analysis_key"
  )
  expected_tuning_rows <- expected_inner_splits *
    length(unique(params$lambda_ratio_grid))
  if (!all(required_tuning %in% names(inner_tuning)) ||
      nrow(inner_tuning) != expected_tuning_rows ||
      anyDuplicated(inner_tuning[
        c("repeat_id", "outer_fold", "inner_fold", "lambda_ratio")
      ]) || any(inner_tuning$analysis_key != analysis_key)) {
    stop("Inner tuning audit table is incomplete or invalid")
  }
  tuning_hash_key <- paste(
    inner_tuning$repeat_id, inner_tuning$outer_fold,
    inner_tuning$inner_fold, sep = "\r"
  )
  inner_hash_key <- paste(
    inner_hashes$repeat_id, inner_hashes$outer_fold,
    inner_hashes$inner_fold, sep = "\r"
  )
  matched_hash <- inner_hashes$probe_map_hash[match(tuning_hash_key, inner_hash_key)]
  expected_candidate_fp <- mapply(
    .combine_probe_model_fingerprint,
    inner_tuning$gene_model_fingerprint,
    inner_tuning$probe_map_hash,
    USE.NAMES = FALSE
  )
  if (anyNA(matched_hash) || any(inner_tuning$probe_map_hash != matched_hash) ||
      any(inner_tuning$model_fingerprint != expected_candidate_fp)) {
    stop("Inner tuning rows are not linked to their frozen probe maps")
  }

  by_outer <- as.data.frame(
    result$inner_probe_frequency_by_outer, stringsAsFactors = FALSE
  )
  required_by_outer <- c(
    "repeat_id", "outer_fold", "SYMBOL", "PROBEID", "selected_splits",
    "denominator", "frequency", "analysis_key"
  )
  if (!all(required_by_outer %in% names(by_outer)) || !nrow(by_outer) ||
      any(by_outer$analysis_key != analysis_key) ||
      anyDuplicated(by_outer[
        c("repeat_id", "outer_fold", "SYMBOL", "PROBEID")
      ])) {
    stop("Per-outer-model inner probe frequency table is invalid")
  }
  by_outer_groups <- split(
    by_outer, paste(by_outer$repeat_id, by_outer$outer_fold, sep = "\r")
  )
  if (length(by_outer_groups) != total_models) {
    stop("Per-outer-model inner probe frequencies omit a model")
  }
  for (tab in by_outer_groups) {
    .validate_probe_frequency(
      tab, params$inner_folds, "Per-outer-model inner probe frequency",
      analysis_key
    )
  }
  expected_inner_mappings <- data.table::rbindlist(
    expected_inner_mapping_rows, fill = TRUE
  )
  expected_by_outer <- expected_inner_mappings[, .(
    selected_splits = .N
  ), by = .(repeat_id, outer_fold, SYMBOL, PROBEID)]
  expected_by_outer$denominator <- as.integer(params$inner_folds)
  expected_by_outer$frequency <- expected_by_outer$selected_splits /
    expected_by_outer$denominator
  expected_by_outer$analysis_key <- analysis_key
  compare_by_outer_columns <- c(
    "repeat_id", "outer_fold", "SYMBOL", "PROBEID", "selected_splits",
    "denominator", "frequency", "analysis_key"
  )
  actual_by_outer <- by_outer[, compare_by_outer_columns, drop = FALSE]
  expected_by_outer <- as.data.frame(
    expected_by_outer[, ..compare_by_outer_columns]
  )
  actual_by_outer <- actual_by_outer[do.call(
    order, actual_by_outer[c("repeat_id", "outer_fold", "SYMBOL", "PROBEID")]
  ), , drop = FALSE]
  expected_by_outer <- expected_by_outer[do.call(
    order, expected_by_outer[c("repeat_id", "outer_fold", "SYMBOL", "PROBEID")]
  ), , drop = FALSE]
  rownames(actual_by_outer) <- rownames(expected_by_outer) <- NULL
  if (!isTRUE(all.equal(
      actual_by_outer, expected_by_outer,
      tolerance = 1e-12, check.attributes = FALSE))) {
    stop(
      "Per-outer-model inner probe frequencies do not match source-derived maps"
    )
  }
  inner_denominator <- expected_inner_splits
  .validate_probe_frequency(
    result$inner_probe_frequency, inner_denominator,
    "Inner probe frequency", analysis_key
  )
  aggregated_inner <- data.table::as.data.table(by_outer)[, .(
    selected_splits = sum(selected_splits)
  ), by = .(SYMBOL, PROBEID)]
  aggregated_inner$denominator <- as.integer(inner_denominator)
  aggregated_inner$frequency <- aggregated_inner$selected_splits /
    aggregated_inner$denominator
  aggregated_inner$analysis_key <- analysis_key
  .same_probe_frequency(
    result$inner_probe_frequency, aggregated_inner, "Inner probe frequency"
  )
  expected_inner_frequency <- .summarize_probe_frequency(
    expected_inner_mappings, inner_denominator, analysis_key
  )
  .same_probe_frequency(
    result$inner_probe_frequency, expected_inner_frequency,
    "Source-derived inner probe frequency"
  )

  full_map <- .canonical_probe_mapping(result$full_development_probe_map)
  full_map_hash <- hash_probe_map(full_map)
  expected_full_map <- learn_unique_highest_probe_map(
    probe_expr,
    probe_annotation,
    training_sample_ids = sample_ids
  )
  if (!identical(full_map, .canonical_probe_mapping(expected_full_map$mapping)) ||
      !identical(full_map_hash, expected_full_map$map_hash)) {
    stop(
      "Full-development probe map does not match source expression and annotation"
    )
  }
  if (!identical(result$full_development_probe_map_hash, full_map_hash) ||
      !identical(
        result$full_development_training_sample_hash,
        .hash_probe_sample_ids(sample_ids)
      )) {
    stop("Full-development probe-map metadata mismatch")
  }
  full_model <- result$full_development_model
  if (!is.list(full_model) ||
      !identical(.canonical_probe_mapping(full_model$probe_map), full_map) ||
      !identical(full_model$probe_map_hash, full_map_hash) ||
      !identical(full_model$probe_map$map_hash, full_map_hash) ||
      !identical(full_model$probe_map$training_sample_hash,
                 .hash_probe_sample_ids(sample_ids)) ||
      !identical(full_model$probe_map$training_n, length(sample_ids))) {
    stop("Full-development model does not carry the certified frozen map")
  }
  recomputed_gene_hash <- .hash_gene_signature_model(full_model)
  recomputed_model_hash <- .combine_probe_model_fingerprint(
    recomputed_gene_hash, full_map_hash
  )
  if (!identical(full_model$gene_model_fingerprint, recomputed_gene_hash) ||
      !identical(full_model$model_fingerprint, recomputed_model_hash)) {
    stop("Full-development model fingerprint mismatch")
  }
  full_inner_hashes <- as.data.frame(full_model$inner_probe_map_hashes)
  if (nrow(full_inner_hashes) != params$inner_folds ||
      anyDuplicated(full_inner_hashes$inner_fold)) {
    stop("Full-development inner probe-map hashes are incomplete")
  }
  expected_full_inner_mapping_rows <- vector(
    "list", nrow(full_inner_hashes)
  )
  for (i in seq_len(nrow(full_inner_hashes))) {
    fold_id <- full_inner_hashes$inner_fold[[i]]
    train_samples <- full_model$inner_fold_assignments$sample[
      full_model$inner_fold_assignments$inner_fold != fold_id
    ]
    train_idx <- match(train_samples, colnames(probe_expr))
    if (anyNA(train_idx)) {
      stop(
        "Full-development inner training samples are absent from source expression"
      )
    }
    expected_inner_map <- learn_unique_highest_probe_map(
      probe_expr[, train_idx, drop = FALSE],
      probe_annotation,
      training_sample_ids = train_samples
    )
    if (full_inner_hashes$training_n[[i]] != length(train_samples) ||
        full_inner_hashes$training_sample_hash[[i]] !=
          .hash_probe_sample_ids(train_samples) ||
        full_inner_hashes$probe_map_hash[[i]] != expected_inner_map$map_hash) {
      stop("Full-development inner probe-map training identity mismatch")
    }
    expected_full_inner_mapping_rows[[i]] <- transform(
      expected_inner_map$mapping, inner_fold = fold_id
    )
  }
  full_inner_tuning <- as.data.frame(full_model$inner_fold_performance)
  full_inner_hash_key <- as.character(full_inner_hashes$inner_fold)
  matched_full_inner_hash <- full_inner_hashes$probe_map_hash[match(
    as.character(full_inner_tuning$inner_fold), full_inner_hash_key
  )]
  expected_full_candidate_fp <- mapply(
    .combine_probe_model_fingerprint,
    full_inner_tuning$gene_model_fingerprint,
    full_inner_tuning$probe_map_hash,
    USE.NAMES = FALSE
  )
  if (anyNA(matched_full_inner_hash) ||
      any(full_inner_tuning$probe_map_hash != matched_full_inner_hash) ||
      any(full_inner_tuning$model_fingerprint != expected_full_candidate_fp)) {
    stop(
      "Full-development inner tuning rows are not linked to source-derived maps"
    )
  }
  .validate_probe_frequency(
    full_model$inner_probe_frequency, params$inner_folds,
    "Full-development inner probe frequency"
  )
  expected_full_inner_frequency <- .summarize_probe_frequency(
    data.table::rbindlist(expected_full_inner_mapping_rows, fill = TRUE),
    params$inner_folds
  )
  .same_probe_frequency(
    full_model$inner_probe_frequency,
    expected_full_inner_frequency,
    "Source-derived full-development inner probe frequency"
  )
  invisible(TRUE)
}
