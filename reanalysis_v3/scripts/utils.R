suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(timeROC)
  library(data.table)
  library(digest)
  library(ggplot2)
  library(AnnotationDbi)
  library(hgu133plus2.db)
  library(org.Hs.eg.db)
  library(metafor)
})

source_env <- Sys.getenv("CRC_SOURCE_ROOT", unset = "")
work_env <- Sys.getenv("CRC_WORK_ROOT", unset = "")
PKG_ROOT <- gsub("\\\\", "/", if (nzchar(source_env)) source_env else getwd())
WORK_ROOT <- gsub("\\\\", "/", if (nzchar(work_env)) work_env else getwd())
V3_ROOT <- file.path(WORK_ROOT, "reanalysis_v3")
SCRIPT_DIR <- file.path(V3_ROOT, "scripts")
INPUT_DIR <- file.path(PKG_ROOT, "data")
PREPARED_INPUT_DIR <- gsub(
  "\\\\", "/",
  Sys.getenv("CRC_PREPARED_INPUT_DIR", unset = file.path(V3_ROOT, "prepared_inputs"))
)
EXTERNAL_GEO_DIR <- file.path(INPUT_DIR, "external_GEO")
runs_env <- Sys.getenv("CRC_RUNS_ROOT", unset = "")
RUNS_ROOT <- gsub(
  "\\\\", "/",
  if (nzchar(runs_env)) runs_env else file.path(V3_ROOT, "runs")
)
CACHE_ROOT <- file.path(V3_ROOT, "cache")

SIGNATURE_GENES <- c(
  "CCDC134", "EIF4A2", "FGF19", "GJB6", "INHBB",
  "JAGN1", "LGALS9", "MSLN", "TAPBPL"
)

read_expr_csv <- function(path) {
  dt <- fread(path, check.names = FALSE)
  genes <- as.character(dt[[1]])
  samples <- names(dt)[-1L]
  if (anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes)) {
    stop("Expression input has missing, blank, or duplicated gene identifiers: ", path)
  }
  if (anyNA(samples) || any(!nzchar(samples)) || anyDuplicated(samples)) {
    stop("Expression input has missing, blank, or duplicated sample identifiers: ", path)
  }
  mat <- as.matrix(dt[, -1, with = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- genes
  if (any(!is.finite(mat))) stop("Expression input contains non-finite values: ", path)
  mat
}

validate_required_gene_ids <- function(required_genes) {
  required_genes <- as.character(required_genes)
  if (!length(required_genes) || anyNA(required_genes) ||
      any(!nzchar(required_genes)) || anyDuplicated(required_genes)) {
    stop("Required genes must be a non-empty, ordered vector of unique identifiers")
  }
  required_genes
}

# Build an outcome-blind, annotation-only crosswalk for current human SYMBOLs.
# Historical aliases are retained in the audit table even when they are not
# eligible. An alias is eligible only when its reverse lookup identifies one
# current SYMBOL and one ENTREZID, both identical to the required gene record.
build_human_symbol_crosswalk <- function(required_genes) {
  required_genes <- validate_required_gene_ids(required_genes)
  annotation_package <- "org.Hs.eg.db"
  annotation_version <- as.character(utils::packageVersion(annotation_package))

  forward <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db,
    keys = required_genes,
    columns = c("ENTREZID", "ALIAS"),
    keytype = "SYMBOL"
  ))
  forward$SYMBOL <- as.character(forward$SYMBOL)
  forward$ENTREZID <- as.character(forward$ENTREZID)
  forward$ALIAS <- as.character(forward$ALIAS)

  entrez_by_gene <- lapply(required_genes, function(gene) {
    values <- unique(forward$ENTREZID[
      forward$SYMBOL == gene & !is.na(forward$ENTREZID) &
        nzchar(forward$ENTREZID)
    ])
    if (length(values) != 1L) {
      stop(
        "Current SYMBOL does not map to exactly one ENTREZID in ",
        annotation_package, ": ", gene
      )
    }
    values[[1L]]
  })
  names(entrez_by_gene) <- required_genes

  aliases_by_gene <- lapply(required_genes, function(gene) {
    values <- unique(forward$ALIAS[
      forward$SYMBOL == gene & !is.na(forward$ALIAS) &
        nzchar(forward$ALIAS) & forward$ALIAS != gene
    ])
    sort(values)
  })
  names(aliases_by_gene) <- required_genes
  all_aliases <- sort(unique(unlist(aliases_by_gene, use.names = FALSE)))
  reverse <- if (length(all_aliases)) {
    suppressMessages(AnnotationDbi::select(
      org.Hs.eg.db,
      keys = all_aliases,
      columns = c("SYMBOL", "ENTREZID"),
      keytype = "ALIAS"
    ))
  } else {
    data.frame(
      ALIAS = character(), SYMBOL = character(), ENTREZID = character(),
      stringsAsFactors = FALSE
    )
  }
  reverse$ALIAS <- as.character(reverse$ALIAS)
  reverse$SYMBOL <- as.character(reverse$SYMBOL)
  reverse$ENTREZID <- as.character(reverse$ENTREZID)

  rows <- lapply(required_genes, function(gene) {
    required_entrez <- entrez_by_gene[[gene]]
    aliases <- aliases_by_gene[[gene]]
    if (!length(aliases)) {
      return(data.frame(
        required_gene = gene,
        required_entrez_id = required_entrez,
        alias = NA_character_,
        reverse_symbols = "",
        reverse_entrez_ids = "",
        reverse_symbol_count = 0L,
        reverse_entrez_count = 0L,
        reverse_unique = FALSE,
        symbol_consistent = FALSE,
        entrez_consistent = FALSE,
        eligible_alias = FALSE,
        annotation_package = annotation_package,
        annotation_version = annotation_version,
        stringsAsFactors = FALSE
      ))
    }
    do.call(rbind, lapply(aliases, function(alias) {
      alias_rows <- reverse[reverse$ALIAS == alias, , drop = FALSE]
      reverse_symbols <- sort(unique(alias_rows$SYMBOL[
        !is.na(alias_rows$SYMBOL) & nzchar(alias_rows$SYMBOL)
      ]))
      reverse_entrez <- sort(unique(alias_rows$ENTREZID[
        !is.na(alias_rows$ENTREZID) & nzchar(alias_rows$ENTREZID)
      ]))
      reverse_unique <- length(reverse_symbols) == 1L &&
        length(reverse_entrez) == 1L
      symbol_consistent <- reverse_unique &&
        identical(reverse_symbols[[1L]], gene)
      entrez_consistent <- reverse_unique &&
        identical(reverse_entrez[[1L]], required_entrez)
      data.frame(
        required_gene = gene,
        required_entrez_id = required_entrez,
        alias = alias,
        reverse_symbols = paste(reverse_symbols, collapse = ";"),
        reverse_entrez_ids = paste(reverse_entrez, collapse = ";"),
        reverse_symbol_count = length(reverse_symbols),
        reverse_entrez_count = length(reverse_entrez),
        reverse_unique = reverse_unique,
        symbol_consistent = symbol_consistent,
        entrez_consistent = entrez_consistent,
        eligible_alias = reverse_unique && symbol_consistent &&
          entrez_consistent,
        annotation_package = annotation_package,
        annotation_version = annotation_version,
        stringsAsFactors = FALSE
      )
    }))
  })
  crosswalk <- do.call(rbind, rows)
  crosswalk$required_gene <- factor(
    crosswalk$required_gene, levels = required_genes, ordered = TRUE
  )
  crosswalk <- crosswalk[
    order(crosswalk$required_gene, crosswalk$alias, na.last = TRUE),
    , drop = FALSE
  ]
  crosswalk$required_gene <- as.character(crosswalk$required_gene)
  rownames(crosswalk) <- NULL
  crosswalk
}

# Resolve expression rows without consulting expression values or outcomes.
# Exact current SYMBOLs always take precedence. If absent, the only permitted
# fallback is exactly one expression-row alias that is reverse-unique and
# ENTREZ-consistent in the supplied frozen annotation crosswalk.
resolve_expression_gene_rows <- function(
    expr, required_genes, crosswalk = NULL, cohort = "expression cohort") {
  required_genes <- validate_required_gene_ids(required_genes)
  expr <- as.matrix(expr)
  storage.mode(expr) <- "double"
  row_ids <- rownames(expr)
  if (is.null(row_ids) || anyNA(row_ids) || any(!nzchar(row_ids)) ||
      anyDuplicated(row_ids)) {
    stop(cohort, " expression has invalid or duplicated gene identifiers")
  }
  if (is.null(colnames(expr)) || anyNA(colnames(expr)) ||
      any(!nzchar(colnames(expr))) || anyDuplicated(colnames(expr))) {
    stop(cohort, " expression has invalid or duplicated sample identifiers")
  }
  if (any(!is.finite(expr))) stop(cohort, " expression contains non-finite values")
  if (is.null(crosswalk)) {
    crosswalk <- build_human_symbol_crosswalk(required_genes)
  }
  crosswalk <- as.data.frame(crosswalk, stringsAsFactors = FALSE)
  required_columns <- c(
    "required_gene", "required_entrez_id", "alias", "reverse_symbols",
    "reverse_entrez_ids", "reverse_symbol_count", "reverse_entrez_count",
    "reverse_unique", "symbol_consistent", "entrez_consistent",
    "eligible_alias", "annotation_package", "annotation_version"
  )
  if (!all(required_columns %in% names(crosswalk))) {
    stop("Human SYMBOL crosswalk lacks required audit columns")
  }
  crosswalk$required_gene <- as.character(crosswalk$required_gene)
  crosswalk$required_entrez_id <- as.character(crosswalk$required_entrez_id)
  crosswalk$alias <- as.character(crosswalk$alias)
  for (field in c(
      "reverse_unique", "symbol_consistent", "entrez_consistent",
      "eligible_alias")) {
    crosswalk[[field]] <- as.logical(crosswalk[[field]])
    if (anyNA(crosswalk[[field]])) stop("Crosswalk has missing ", field)
  }
  if (anyNA(crosswalk$required_gene) ||
      any(!nzchar(crosswalk$required_gene)) ||
      anyNA(crosswalk$required_entrez_id) ||
      any(!nzchar(crosswalk$required_entrez_id))) {
    stop("Crosswalk has invalid required-gene or ENTREZ identifiers")
  }
  missing_crosswalk <- setdiff(required_genes, unique(crosswalk$required_gene))
  if (length(missing_crosswalk)) {
    stop("Crosswalk is missing required genes: ",
         paste(missing_crosswalk, collapse = ", "))
  }

  source_rows <- character(length(required_genes))
  audit_rows <- vector("list", length(required_genes))
  for (i in seq_along(required_genes)) {
    gene <- required_genes[[i]]
    gene_crosswalk <- crosswalk[
      crosswalk$required_gene == gene, , drop = FALSE
    ]
    required_entrez <- unique(gene_crosswalk$required_entrez_id)
    if (length(required_entrez) != 1L) {
      stop("Crosswalk does not define one required ENTREZID for ", gene)
    }
    alias_rows <- gene_crosswalk[
      !is.na(gene_crosswalk$alias) & nzchar(gene_crosswalk$alias),
      , drop = FALSE
    ]
    present_alias_rows <- alias_rows[
      alias_rows$alias %in% row_ids, , drop = FALSE
    ]
    exact_present <- gene %in% row_ids
    if (exact_present) {
      source_row <- gene
      method <- "exact_SYMBOL"
      reverse_unique <- NA
      symbol_consistent <- NA
      entrez_consistent <- TRUE
      reverse_symbols <- ""
      reverse_entrez_ids <- ""
      reverse_symbol_count <- NA_integer_
      reverse_entrez_count <- NA_integer_
    } else {
      if (!nrow(present_alias_rows)) {
        stop(
          cohort, " missing required gene ", gene,
          ": no exact SYMBOL or annotated historical alias row"
        )
      }
      invalid_present <- present_alias_rows[!present_alias_rows$eligible_alias, , drop = FALSE]
      if (nrow(invalid_present)) {
        stop(
          cohort, " has ambiguous or ENTREZ-inconsistent alias candidate(s) for ",
          gene, ": ", paste(sort(unique(invalid_present$alias)), collapse = ", ")
        )
      }
      eligible_present <- present_alias_rows[present_alias_rows$eligible_alias, , drop = FALSE]
      eligible_aliases <- sort(unique(eligible_present$alias))
      if (length(eligible_aliases) != 1L) {
        stop(
          cohort, " requires exactly one eligible historical alias for ", gene,
          "; found ", length(eligible_aliases), ": ",
          paste(eligible_aliases, collapse = ", ")
        )
      }
      source_row <- eligible_aliases[[1L]]
      selected <- eligible_present[
        eligible_present$alias == source_row, , drop = FALSE
      ]
      if (nrow(selected) != 1L) {
        stop("Crosswalk contains duplicated eligible alias records for ", gene)
      }
      method <- "unique_reverse_unique_ENTREZ_alias"
      reverse_unique <- selected$reverse_unique[[1L]]
      symbol_consistent <- selected$symbol_consistent[[1L]]
      entrez_consistent <- selected$entrez_consistent[[1L]]
      reverse_symbols <- selected$reverse_symbols[[1L]]
      reverse_entrez_ids <- selected$reverse_entrez_ids[[1L]]
      reverse_symbol_count <- selected$reverse_symbol_count[[1L]]
      reverse_entrez_count <- selected$reverse_entrez_count[[1L]]
    }
    source_rows[[i]] <- source_row
    audit_rows[[i]] <- data.frame(
      required_gene = gene,
      required_entrez_id = required_entrez,
      source_row_id = source_row,
      resolution_method = method,
      exact_symbol_present = exact_present,
      expression_alias_candidates = paste(
        sort(unique(present_alias_rows$alias)), collapse = ";"
      ),
      expression_alias_candidate_count = length(unique(present_alias_rows$alias)),
      reverse_unique = reverse_unique,
      symbol_consistent = symbol_consistent,
      entrez_consistent = entrez_consistent,
      reverse_symbols = reverse_symbols,
      reverse_entrez_ids = reverse_entrez_ids,
      reverse_symbol_count = reverse_symbol_count,
      reverse_entrez_count = reverse_entrez_count,
      annotation_package = unique(gene_crosswalk$annotation_package)[[1L]],
      annotation_version = unique(gene_crosswalk$annotation_version)[[1L]],
      outcome_or_expression_values_used_for_resolution = FALSE,
      stringsAsFactors = FALSE
    )
  }
  if (anyDuplicated(source_rows)) {
    stop(cohort, " gene resolution reuses one expression row for multiple genes")
  }
  resolved <- expr[source_rows, , drop = FALSE]
  rownames(resolved) <- required_genes
  audit <- do.call(rbind, audit_rows)
  rownames(audit) <- NULL
  if (!identical(rownames(resolved), required_genes) ||
      !identical(audit$required_gene, required_genes)) {
    stop("Resolved expression rows do not preserve required-gene order")
  }
  list(expr = resolved, audit = audit, crosswalk = crosswalk)
}

zscore_rows <- function(mat, center = NULL, scale = NULL) {
  mat <- as.matrix(mat)
  if (is.null(center)) center <- rowMeans(mat, na.rm = TRUE)
  if (is.null(scale)) scale <- apply(mat, 1, sd, na.rm = TRUE)
  scale[!is.finite(scale) | scale == 0] <- 1
  z <- sweep(sweep(mat, 1, center, "-"), 1, scale, "/")
  list(z = z, center = center, scale = scale)
}

read_geo_matrix <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  skip <- 0L
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) stop("GEO matrix begin marker not found: ", path)
    skip <- skip + 1L
    if (identical(line, "!series_matrix_table_begin")) break
  }
  close(con)
  on.exit(NULL, add = FALSE)
  tab <- read.delim(
    gzfile(path), skip = skip, header = TRUE, check.names = FALSE,
    comment.char = "!", quote = "\"", stringsAsFactors = FALSE
  )
  probe <- as.character(tab[[1]])
  mat <- as.matrix(tab[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- probe
  colnames(mat) <- sub('^"|"$', "", colnames(mat))
  mat
}

map_probes <- function(
    probe_mat,
    method = c("unique_highest_mean", "unique_mean"),
    required_genes = NULL) {
  method <- match.arg(method)
  probe_mat <- as.matrix(probe_mat)
  storage.mode(probe_mat) <- "double"
  probe_ids <- rownames(probe_mat)
  if (is.null(probe_ids) || anyNA(probe_ids) || any(!nzchar(probe_ids)) ||
      anyDuplicated(probe_ids)) {
    stop("Probe matrix row names must be nonempty, nonmissing, and unique")
  }
  if (is.null(colnames(probe_mat)) || anyNA(colnames(probe_mat)) ||
      any(!nzchar(colnames(probe_mat))) || anyDuplicated(colnames(probe_mat))) {
    stop("Probe matrix sample names must be nonempty, nonmissing, and unique")
  }
  if (any(!is.finite(probe_mat))) stop("Probe matrix contains non-finite values")

  ann <- suppressMessages(AnnotationDbi::select(
    hgu133plus2.db,
    keys = probe_ids,
    columns = "SYMBOL",
    keytype = "PROBEID"
  ))
  ann <- ann[
    !is.na(ann$SYMBOL) & nzchar(ann$SYMBOL) & ann$PROBEID %in% probe_ids,
    c("PROBEID", "SYMBOL"),
    drop = FALSE
  ]
  ann$PROBEID <- as.character(ann$PROBEID)
  ann$SYMBOL <- as.character(ann$SYMBOL)
  ann <- unique(ann)

  symbols_by_probe <- split(ann$SYMBOL, ann$PROBEID)
  symbols_by_probe <- lapply(symbols_by_probe, function(x) sort(unique(x)))
  n_symbols <- lengths(symbols_by_probe)[probe_ids]
  n_symbols[is.na(n_symbols)] <- 0L
  symbol_text <- vapply(probe_ids, function(id) {
    z <- symbols_by_probe[[id]]
    if (is.null(z)) "" else paste(z, collapse = ";")
  }, character(1))

  probe_audit <- data.frame(
    PROBEID = probe_ids,
    symbols = unname(symbol_text),
    n_symbols_per_probe = as.integer(n_symbols),
    eligible = as.integer(n_symbols) == 1L,
    stringsAsFactors = FALSE
  )
  probe_audit$SYMBOL <- ifelse(
    probe_audit$eligible,
    probe_audit$symbols,
    NA_character_
  )
  probe_audit$probe_mean <- rowMeans(probe_mat, na.rm = TRUE)

  eligible <- probe_audit[probe_audit$eligible, , drop = FALSE]
  if (!nrow(eligible)) stop("No observed probe maps to exactly one nonempty SYMBOL")
  eligible$idx <- match(eligible$PROBEID, probe_ids)
  gene_probe_counts <- table(eligible$SYMBOL)
  eligible$n_probes_gene <- as.integer(gene_probe_counts[eligible$SYMBOL])
  probe_audit$n_probes_gene <- as.integer(gene_probe_counts[probe_audit$SYMBOL])

  if (method == "unique_highest_mean") {
    eligible <- eligible[
      order(eligible$SYMBOL, -eligible$probe_mean, eligible$PROBEID),
      ,
      drop = FALSE
    ]
    chosen <- eligible[!duplicated(eligible$SYMBOL), , drop = FALSE]
    out <- probe_mat[chosen$idx, , drop = FALSE]
    rownames(out) <- chosen$SYMBOL
    mapping <- data.frame(
      SYMBOL = chosen$SYMBOL,
      method = method,
      n_symbols_per_probe = chosen$n_symbols_per_probe,
      eligible = TRUE,
      n_probes_gene = chosen$n_probes_gene,
      actual_probe = chosen$PROBEID,
      PROBEID = chosen$PROBEID,
      probe_mean = chosen$probe_mean,
      stringsAsFactors = FALSE
    )
    used_probe_ids <- chosen$PROBEID
  } else {
    symbols <- sort(unique(eligible$SYMBOL))
    out <- vapply(symbols, function(sym) {
      idx <- eligible$idx[eligible$SYMBOL == sym]
      if (length(idx) == 1L) {
        probe_mat[idx, ]
      } else {
        colMeans(probe_mat[idx, , drop = FALSE])
      }
    }, numeric(ncol(probe_mat)))
    out <- t(out)
    colnames(out) <- colnames(probe_mat)
    rownames(out) <- symbols
    actual_probes <- vapply(symbols, function(sym) {
      paste(sort(eligible$PROBEID[eligible$SYMBOL == sym]), collapse = ";")
    }, character(1))
    mapping <- data.frame(
      SYMBOL = symbols,
      method = method,
      n_symbols_per_probe = 1L,
      eligible = TRUE,
      n_probes_gene = as.integer(gene_probe_counts[symbols]),
      actual_probe = actual_probes,
      PROBEID = actual_probes,
      probe_mean = vapply(symbols, function(sym) {
        mean(eligible$probe_mean[eligible$SYMBOL == sym])
      }, numeric(1)),
      stringsAsFactors = FALSE
    )
    used_probe_ids <- eligible$PROBEID
  }

  missing_required <- setdiff(unique(as.character(required_genes)), rownames(out))
  if (length(missing_required)) {
    stop(
      "Required signature genes missing after ", method,
      " one-to-one probe filtering: ", paste(missing_required, collapse = ", ")
    )
  }

  probe_audit$method <- method
  probe_audit$used_in_expression <- probe_audit$PROBEID %in% used_probe_ids
  probe_audit$actual_probe <- ifelse(
    probe_audit$used_in_expression,
    probe_audit$PROBEID,
    NA_character_
  )
  probe_audit$selection_reason <- if (method == "unique_mean") {
    ifelse(
      probe_audit$n_symbols_per_probe == 0L,
      "ineligible_no_SYMBOL",
      ifelse(
        probe_audit$n_symbols_per_probe > 1L,
        "ineligible_multiple_SYMBOLs",
        "mean_all_eligible"
      )
    )
  } else {
    ifelse(
      probe_audit$n_symbols_per_probe == 0L,
      "ineligible_no_SYMBOL",
      ifelse(
        probe_audit$n_symbols_per_probe > 1L,
        "ineligible_multiple_SYMBOLs",
        ifelse(
          probe_audit$used_in_expression,
          "highest_mean_then_PROBEID",
          "eligible_not_highest_mean"
        )
      )
    )
  }
  probe_audit <- probe_audit[, c(
    "PROBEID", "symbols", "n_symbols_per_probe", "eligible", "SYMBOL",
    "n_probes_gene", "probe_mean", "method", "used_in_expression",
    "actual_probe", "selection_reason"
  )]

  list(expr = out, mapping = mapping, probe_audit = probe_audit)
}

read_geo_metadata <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con))
  header <- character()
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line) || identical(line, "!series_matrix_table_begin")) break
    header <- c(header, line)
  }
  title_line <- grep("^!Sample_geo_accession", header, value = TRUE)
  if (!length(title_line)) stop("Sample accessions missing: ", path)
  samples <- strsplit(title_line[[1]], "\t", fixed = TRUE)[[1]][-1]
  samples <- gsub('^"|"$', "", samples)
  char_lines <- grep("^!Sample_characteristics_ch1", header, value = TRUE)
  chars <- lapply(char_lines, function(x) {
    v <- strsplit(x, "\t", fixed = TRUE)[[1]][-1]
    v <- gsub('^"|"$', "", v)
    length(v) <- length(samples)
    v
  })
  per_sample <- lapply(seq_along(samples), function(j) {
    vals <- vapply(chars, function(x) x[[j]], character(1))
    vals[!is.na(vals) & nzchar(vals)]
  })
  names(per_sample) <- samples
  list(samples = samples, fields = per_sample)
}

extract_numeric <- function(fields, pattern) {
  hit <- grep(pattern, fields, ignore.case = TRUE, value = TRUE)
  if (!length(hit)) return(NA_real_)
  x <- trimws(sub("^[^:]*:[[:space:]]*", "", hit[[1]]))
  suppressWarnings(as.numeric(sub("[^0-9.].*$", "", x)))
}

extract_text <- function(fields, pattern) {
  hit <- grep(pattern, fields, ignore.case = TRUE, value = TRUE)
  if (!length(hit)) return(NA_character_)
  trimws(sub("^[^:]*:[[:space:]]*", "", hit[[1]]))
}

parse_gse14333_clin <- function(meta) {
  rows <- lapply(meta$samples, function(s) {
    fields <- meta$fields[[s]]
    joined <- paste(fields, collapse = "; ")
    get_num <- function(key) {
      m <- regexec(paste0(key, ":[ ]*([0-9.]+)"), joined, ignore.case = TRUE)
      z <- regmatches(joined, m)[[1]]
      if (length(z) < 2L) NA_real_ else as.numeric(z[[2]])
    }
    get_txt <- function(key, value_pattern = "[^;]+") {
      m <- regexec(paste0(key, ":[ ]*", "(", value_pattern, ")"), joined, ignore.case = TRUE)
      z <- regmatches(joined, m)[[1]]
      if (length(z) < 2L) NA_character_ else trimws(z[[2]])
    }
    cens <- get_num("DFS_Cens")
    data.frame(
      sample = s,
      time = get_num("DFS_Time") * 30.4375,
      status = ifelse(is.na(cens), NA_real_, 1 - cens),
      stage = get_txt("DukesStage", "[A-D]"),
      age = get_num("Age_Diag"),
      endpoint = "DFS",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

parse_gse17x_clin <- function(meta) {
  rows <- lapply(meta$samples, function(s) {
    fields <- meta$fields[[s]]
    event_text <- extract_text(fields, "^overall_event")
    status <- if (is.na(event_text)) NA_real_ else if (grepl("no death", event_text, ignore.case = TRUE)) 0 else if (grepl("death", event_text, ignore.case = TRUE)) 1 else NA_real_
    stage_text <- extract_text(fields, "^ajcc_stage|^stage")
    age <- extract_numeric(fields, "^age")
    data.frame(
      sample = s,
      time = extract_numeric(fields, "^overall survival follow-up time") * 30.4375,
      status = status,
      stage = stage_text,
      age = age,
      endpoint = "OS",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

apply_signature <- function(expr, clin, coefs, cohort, mapping_method) {
  expr <- as.matrix(expr)
  storage.mode(expr) <- "double"
  if (is.null(rownames(expr)) || anyNA(rownames(expr)) ||
      any(!nzchar(rownames(expr))) || anyDuplicated(rownames(expr))) {
    stop(cohort, " expression has invalid or duplicated gene identifiers")
  }
  if (is.null(colnames(expr)) || anyNA(colnames(expr)) ||
      any(!nzchar(colnames(expr))) || anyDuplicated(colnames(expr))) {
    stop(cohort, " expression has invalid or duplicated sample identifiers")
  }
  if (any(!is.finite(expr))) stop(cohort, " expression contains non-finite values")
  if (!all(c("sample", "time", "status", "endpoint") %in% names(clin))) {
    stop(cohort, " clinical table lacks required columns")
  }
  if (anyNA(clin$sample) || any(!nzchar(as.character(clin$sample))) ||
      anyDuplicated(as.character(clin$sample))) {
    stop(cohort, " clinical table has invalid or duplicated sample identifiers")
  }
  if (is.null(names(coefs)) || anyNA(names(coefs)) ||
      any(!nzchar(names(coefs))) || anyDuplicated(names(coefs)) ||
      any(!is.finite(coefs))) {
    stop(cohort, " signature coefficients are invalid")
  }
  missing_genes <- setdiff(names(coefs), rownames(expr))
  if (length(missing_genes)) {
    stop(cohort, " missing signature genes: ", paste(missing_genes, collapse = ", "))
  }
  genes <- names(coefs)
  clin <- clin[clin$sample %in% colnames(expr) & is.finite(clin$time) & clin$time > 0 & clin$status %in% 0:1, , drop = FALSE]
  if (!nrow(clin) || length(unique(clin$status)) != 2L ||
      length(unique(clin$endpoint)) != 1L) {
    stop(cohort, " has no valid two-outcome survival analysis set")
  }
  expr <- expr[genes, clin$sample, drop = FALSE]
  if (!identical(colnames(expr), as.character(clin$sample))) {
    stop(cohort, " expression and clinical sample orders differ")
  }
  z <- zscore_rows(expr)$z
  risk <- as.numeric(t(z) %*% coefs[genes])
  risk_scale <- sd(risk)
  if (any(!is.finite(risk)) || !is.finite(risk_scale) || risk_scale <= 0) {
    stop(cohort, " produced a non-finite or zero-variance risk score")
  }
  risk_sd <- as.numeric((risk - mean(risk)) / risk_scale)
  data.frame(
    clin,
    cohort = cohort,
    mapping_method = mapping_method,
    preprocessing_scope = "cohort_adaptive_standardization",
    estimand = "within_cohort_prognostic_association",
    risk = risk,
    risk_sd = risk_sd,
    group = factor(ifelse(risk > median(risk, na.rm = TRUE), "High", "Low"), levels = c("Low", "High")),
    stringsAsFactors = FALSE
  )
}

cohort_stats <- function(d) {
  required <- c("time", "status", "risk_sd", "group", "cohort", "endpoint", "mapping_method")
  d <- as.data.frame(d)
  if (!all(required %in% names(d)) || any(!complete.cases(d[, required, drop = FALSE]))) {
    stop("Cohort performance input is incomplete")
  }
  fit <- coxph(
    Surv(time, status) ~ risk_sd,
    data = d,
    ties = "efron",
    x = TRUE,
    y = TRUE,
    na.action = na.fail
  )
  if (fit$n != nrow(d) || fit$nevent != sum(d$status)) {
    stop("Cox model frame does not match the reported cohort denominator")
  }
  s <- summary(fit)
  km <- survdiff(Surv(time, status) ~ group, data = d)
  c_obj <- concordance(
    Surv(time, status) ~ risk_sd,
    data = d,
    reverse = TRUE,
    timewt = "n"
  )
  out <- data.frame(
    cohort = unique(d$cohort),
    endpoint = unique(d$endpoint),
    mapping_method = unique(d$mapping_method),
    n = nrow(d),
    events = sum(d$status),
    HR_per_SD = unname(exp(coef(fit))),
    lower95 = unname(exp(confint(fit)[1])),
    upper95 = unname(exp(confint(fit)[2])),
    cox_p = s$coefficients[1, "Pr(>|z|)"],
    logrank_p = pchisq(km$chisq, df = 1, lower.tail = FALSE),
    c_index = unname(c_obj$concordance),
    c_index_method = "Harrell C; timewt=n; higher score=higher hazard",
    stringsAsFactors = FALSE
  )
  attr(out, "fit") <- fit
  out
}

eligible_auc <- function(d, years = c(1, 3, 5)) {
  min_events <- 10L
  min_at_risk <- 20L
  rows <- lapply(years, function(y) {
    tt <- y * 365.25
    events_before <- sum(d$status == 1 & d$time <= tt)
    at_risk <- sum(d$time >= tt)
    if (events_before < min_events || at_risk < min_at_risk || tt >= max(d$time)) {
      reason <- if (tt >= max(d$time)) "time_beyond_followup" else "insufficient_events_or_at_risk"
      return(data.frame(
        year = y, AUC = NA_real_, lower95 = NA_real_, upper95 = NA_real_,
        eligible = FALSE, auc_estimable = FALSE, ci_estimable = FALSE,
        reason = reason, failure_class = "eligibility",
        error_message = "", events_before = events_before, at_risk = at_risk,
        min_events_required = min_events, min_at_risk_required = min_at_risk,
        weighting = "marginal"
      ))
    }
    roc_error <- ""
    roc <- tryCatch(
      timeROC(
        T = d$time, delta = d$status, marker = d$risk_sd,
        cause = 1, times = tt, iid = TRUE, weighting = "marginal"
      ),
      error = function(e) {
        roc_error <<- conditionMessage(e)
        NULL
      }
    )
    if (is.null(roc)) {
      return(data.frame(
        year = y, AUC = NA_real_, lower95 = NA_real_, upper95 = NA_real_,
        eligible = FALSE, auc_estimable = FALSE, ci_estimable = FALSE,
        reason = "auc_estimation_failed", failure_class = "timeROC",
        error_message = roc_error, events_before = events_before, at_risk = at_risk,
        min_events_required = min_events, min_at_risk_required = min_at_risk,
        weighting = "marginal"
      ))
    }

    # timeROC prepends t=0 when only one target time is supplied. Select the
    # requested time explicitly rather than taking AUC[1], which is NA at t=0.
    auc_idx <- which.min(abs(roc$times - tt))
    auc_value <- unname(roc$AUC[auc_idx])
    ci_error <- ""
    ci <- tryCatch(
      confint(roc, level = 0.95)$CI_AUC,
      error = function(e) {
        ci_error <<- conditionMessage(e)
        NULL
      }
    )
    ci_lower <- if (is.null(ci)) NA_real_ else unname(ci[1, "2.5%"]) / 100
    ci_upper <- if (is.null(ci)) NA_real_ else unname(ci[1, "97.5%"]) / 100
    auc_ok <- is.finite(auc_value)
    ci_ok <- is.finite(ci_lower) && is.finite(ci_upper)
    data.frame(
      year = y, AUC = auc_value, lower95 = ci_lower, upper95 = ci_upper,
      eligible = auc_ok && ci_ok, auc_estimable = auc_ok, ci_estimable = ci_ok,
      reason = if (!auc_ok) "auc_estimation_failed" else if (!ci_ok) "ci_estimation_failed" else "ok",
      failure_class = if (!auc_ok) "timeROC" else if (!ci_ok) "confidence_interval" else "",
      error_message = if (!auc_ok) roc_error else if (!ci_ok) ci_error else "",
      events_before = events_before, at_risk = at_risk,
      min_events_required = min_events, min_at_risk_required = min_at_risk,
      weighting = "marginal"
    )
  })
  do.call(rbind, rows)
}

p_format <- function(x) ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1),
      axis.title = element_text(face = "plain"),
      legend.title = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey40"),
      strip.text = element_text(face = "bold")
    )
}

save_pub_plot <- function(plot, stem, width, height, dpi = 400) {
  ggsave(
    paste0(stem, ".pdf"), plot, width = width, height = height,
    units = "in", device = cairo_pdf
  )
  # ggplot2 may select ragg on this Windows host; ragg 1.5 can finish with a
  # warning but leave no file. Use the Cairo PNG device explicitly.
  png(
    paste0(stem, ".png"), width = width, height = height, units = "in",
    res = dpi, type = "cairo-png"
  )
  print(plot)
  dev.off()
  invisible(stem)
}

# -------------------------------------------------------------------------
# Strict nested-CV survival-signature functions.
#
# These functions deliberately accept training data only. In particular,
# fit_signature_pipeline() has no argument through which an outer validation
# partition can enter feature filtering, screening, scaling, tuning, or
# refitting.
# -------------------------------------------------------------------------

validate_pipeline_params <- function(p) {
  required <- c(
    "pipeline_version", "outer_folds", "outer_repeats", "inner_folds", "base_seed",
    "sd_cutoff", "univ_p_cutoff", "min_candidates",
    "candidate_fallback_n", "coef_cutoff", "min_selected",
    "selected_fallback_n", "max_selected", "lambda_ratio_grid",
    "lambda_rule", "ties"
  )
  missing <- setdiff(required, names(p))
  if (length(missing)) {
    stop("Pipeline parameters missing: ", paste(missing, collapse = ", "))
  }
  if (!is.character(p$pipeline_version) ||
      length(p$pipeline_version) != 1L || is.na(p$pipeline_version) ||
      !nzchar(trimws(p$pipeline_version))) {
    stop("pipeline_version must be one non-empty character value")
  }
  integer_fields <- c(
    "outer_folds", "outer_repeats", "inner_folds", "base_seed",
    "min_candidates", "candidate_fallback_n", "min_selected",
    "selected_fallback_n", "max_selected"
  )
  if (any(vapply(p[integer_fields], function(x) {
    length(x) != 1L || !is.finite(x) || x < 1 || x != as.integer(x)
  }, logical(1)))) {
    stop("Count and seed parameters must be positive finite integers")
  }
  if (p$outer_folds < 2L || p$inner_folds < 2L) {
    stop("outer_folds and inner_folds must both be at least 2")
  }
  if (!is.finite(p$sd_cutoff) || p$sd_cutoff < 0 ||
      !is.finite(p$univ_p_cutoff) || p$univ_p_cutoff < 0 ||
      p$univ_p_cutoff > 1 ||
      !is.finite(p$coef_cutoff) || p$coef_cutoff < 0) {
    stop("Filtering and coefficient cutoffs are outside their valid ranges")
  }
  if (p$candidate_fallback_n < p$min_candidates ||
      p$selected_fallback_n < p$min_selected ||
      p$max_selected < p$min_selected) {
    stop("Fallback and maximum feature counts are internally inconsistent")
  }
  if (!is.numeric(p$lambda_ratio_grid) || !length(p$lambda_ratio_grid) ||
      any(!is.finite(p$lambda_ratio_grid)) ||
      any(p$lambda_ratio_grid <= 0 | p$lambda_ratio_grid > 1)) {
    stop("lambda_ratio_grid must contain finite values in (0, 1]")
  }
  if (anyDuplicated(p$lambda_ratio_grid)) {
    stop("lambda_ratio_grid must not contain duplicates")
  }
  if (!identical(p$lambda_rule, "one_se_larger_penalty")) {
    stop("Unsupported lambda_rule: ", p$lambda_rule)
  }
  if (!p$ties %in% c("efron", "breslow")) {
    stop("Unsupported Cox ties method: ", p$ties)
  }
  invisible(TRUE)
}

validate_survival_training_data <- function(expr, time, status, context = "training") {
  expr <- as.matrix(expr)
  storage.mode(expr) <- "double"
  time <- as.numeric(time)
  status <- as.numeric(status)
  if (ncol(expr) != length(time) || length(time) != length(status)) {
    stop(context, ": expression columns, time, and status have different lengths")
  }
  if (nrow(expr) < 2L || ncol(expr) < 2L) {
    stop(context, ": at least two genes and two samples are required")
  }
  if (is.null(rownames(expr)) || anyNA(rownames(expr)) ||
      any(!nzchar(rownames(expr))) || anyDuplicated(rownames(expr))) {
    stop(context, ": gene row names must be nonempty and unique")
  }
  if (any(!is.finite(expr))) stop(context, ": non-finite expression value")
  if (any(!is.finite(time)) || any(time <= 0)) {
    stop(context, ": survival times must be finite and positive")
  }
  if (any(!status %in% 0:1) || length(unique(status)) != 2L) {
    stop(context, ": status must contain both 0 and 1 only")
  }
  invisible(list(expr = expr, time = time, status = status))
}

make_stratified_folds <- function(status, k, seed) {
  status <- as.numeric(status)
  k <- as.integer(k)
  if (length(status) < k || k < 2L || any(!status %in% 0:1) ||
      length(unique(status)) != 2L) {
    stop("Cannot create stratified folds from the supplied status vector")
  }
  counts <- table(factor(status, levels = 0:1))
  if (any(counts < k)) {
    stop("Each outcome class must contain at least k observations")
  }
  set.seed(as.integer(seed))
  out <- integer(length(status))
  for (s in 0:1) {
    idx <- sample(which(status == s), replace = FALSE)
    labels <- rep(seq_len(k), length.out = length(idx))
    out[idx] <- sample(labels, replace = FALSE)
  }
  if (any(tabulate(out, nbins = k) == 0L)) stop("An empty fold was generated")
  for (fold_id in seq_len(k)) {
    if (!setequal(unique(status[out == fold_id]), 0:1)) {
      stop("A stratified validation fold lacks events or censoring")
    }
  }
  out
}

univ_cox_table <- function(expr, time, status, ties = "efron") {
  surv_y <- Surv(time, status)
  n <- length(time)
  offsets <- rep(0, n)
  weights <- rep(1, n)
  row_ids <- as.character(seq_len(n))
  rows <- lapply(seq_len(nrow(expr)), function(i) {
    x <- matrix(as.numeric(expr[i, ]), ncol = 1L)
    fit <- tryCatch(
      suppressWarnings(coxph.fit(
        x = x,
        y = surv_y,
        strata = NULL,
        offset = offsets,
        init = NULL,
        control = coxph.control(),
        weights = weights,
        method = ties,
        rownames = row_ids,
        resid = FALSE
      )),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(c(HR = NA_real_, p = NA_real_))
    }
    beta <- unname(fit$coefficients[1L])
    variance <- unname(fit$var[1L, 1L])
    if (!is.finite(beta) || !is.finite(variance) || variance <= 0) {
      return(c(HR = NA_real_, p = NA_real_))
    }
    z <- beta / sqrt(variance)
    c(
      HR = suppressWarnings(exp(beta)),
      p = 2 * pnorm(abs(z), lower.tail = FALSE)
    )
  })
  tab <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  tab$gene <- rownames(expr)
  tab$fdr <- p.adjust(tab$p, method = "BH")
  tab <- tab[, c("gene", "HR", "p", "fdr"), drop = FALSE]
  tab[order(tab$p, tab$gene, na.last = TRUE), , drop = FALSE]
}

prepare_training <- function(expr, time, status, params) {
  validate_pipeline_params(params)
  checked <- validate_survival_training_data(expr, time, status)
  expr <- checked$expr
  time <- checked$time
  status <- checked$status

  sdev <- apply(expr, 1L, sd)
  keep <- is.finite(sdev) & sdev > params$sd_cutoff
  filtered <- expr[keep, , drop = FALSE]
  if (nrow(filtered) < params$min_candidates) {
    stop("Too few genes after training-partition variance filtering")
  }

  univ <- univ_cox_table(filtered, time, status, ties = params$ties)
  finite_ranked <- univ$gene[is.finite(univ$p)]
  candidates <- univ$gene[is.finite(univ$p) & univ$p < params$univ_p_cutoff]
  candidate_fallback <- length(candidates) < params$min_candidates
  if (candidate_fallback) {
    candidates <- head(finite_ranked, params$candidate_fallback_n)
  }
  if (length(candidates) < 2L) {
    stop("Too few finite univariable candidates for Cox LASSO")
  }

  candidate_expr <- filtered[candidates, , drop = FALSE]
  center <- rowMeans(candidate_expr)
  scale <- apply(candidate_expr, 1L, sd)
  if (any(!is.finite(center)) || any(!is.finite(scale) | scale <= 0)) {
    stop("Invalid training-partition standardization parameters")
  }
  z <- sweep(sweep(candidate_expr, 1L, center, "-"), 1L, scale, "/")
  if (any(!is.finite(z))) stop("Non-finite standardized training expression")

  list(
    z = z,
    center = center,
    scale = scale,
    univ_table = univ,
    univ_rank = finite_ranked,
    candidates = candidates,
    variance_genes = rownames(filtered),
    candidate_fallback = candidate_fallback
  )
}

select_refit_genes <- function(beta, candidate_rank, params) {
  beta <- beta[is.finite(beta)]
  if (is.null(names(beta)) || any(!nzchar(names(beta)))) {
    stop("LASSO coefficients must have gene names")
  }
  nonzero <- names(beta)[abs(beta) > 0]
  ordered_nonzero <- nonzero[order(-abs(beta[nonzero]), nonzero)]
  chosen <- ordered_nonzero[abs(beta[ordered_nonzero]) >= params$coef_cutoff]
  used_fallback <- length(chosen) < params$min_selected
  if (used_fallback) {
    ranked_candidates <- candidate_rank[candidate_rank %in% names(beta)]
    pool <- unique(c(ordered_nonzero, ranked_candidates))
    chosen <- head(pool, min(params$selected_fallback_n, length(pool)))
  }
  chosen <- head(chosen, params$max_selected)
  if (!length(chosen)) stop("Feature-selection fallback produced no genes")
  attr(chosen, "used_fallback") <- used_fallback
  chosen
}

fit_lasso_path <- function(prep, time, status, params) {
  suppressWarnings(glmnet(
    x = t(prep$z),
    y = Surv(time, status),
    family = "cox",
    alpha = 1,
    standardize = FALSE,
    nlambda = max(100L, length(params$lambda_ratio_grid)),
    lambda.min.ratio = min(params$lambda_ratio_grid)
  ))
}

restore_cox_gene_names <- function(fit, selected) {
  selected <- as.character(selected)
  fit_coef <- stats::coef(fit)
  fit_terms <- stats::terms(fit)
  term_labels <- attr(fit_terms, "term.labels")
  formula_genes <- all.vars(stats::delete.response(fit_terms))

  if (length(fit_coef) != length(selected) ||
      !identical(formula_genes, selected) ||
      !identical(names(fit_coef), term_labels)) {
    stop("Cox refit coefficient terms differ from selected genes")
  }

  # Formula processing encloses non-syntactic names such as KRTAP2-3 in
  # backticks.  The fit is unchanged; restore the original SYMBOL metadata
  # only after the one-to-one term/order checks above have passed.
  names(fit$coefficients) <- selected
  if (!is.null(fit$means)) {
    if (length(fit$means) != length(selected)) {
      stop("Cox refit means differ from selected genes")
    }
    names(fit$means) <- selected
  }
  if (!is.null(fit$x)) {
    if (ncol(fit$x) != length(selected)) {
      stop("Cox refit design matrix differs from selected genes")
    }
    colnames(fit$x) <- selected
  }
  if (!is.null(fit$assign)) {
    if (length(fit$assign) != length(selected)) {
      stop("Cox refit term assignments differ from selected genes")
    }
    names(fit$assign) <- selected
  }
  fit
}

fit_at_ratio <- function(prep, time, status, ratio, params, path = NULL) {
  if (is.null(path)) path <- fit_lasso_path(prep, time, status, params)
  lambda_max <- max(path$lambda)
  lambda <- lambda_max * ratio
  beta_matrix <- as.matrix(coef(path, s = lambda))
  beta <- as.numeric(beta_matrix[, 1L])
  beta_names <- rownames(beta_matrix)
  if (is.null(beta_names) && length(beta) == length(prep$candidates)) {
    beta_names <- prep$candidates
  }
  if (is.null(beta_names) || length(beta_names) != length(beta)) {
    stop("Unable to recover gene names from the Cox LASSO path")
  }
  names(beta) <- beta_names
  beta <- beta[prep$candidates]
  selected <- select_refit_genes(beta, prep$univ_rank, params)
  selected_fallback <- isTRUE(attr(selected, "used_fallback"))
  selected <- as.character(selected)

  model_dat <- data.frame(
    time = time,
    status = status,
    t(prep$z[selected, , drop = FALSE]),
    check.names = FALSE
  )
  refit <- suppressWarnings(coxph(
    Surv(time, status) ~ .,
    data = model_dat,
    ties = params$ties,
    x = TRUE,
    y = TRUE,
    singular.ok = FALSE
  ))
  refit <- restore_cox_gene_names(refit, selected)
  refit_coef <- coef(refit)
  if (!length(refit_coef) || any(!is.finite(refit_coef))) {
    stop("Non-finite coefficient in the unpenalized Cox refit")
  }
  if (!identical(names(refit_coef), selected)) {
    stop("Cox refit coefficient names differ from selected genes")
  }
  lp_train <- as.numeric(t(prep$z[selected, , drop = FALSE]) %*% refit_coef)
  lp_center <- mean(lp_train)
  lp_scale <- sd(lp_train)
  if (!is.finite(lp_center) || !is.finite(lp_scale) || lp_scale <= 0) {
    stop("Invalid training linear-predictor scale")
  }
  fingerprint <- digest(list(
    genes = selected,
    coefficients = unname(refit_coef),
    center = unname(prep$center[selected]),
    scale = unname(prep$scale[selected]),
    lp_center = lp_center,
    lp_scale = lp_scale,
    lambda_ratio = ratio,
    lambda = lambda
  ), algo = "sha256", serialize = TRUE)

  list(
    fit = refit,
    coefs = refit_coef,
    genes = selected,
    selected = selected,
    center = prep$center[selected],
    scale = prep$scale[selected],
    train_lp = lp_train,
    train_lp_mean = lp_center,
    train_lp_sd = lp_scale,
    lambda_ratio = ratio,
    lambda = lambda,
    lambda_max = lambda_max,
    lasso_coefficients = beta,
    lasso_nonzero = beta[abs(beta) > 0],
    selected_fallback = selected_fallback,
    model_fingerprint = fingerprint
  )
}

predict_signature_pipeline <- function(model, expr_new) {
  expr_new <- as.matrix(expr_new)
  storage.mode(expr_new) <- "double"
  missing_genes <- setdiff(model$genes, rownames(expr_new))
  if (length(missing_genes)) {
    stop("Prediction data missing genes: ", paste(missing_genes, collapse = ", "))
  }
  x <- expr_new[model$genes, , drop = FALSE]
  if (any(!is.finite(x))) stop("Non-finite expression value in prediction data")
  z <- sweep(sweep(x, 1L, model$center[model$genes], "-"),
             1L, model$scale[model$genes], "/")
  lp <- as.numeric(t(z) %*% model$coefs[model$genes])
  lp_train_sd <- (lp - model$train_lp_mean) / model$train_lp_sd
  if (any(!is.finite(lp_train_sd))) stop("Non-finite standardized prediction")
  list(lp = lp, lp_train_sd = lp_train_sd)
}

fit_signature_pipeline <- function(expr, time, status, seed, params) {
  validate_pipeline_params(params)
  checked <- validate_survival_training_data(expr, time, status)
  expr <- checked$expr
  time <- checked$time
  status <- checked$status
  if (any(table(factor(status, levels = 0:1)) < params$inner_folds)) {
    stop("Too few events or censored observations for the requested inner folds")
  }

  inner_foldid <- make_stratified_folds(status, params$inner_folds, seed)
  ratios <- sort(unique(params$lambda_ratio_grid), decreasing = TRUE)
  perf_rows <- vector("list", params$inner_folds * length(ratios))
  counter <- 0L
  for (inner_fold in seq_len(params$inner_folds)) {
    inner_train <- which(inner_foldid != inner_fold)
    inner_valid <- which(inner_foldid == inner_fold)
    prep <- tryCatch(
      prepare_training(
        expr[, inner_train, drop = FALSE],
        time[inner_train],
        status[inner_train],
        params
      ),
      error = function(e) {
        stop(
          "Inner fold ", inner_fold, " preprocessing failed: ",
          conditionMessage(e), call. = FALSE
        )
      }
    )
    path <- tryCatch(
      fit_lasso_path(prep, time[inner_train], status[inner_train], params),
      error = function(e) {
        stop(
          "Inner fold ", inner_fold, " Cox LASSO path failed: ",
          conditionMessage(e), call. = FALSE
        )
      }
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
      pred <- predict_signature_pipeline(
        candidate_model, expr[, inner_valid, drop = FALSE]
      )$lp_train_sd
      cidx <- unname(concordance(
        Surv(time[inner_valid], status[inner_valid]) ~ pred,
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
        model_fingerprint = candidate_model$model_fingerprint,
        stringsAsFactors = FALSE
      )
    }
  }
  inner_perf <- rbindlist(perf_rows, fill = TRUE)
  tuning <- inner_perf[, .(
    mean_c_index = mean(c_index),
    sd_c_index = sd(c_index),
    se_c_index = sd(c_index) / sqrt(.N),
    min_c_index = min(c_index),
    max_c_index = max(c_index),
    inner_folds = .N
  ), by = lambda_ratio]
  setorder(tuning, -lambda_ratio)
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

  full_prep <- prepare_training(expr, time, status, params)
  full_path <- fit_lasso_path(full_prep, time, status, params)
  final <- fit_at_ratio(
    full_prep, time, status, chosen_ratio, params, path = full_path
  )
  final$lasso_path <- full_path
  final$candidate_genes <- full_prep$candidates
  final$variance_genes <- full_prep$variance_genes
  final$univ_df <- full_prep$univ_table
  final$tuning_curve <- as.data.frame(tuning)
  final$inner_fold_performance <- as.data.frame(inner_perf)
  final$inner_fold_assignments <- data.frame(
    sample = if (is.null(colnames(expr))) as.character(seq_len(ncol(expr))) else colnames(expr),
    inner_fold = inner_foldid,
    stringsAsFactors = FALSE
  )
  final$params <- params
  final$seed <- as.integer(seed)
  final$audit <- data.frame(
    seed = as.integer(seed),
    input_genes = nrow(expr),
    input_samples = ncol(expr),
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
    stringsAsFactors = FALSE
  )
  final
}

model_set_jaccard <- function(a, b) {
  union_n <- length(union(a, b))
  if (!union_n) return(NA_real_)
  length(intersect(a, b)) / union_n
}

run_nested_cv <- function(expr, time, status, sample_ids, params, analysis_key) {
  validate_pipeline_params(params)
  checked <- validate_survival_training_data(expr, time, status)
  expr <- checked$expr
  time <- checked$time
  status <- checked$status
  sample_ids <- as.character(sample_ids)
  if (length(sample_ids) != ncol(expr) || anyNA(sample_ids) ||
      any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
    stop("sample_ids must be nonempty, unique, and match expression columns")
  }
  if (!is.character(analysis_key) || length(analysis_key) != 1L || !nzchar(analysis_key)) {
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
  fingerprint_rows <- list()
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
          !setequal(c(train_idx, test_idx), seq_len(ncol(expr)))) {
        stop("Outer training/test partitions are not a disjoint exhaustive split")
      }
      if (!setequal(unique(status[test_idx]), 0:1) ||
          !setequal(unique(status[train_idx]), 0:1)) {
        stop("Outer fold ", outer_fold, " lacks events or censoring")
      }
      fold_seed <- params$base_seed + repeat_id * 1000L + outer_fold
      model <- tryCatch(
        fit_signature_pipeline(
          expr[, train_idx, drop = FALSE], time[train_idx], status[train_idx],
          seed = fold_seed, params = params
        ),
        error = function(e) {
          stop(
            "Nested CV failed at repeat ", repeat_id, ", outer fold ", outer_fold,
            ": ", conditionMessage(e), call. = FALSE
          )
        }
      )
      pred <- predict_signature_pipeline(model, expr[, test_idx, drop = FALSE])
      cidx <- unname(concordance(
        Surv(time[test_idx], status[test_idx]) ~ pred$lp_train_sd,
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
          model$coefs[model$genes] > 0,
          "positive",
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
      fingerprint_rows[[counter]] <- data.frame(
        repeat_id = repeat_id,
        outer_fold = outer_fold,
        model_fingerprint = model$model_fingerprint,
        selected_lambda_ratio = model$lambda_ratio,
        analysis_key = analysis_key,
        stringsAsFactors = FALSE
      )
    }
  }

  fold_performance <- rbindlist(fold_rows, fill = TRUE)
  oof <- rbindlist(oof_rows, fill = TRUE)
  selection <- rbindlist(selection_rows, fill = TRUE)
  outer_assignments <- rbindlist(assignment_rows, fill = TRUE)
  inner_assignments <- rbindlist(inner_assignment_rows, fill = TRUE)
  tuning_performance <- rbindlist(tuning_rows, fill = TRUE)
  model_fingerprints <- rbindlist(fingerprint_rows, fill = TRUE)

  repeat_performance <- oof[, .(
    n = .N,
    events = sum(status),
    c_index = unname(concordance(
      Surv(time, status) ~ lp_oof, reverse = TRUE
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
    modal_sign = if (sum(coefficient > 0) >= sum(coefficient < 0)) "positive" else "negative",
    sign_consistency = max(sum(coefficient > 0), sum(coefficient < 0)) / .N,
    coefficient_mean = mean(coefficient),
    coefficient_sd = if (.N > 1L) sd(coefficient) else NA_real_,
    coefficient_median = median(coefficient),
    coefficient_min = min(coefficient),
    coefficient_max = max(coefficient)
  ), by = gene]
  setorder(selection_frequency, -frequency, -sign_consistency, gene)
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
  jaccard <- rbindlist(jaccard_rows, fill = TRUE)
  summary <- data.frame(
    pipeline_version = if (!is.null(params$pipeline_version)) {
      as.character(params$pipeline_version)
    } else {
      NA_character_
    },
    analysis_key = analysis_key,
    repeats = params$outer_repeats,
    outer_folds = params$outer_folds,
    inner_folds = params$inner_folds,
    completed_folds = nrow(fold_performance),
    median_repeat_oof_c_index = median(repeat_performance$c_index),
    repeat_oof_c_index_q1 = unname(quantile(repeat_performance$c_index, 0.25)),
    repeat_oof_c_index_q3 = unname(quantile(repeat_performance$c_index, 0.75)),
    min_repeat_oof_c_index = min(repeat_performance$c_index),
    max_repeat_oof_c_index = max(repeat_performance$c_index),
    median_within_repeat_jaccard = median(jaccard$jaccard),
    jaccard_q1 = unname(quantile(jaccard$jaccard, 0.25)),
    jaccard_q3 = unname(quantile(jaccard$jaccard, 0.75)),
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
    outer_fold_assignments = as.data.frame(outer_assignments),
    inner_fold_assignments = as.data.frame(inner_assignments),
    oof = as.data.frame(oof),
    fold_performance = as.data.frame(fold_performance),
    repeat_performance = as.data.frame(repeat_performance),
    selection = as.data.frame(selection),
    selection_frequency = as.data.frame(selection_frequency),
    jaccard = as.data.frame(jaccard),
    tuning_performance = as.data.frame(tuning_performance),
    model_fingerprints = as.data.frame(model_fingerprints),
    summary = summary,
    failure_log = failure_log
  )
  validate_nested_result(result, expr, time, status, sample_ids, params, analysis_key)
  result
}

summarize_repeated_oof <- function(oof, expected_repeats) {
  if (!is.numeric(expected_repeats) || length(expected_repeats) != 1L ||
      !is.finite(expected_repeats) || expected_repeats < 1L ||
      expected_repeats != as.integer(expected_repeats)) {
    stop("expected_repeats must be a positive integer")
  }
  required <- c("sample", "repeat_id", "time", "status", "lp_oof")
  missing <- setdiff(required, names(oof))
  if (length(missing)) {
    stop("OOF table missing columns: ", paste(missing, collapse = ", "))
  }
  dt <- as.data.table(oof)
  if (!nrow(dt) || anyNA(dt[, ..required]) ||
      any(!is.finite(dt$time)) || any(dt$time <= 0) ||
      any(!dt$status %in% 0:1) || any(!is.finite(dt$lp_oof))) {
    stop("OOF table contains invalid identifiers, outcomes, or predictions")
  }
  integrity <- dt[, .(
    n_rows = .N,
    n_unique_time = uniqueN(time),
    n_unique_status = uniqueN(status),
    n_unique_repeats = uniqueN(repeat_id)
  ), by = sample]
  if (any(integrity$n_rows != as.integer(expected_repeats)) ||
      any(integrity$n_unique_time != 1L) ||
      any(integrity$n_unique_status != 1L) ||
      any(integrity$n_unique_repeats != as.integer(expected_repeats))) {
    stop("Repeated OOF predictions are inconsistent within patient")
  }
  ensemble <- dt[, .(
    time = first(time),
    status = first(status),
    lp_oof_mean = mean(lp_oof),
    oof_repeats = .N
  ), by = sample]
  if (anyDuplicated(ensemble$sample) || nrow(ensemble) != nrow(integrity) ||
      any(!is.finite(ensemble$lp_oof_mean))) {
    stop("Repeated OOF predictions could not be reduced to one row per patient")
  }
  as.data.frame(ensemble)
}

validate_nested_result <- function(result, expr, time, status, sample_ids, params,
                                   analysis_key) {
  required <- c(
    "pipeline_version", "analysis_key", "params", "outer_fold_assignments", "inner_fold_assignments",
    "oof", "fold_performance", "repeat_performance", "selection",
    "selection_frequency", "jaccard", "tuning_performance",
    "model_fingerprints", "summary", "failure_log"
  )
  missing <- setdiff(required, names(result))
  if (length(missing)) stop("Nested-CV result missing fields: ", paste(missing, collapse = ", "))
  if (!identical(result$pipeline_version, params$pipeline_version)) {
    stop("Cached pipeline_version does not match current parameters")
  }
  if (!identical(result$analysis_key, analysis_key)) stop("Cached analysis key mismatch")
  if (!identical(result$params, params)) stop("Cached parameter object mismatch")
  expected_oof <- length(sample_ids) * params$outer_repeats
  expected_models <- params$outer_repeats * params$outer_folds
  if (nrow(result$oof) != expected_oof) stop("Incorrect number of OOF rows")
  key <- paste(result$oof$repeat_id, result$oof$sample, sep = "\r")
  if (anyDuplicated(key)) stop("A patient has multiple OOF predictions in one repeat")
  counts <- table(result$oof$repeat_id)
  if (length(counts) != params$outer_repeats || any(counts != length(sample_ids))) {
    stop("At least one repeat lacks complete patient-level OOF coverage")
  }
  if (any(!is.finite(result$oof$lp_oof))) stop("OOF predictions contain non-finite values")
  if (any(result$oof$analysis_key != analysis_key)) {
    stop("OOF rows contain an incorrect analysis key")
  }
  if (nrow(result$fold_performance) != expected_models) {
    stop("Incorrect number of completed outer folds")
  }
  if (nrow(result$repeat_performance) != params$outer_repeats ||
      anyDuplicated(result$repeat_performance$repeat_id) ||
      !setequal(result$repeat_performance$repeat_id, seq_len(params$outer_repeats)) ||
      any(!is.finite(result$repeat_performance$c_index))) {
    stop("Repeat-level performance table is incomplete or invalid")
  }
  if (nrow(result$model_fingerprints) != expected_models ||
      any(!nzchar(result$model_fingerprints$model_fingerprint))) {
    stop("Missing outer-fold model fingerprints")
  }
  expected_jaccard_rows <- params$outer_repeats * choose(params$outer_folds, 2L)
  required_jaccard <- c("repeat_id", "fold_a", "fold_b", "jaccard")
  if (!all(required_jaccard %in% names(result$jaccard)) ||
      nrow(result$jaccard) != expected_jaccard_rows ||
      anyDuplicated(result$jaccard[c("repeat_id", "fold_a", "fold_b")]) ||
      any(!is.finite(result$jaccard$jaccard)) ||
      any(result$jaccard$jaccard < 0 | result$jaccard$jaccard > 1)) {
    stop("Jaccard stability table is incomplete or invalid")
  }
  keyed_tables <- c(
    "outer_fold_assignments", "inner_fold_assignments", "oof",
    "fold_performance", "repeat_performance", "selection",
    "selection_frequency", "jaccard", "tuning_performance",
    "model_fingerprints"
  )
  for (table_name in keyed_tables) {
    tab <- result[[table_name]]
    if (!"analysis_key" %in% names(tab) || any(tab$analysis_key != analysis_key)) {
      stop("Incorrect analysis key in cached table: ", table_name)
    }
  }
  if (nrow(result$outer_fold_assignments) != expected_oof) {
    stop("Incorrect number of outer-fold assignment rows")
  }
  assignment_key <- paste(
    result$outer_fold_assignments$repeat_id,
    result$outer_fold_assignments$sample,
    sep = "\r"
  )
  if (anyDuplicated(assignment_key)) stop("Duplicate outer-fold assignment")
  for (repeat_id in seq_len(params$outer_repeats)) {
    assignment <- result$outer_fold_assignments[
      result$outer_fold_assignments$repeat_id == repeat_id, , drop = FALSE
    ]
    oof_repeat <- result$oof[result$oof$repeat_id == repeat_id, , drop = FALSE]
    if (!setequal(assignment$sample, sample_ids)) {
      stop("Outer-fold assignment omits or adds a source sample")
    }
    if (!setequal(oof_repeat$sample, sample_ids)) {
      stop("OOF predictions omit or add a source sample")
    }
    source_idx <- match(oof_repeat$sample, sample_ids)
    if (anyNA(source_idx)) stop("OOF sample identity is not present in source data")
    if (!isTRUE(all.equal(
      as.numeric(oof_repeat$time), as.numeric(time[source_idx]),
      tolerance = 0, check.attributes = FALSE
    ))) {
      stop("OOF survival times do not match source data")
    }
    if (!identical(as.numeric(oof_repeat$status), as.numeric(status[source_idx]))) {
      stop("OOF event indicators do not match source data")
    }
    assignment_idx <- match(oof_repeat$sample, assignment$sample)
    if (anyNA(assignment_idx) || any(
      oof_repeat$outer_fold != assignment$outer_fold[assignment_idx]
    )) {
      stop("OOF fold labels do not match outer-fold assignments")
    }
    repeat_row <- result$repeat_performance[
      result$repeat_performance$repeat_id == repeat_id, , drop = FALSE
    ]
    recomputed_c <- unname(concordance(
      Surv(oof_repeat$time, oof_repeat$status) ~ oof_repeat$lp_oof,
      reverse = TRUE
    )$concordance)
    if (nrow(repeat_row) != 1L || repeat_row$n != nrow(oof_repeat) ||
        repeat_row$events != sum(oof_repeat$status) ||
        !isTRUE(all.equal(
          as.numeric(repeat_row$c_index), recomputed_c,
          tolerance = 1e-12, check.attributes = FALSE
        ))) {
      stop("Repeat-level performance does not match patient-level OOF predictions")
    }
    for (outer_fold in seq_len(params$outer_folds)) {
      test_samples <- assignment$sample[assignment$outer_fold == outer_fold]
      train_samples <- setdiff(sample_ids, test_samples)
      idx <- match(test_samples, sample_ids)
      if (!setequal(unique(status[idx]), 0:1)) {
        stop("Cached outer test fold lacks events or censoring")
      }
      oof_fold <- oof_repeat[oof_repeat$outer_fold == outer_fold, , drop = FALSE]
      expected_events_test <- sum(status[idx])
      train_idx <- match(train_samples, sample_ids)
      expected_events_train <- sum(status[train_idx])
      if (!setequal(oof_fold$sample, test_samples) ||
          any(oof_fold$n_test != length(test_samples)) ||
          any(oof_fold$n_train != length(train_samples)) ||
          any(oof_fold$events_test != expected_events_test) ||
          any(oof_fold$events_train != expected_events_train) ||
          any(oof_fold$n_genes < 1L)) {
        stop("OOF fold audit counts do not match source partitions")
      }
      fold_row <- result$fold_performance[
        result$fold_performance$repeat_id == repeat_id &
          result$fold_performance$outer_fold == outer_fold,
        , drop = FALSE
      ]
      recomputed_fold_c <- unname(concordance(
        Surv(oof_fold$time, oof_fold$status) ~ oof_fold$lp_oof,
        reverse = TRUE
      )$concordance)
      if (nrow(fold_row) != 1L ||
          fold_row$n_test != length(test_samples) ||
          fold_row$n_train != length(train_samples) ||
          fold_row$events_test != expected_events_test ||
          fold_row$events_train != expected_events_train ||
          !isTRUE(all.equal(
            as.numeric(fold_row$c_index), recomputed_fold_c,
            tolerance = 1e-12, check.attributes = FALSE
          ))) {
        stop("Fold-level performance audit counts do not match source partitions")
      }
      inner <- result$inner_fold_assignments[
        result$inner_fold_assignments$repeat_id == repeat_id &
          result$inner_fold_assignments$outer_fold == outer_fold,
        , drop = FALSE
      ]
      if (!setequal(inner$sample, train_samples) || anyDuplicated(inner$sample) ||
          !setequal(unique(inner$inner_fold), seq_len(params$inner_folds))) {
        stop("Inner-fold assignments do not cover the outer training partition")
      }
      for (inner_fold in seq_len(params$inner_folds)) {
        inner_idx <- match(inner$sample[inner$inner_fold == inner_fold], sample_ids)
        if (!setequal(unique(status[inner_idx]), 0:1)) {
          stop("Cached inner validation fold lacks events or censoring")
        }
      }
    }
  }
  summary_fields <- c(
    "pipeline_version", "analysis_key", "repeats", "outer_folds", "inner_folds",
    "completed_folds", "median_repeat_oof_c_index",
    "repeat_oof_c_index_q1", "repeat_oof_c_index_q3",
    "min_repeat_oof_c_index", "max_repeat_oof_c_index",
    "median_within_repeat_jaccard", "jaccard_q1", "jaccard_q3"
  )
  if (!is.data.frame(result$summary) || nrow(result$summary) != 1L) {
    stop("Nested-CV summary must contain exactly one row")
  }
  missing_summary <- setdiff(summary_fields, names(result$summary))
  if (length(missing_summary)) {
    stop(
      "Nested-CV summary missing fields: ",
      paste(missing_summary, collapse = ", ")
    )
  }
  if (!identical(result$summary$pipeline_version, result$pipeline_version)) {
    stop("Nested-CV summary pipeline_version does not match result metadata")
  }
  repeat_c <- result$repeat_performance$c_index
  jaccard_values <- result$jaccard$jaccard
  expected_summary <- list(
    pipeline_version = if (!is.null(params$pipeline_version)) {
      as.character(params$pipeline_version)
    } else {
      NA_character_
    },
    analysis_key = analysis_key,
    repeats = params$outer_repeats,
    outer_folds = params$outer_folds,
    inner_folds = params$inner_folds,
    completed_folds = nrow(result$fold_performance),
    median_repeat_oof_c_index = median(repeat_c),
    repeat_oof_c_index_q1 = unname(quantile(repeat_c, 0.25)),
    repeat_oof_c_index_q3 = unname(quantile(repeat_c, 0.75)),
    min_repeat_oof_c_index = min(repeat_c),
    max_repeat_oof_c_index = max(repeat_c),
    median_within_repeat_jaccard = median(jaccard_values),
    jaccard_q1 = unname(quantile(jaccard_values, 0.25)),
    jaccard_q3 = unname(quantile(jaccard_values, 0.75))
  )
  for (field in summary_fields) {
    if (!identical(result$summary[[field]], expected_summary[[field]])) {
      stop("Nested-CV summary field mismatch: ", field)
    }
  }
  if (nrow(result$failure_log)) stop("A supposedly completed result contains failures")
  invisible(TRUE)
}

sha256_file <- function(path) {
  if (!file.exists(path)) stop("Hash input does not exist: ", path)
  digest(file = path, algo = "sha256")
}

make_analysis_key <- function(input_paths, code_paths, params,
                              return_manifest = FALSE) {
  validate_pipeline_params(params)
  input_paths <- normalizePath(input_paths, winslash = "/", mustWork = TRUE)
  code_paths <- normalizePath(code_paths, winslash = "/", mustWork = TRUE)
  input_names <- basename(input_paths)
  code_names <- basename(code_paths)
  if (anyDuplicated(input_names) || anyDuplicated(code_names)) {
    stop("Hash inputs must have unique basenames")
  }
  packages <- c("survival", "glmnet", "data.table", "digest")
  manifest <- list(
    inputs = setNames(vapply(input_paths, sha256_file, character(1)), input_names),
    code = setNames(vapply(code_paths, sha256_file, character(1)), code_names),
    params = params,
    r_version = R.version.string,
    packages = setNames(vapply(
      packages,
      function(x) as.character(utils::packageVersion(x)),
      character(1)
    ), packages)
  )
  key <- digest(manifest, algo = "sha256", serialize = TRUE)
  if (isTRUE(return_manifest)) list(key = key, manifest = manifest) else key
}

canonicalize_for_hash <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    nm <- names(x)
    if (!is.null(nm)) {
      if (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm)) {
        stop("Named hash lists must have unique, non-empty names")
      }
      x <- x[order(nm)]
    }
    return(lapply(x, canonicalize_for_hash))
  }
  x
}

normalize_identity_paths <- function(paths, label) {
  ids <- names(paths)
  paths <- as.character(paths)
  if (!length(paths)) stop(label, " must contain at least one file")
  if (is.null(ids)) {
    ids <- basename(paths)
  } else if (anyNA(ids) || any(!nzchar(ids))) {
    stop(label, " names must all be non-empty when supplied")
  }
  if (anyDuplicated(ids)) stop(label, " identifiers must be unique")
  ord <- order(ids)
  list(
    ids = ids[ord],
    paths = normalizePath(paths[ord], winslash = "/", mustWork = TRUE)
  )
}

make_run_key <- function(input_paths, code_paths, params,
                         return_manifest = FALSE) {
  inputs <- normalize_identity_paths(input_paths, "Run-key inputs")
  code <- normalize_identity_paths(code_paths, "Run-key code")
  packages <- c(
    "survival", "glmnet", "timeROC", "data.table", "digest", "ggplot2",
    "AnnotationDbi", "hgu133plus2.db", "org.Hs.eg.db", "metafor"
  )
  manifest <- list(
    schema_version = "crc9-run-identity-v1",
    inputs = setNames(
      vapply(inputs$paths, sha256_file, character(1)), inputs$ids
    ),
    code = setNames(
      vapply(code$paths, sha256_file, character(1)), code$ids
    ),
    params = canonicalize_for_hash(params),
    r_version = R.version.string,
    r_platform = R.version$platform,
    packages = setNames(vapply(
      packages,
      function(x) as.character(utils::packageVersion(x)),
      character(1)
    ), packages)
  )
  key <- digest(manifest, algo = "sha256", serialize = TRUE)
  if (isTRUE(return_manifest)) list(key = key, manifest = manifest) else key
}

write_dput_manifest <- function(object, path) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("Manifest path must be one non-empty value")
  }
  parent <- dirname(path)
  if (!dir.exists(parent)) stop("Manifest parent directory does not exist: ", parent)
  writeLines(capture.output(dput(object)), path, useBytes = TRUE)
  restored <- dget(path)
  if (!identical(restored, object)) {
    stop("Serialized parameter manifest does not round-trip exactly")
  }
  invisible(path)
}

snapshot_code_files <- function(code_paths, run_root,
                                relative_dir = "code_snapshot") {
  code <- normalize_identity_paths(code_paths, "Code snapshot inputs")
  run_root <- normalizePath(run_root, winslash = "/", mustWork = TRUE)
  relative_dir <- validate_relative_artifact_path(relative_dir)
  if (grepl("/", relative_dir, fixed = TRUE)) {
    stop("Code snapshot directory must be a single relative directory name")
  }
  source_basenames <- basename(code$paths)
  if (anyDuplicated(source_basenames)) {
    stop("Code snapshot source basenames must be unique")
  }
  snapshot_dir <- file.path(run_root, relative_dir)
  if (file.exists(snapshot_dir) || dir.exists(snapshot_dir)) {
    stop("Code snapshot directory already exists: ", snapshot_dir)
  }
  if (!dir.create(snapshot_dir, recursive = FALSE, showWarnings = FALSE)) {
    stop("Could not create code snapshot directory: ", snapshot_dir)
  }
  relative_paths <- file.path(relative_dir, source_basenames)
  destination_paths <- file.path(run_root, relative_paths)
  copied <- file.copy(
    code$paths, destination_paths,
    overwrite = FALSE, copy.mode = TRUE, copy.date = TRUE
  )
  if (length(copied) != length(code$paths) || any(!copied)) {
    stop("Failed to copy every code file into the run snapshot")
  }
  source_sha256 <- vapply(code$paths, sha256_file, character(1))
  snapshot_sha256 <- vapply(destination_paths, sha256_file, character(1))
  if (!identical(unname(source_sha256), unname(snapshot_sha256))) {
    stop("Code snapshot hashes do not match the hashed source files")
  }
  data.frame(
    name = code$ids,
    source_basename = source_basenames,
    relative_path = gsub("\\\\", "/", relative_paths),
    bytes = unname(file.info(destination_paths)$size),
    sha256 = unname(snapshot_sha256),
    stringsAsFactors = FALSE
  )
}

summarize_reml_meta_fit <- function(
    fit, inference, primary, cohort_inputs,
    analysis = "unadjusted external OS association", endpoint = "OS") {
  if (!inherits(fit, "rma")) stop("fit must be a metafor rma object")
  cohort_inputs <- sort(unique(as.character(cohort_inputs)))
  if (length(cohort_inputs) != fit$k || anyNA(cohort_inputs) ||
      any(!nzchar(cohort_inputs))) {
    stop("Meta-analysis cohort identifiers must be unique, complete, and match k")
  }
  pred <- predict(fit)
  inference_df <- if (!is.null(fit$dfs) && length(fit$dfs) == 1L &&
      is.finite(fit$dfs)) as.numeric(fit$dfs) else NA_real_
  data.frame(
    endpoint = endpoint,
    cohorts = fit$k,
    pooled_HR = exp(as.numeric(fit$b)),
    lower95 = exp(fit$ci.lb),
    upper95 = exp(fit$ci.ub),
    p = fit$pval,
    prediction_lower95 = exp(pred$pi.lb),
    prediction_upper95 = exp(pred$pi.ub),
    Q = fit$QE,
    Q_df = as.integer(fit$k - fit$p),
    Q_p = fit$QEp,
    I2 = fit$I2,
    tau2 = fit$tau2,
    test_statistic = as.numeric(fit$zval),
    inference_df = inference_df,
    inference = as.character(inference),
    primary = isTRUE(primary),
    method = "REML random-effects inverse-variance meta-analysis",
    cohort_inputs = paste(cohort_inputs, collapse = ";"),
    analysis = as.character(analysis),
    stringsAsFactors = FALSE
  )
}

fit_reml_meta_dual <- function(data,
                               analysis = "unadjusted external OS association",
                               endpoint = "OS") {
  required <- c("cohort", "yi", "sei")
  if (!is.data.frame(data) || !all(required %in% names(data))) {
    stop("Meta-analysis input must contain cohort, yi, and sei")
  }
  d <- as.data.frame(data[, required], stringsAsFactors = FALSE)
  d$cohort <- as.character(d$cohort)
  d$yi <- as.numeric(d$yi)
  d$sei <- as.numeric(d$sei)
  if (nrow(d) < 2L || anyNA(d$cohort) || any(!nzchar(d$cohort)) ||
      anyDuplicated(d$cohort) || any(!is.finite(d$yi)) ||
      any(!is.finite(d$sei)) || any(d$sei <= 0)) {
    stop("Meta-analysis inputs must be at least two unique cohorts with finite effects and positive SEs")
  }
  fit_kh <- metafor::rma.uni(
    yi = d$yi, sei = d$sei, method = "REML", test = "knha"
  )
  fit_wald <- metafor::rma.uni(
    yi = d$yi, sei = d$sei, method = "REML"
  )
  rbind(
    summarize_reml_meta_fit(
      fit_kh, "Knapp-Hartung/t", TRUE, d$cohort,
      analysis = analysis, endpoint = endpoint
    ),
    summarize_reml_meta_fit(
      fit_wald, "normal/Wald", FALSE, d$cohort,
      analysis = analysis, endpoint = endpoint
    )
  )
}

leave_one_out_reml_meta_dual <- function(
    data, analysis = "leave-one-out unadjusted external OS association",
    endpoint = "OS") {
  if (!is.data.frame(data) || !all(c("cohort", "yi", "sei") %in% names(data))) {
    stop("Leave-one-out input must contain cohort, yi, and sei")
  }
  cohorts <- sort(unique(as.character(data$cohort)))
  if (length(cohorts) < 3L || length(cohorts) != nrow(data)) {
    stop("Leave-one-out meta-analysis requires at least three unique cohort rows")
  }
  out <- lapply(cohorts, function(drop_cohort) {
    z <- data[as.character(data$cohort) != drop_cohort, , drop = FALSE]
    cbind(
      dropped_cohort = drop_cohort,
      fit_reml_meta_dual(z, analysis = analysis, endpoint = endpoint),
      stringsAsFactors = FALSE
    )
  })
  rownames(out) <- NULL
  data.table::rbindlist(out, fill = TRUE)
}

validate_relative_artifact_path <- function(path) {
  path <- gsub("\\\\", "/", as.character(path))
  if (length(path) != 1L || is.na(path) || !nzchar(path) ||
      grepl("^[A-Za-z]:/|^/", path) ||
      grepl("(^|/)\\.\\.(/|$)", path)) {
    stop("Artifact paths must be safe relative paths: ", path)
  }
  path
}

new_artifact_registry <- function(run_root, spec) {
  if (!dir.exists(run_root)) stop("Run root does not exist: ", run_root)
  required_columns <- c("artifact_id", "relative_path", "required")
  if (!is.data.frame(spec) || !all(required_columns %in% names(spec))) {
    stop("Artifact specification lacks required columns")
  }
  spec <- as.data.frame(spec[, required_columns], stringsAsFactors = FALSE)
  spec$artifact_id <- as.character(spec$artifact_id)
  spec$relative_path <- vapply(
    spec$relative_path, validate_relative_artifact_path, character(1)
  )
  spec$required <- as.logical(spec$required)
  if (!nrow(spec) || anyNA(spec$artifact_id) || any(!nzchar(spec$artifact_id)) ||
      anyDuplicated(spec$artifact_id) || anyDuplicated(spec$relative_path) ||
      anyNA(spec$required)) {
    stop("Artifact specification contains invalid or duplicated entries")
  }
  registry <- new.env(parent = emptyenv())
  registry$run_root <- normalizePath(run_root, winslash = "/", mustWork = TRUE)
  registry$records <- transform(
    spec,
    status = "pending",
    reason = "",
    bytes = NA_real_,
    sha256 = NA_character_
  )
  class(registry) <- "crc_artifact_registry"
  registry
}

artifact_path <- function(registry, artifact_id) {
  if (!inherits(registry, "crc_artifact_registry")) {
    stop("Invalid artifact registry")
  }
  idx <- match(artifact_id, registry$records$artifact_id)
  if (is.na(idx)) stop("Artifact is not in the allow-list: ", artifact_id)
  file.path(registry$run_root, registry$records$relative_path[idx])
}

register_artifact_file <- function(registry, artifact_id) {
  idx <- match(artifact_id, registry$records$artifact_id)
  if (is.na(idx)) stop("Artifact is not in the allow-list: ", artifact_id)
  if (!identical(registry$records$status[idx], "pending")) {
    stop("Artifact was already resolved: ", artifact_id)
  }
  path <- artifact_path(registry, artifact_id)
  if (!file.exists(path) || dir.exists(path)) {
    stop("Registered artifact file does not exist: ", path)
  }
  path_norm <- normalizePath(path, winslash = "/", mustWork = TRUE)
  root_prefix <- paste0(registry$run_root, "/")
  if (!startsWith(path_norm, root_prefix)) {
    stop("Artifact resolves outside the run root: ", artifact_id)
  }
  registry$records$status[idx] <- "generated"
  registry$records$bytes[idx] <- unname(file.info(path)$size)
  registry$records$sha256[idx] <- sha256_file(path)
  invisible(path)
}

register_not_generated <- function(registry, artifact_id, reason) {
  idx <- match(artifact_id, registry$records$artifact_id)
  if (is.na(idx)) stop("Artifact is not in the allow-list: ", artifact_id)
  if (isTRUE(registry$records$required[idx])) {
    stop("A required artifact cannot be marked not_generated: ", artifact_id)
  }
  if (!identical(registry$records$status[idx], "pending")) {
    stop("Artifact was already resolved: ", artifact_id)
  }
  reason <- as.character(reason)
  if (length(reason) != 1L || is.na(reason) || !nzchar(trimws(reason))) {
    stop("not_generated artifacts require a non-empty reason")
  }
  registry$records$status[idx] <- "not_generated"
  registry$records$reason[idx] <- reason
  invisible(TRUE)
}

finalize_artifact_registry <- function(registry, manifest_path,
                                       control_relative_paths = character()) {
  if (!inherits(registry, "crc_artifact_registry")) {
    stop("Invalid artifact registry")
  }
  for (artifact_id in registry$records$artifact_id[
    registry$records$status == "pending"
  ]) {
    path <- artifact_path(registry, artifact_id)
    if (file.exists(path) && !dir.exists(path)) {
      register_artifact_file(registry, artifact_id)
    }
  }
  unresolved_required <- registry$records$artifact_id[
    registry$records$required & registry$records$status != "generated"
  ]
  if (length(unresolved_required)) {
    stop("Required artifacts were not generated: ",
         paste(unresolved_required, collapse = ", "))
  }
  unresolved_optional <- registry$records$artifact_id[
    !registry$records$required & registry$records$status == "pending"
  ]
  if (length(unresolved_optional)) {
    stop("Optional artifacts lack generated/not_generated status: ",
         paste(unresolved_optional, collapse = ", "))
  }
  control_relative_paths <- vapply(
    control_relative_paths, validate_relative_artifact_path, character(1)
  )
  actual <- list.files(
    registry$run_root, recursive = TRUE, full.names = TRUE,
    all.files = TRUE, no.. = TRUE
  )
  actual <- actual[file.exists(actual) & !dir.exists(actual)]
  root_prefix <- paste0(registry$run_root, "/")
  actual_norm <- normalizePath(actual, winslash = "/", mustWork = TRUE)
  actual_relative <- substring(actual_norm, nchar(root_prefix) + 1L)
  allowed_relative <- c(
    registry$records$relative_path[registry$records$status == "generated"],
    control_relative_paths
  )
  unexpected <- setdiff(actual_relative, allowed_relative)
  if (length(unexpected)) {
    stop("Unregistered files exist in the run directory: ",
         paste(unexpected, collapse = ", "))
  }
  manifest_path <- gsub("\\\\", "/", manifest_path)
  manifest_parent <- normalizePath(
    dirname(manifest_path), winslash = "/", mustWork = TRUE
  )
  if (!startsWith(paste0(manifest_parent, "/"), root_prefix)) {
    stop("Output manifest must be written inside the run root")
  }
  manifest <- registry$records[, c(
    "artifact_id", "relative_path", "required", "status", "reason",
    "bytes", "sha256"
  )]
  data.table::fwrite(manifest, manifest_path)
  manifest
}

write_run_status_file <- function(
    path, status, run_key, nested_analysis_key = NA_character_,
    input_manifest_sha256 = NA_character_,
    parameter_manifest_sha256 = NA_character_,
    code_manifest_sha256 = NA_character_,
    output_manifest_sha256 = NA_character_,
    started_at_utc = NA_character_, detail = "") {
  status <- match.arg(status, c("running", "failed_or_interrupted", "complete"))
  sha_pattern <- "^[0-9a-f]{64}$"
  if (length(run_key) != 1L || is.na(run_key) || !grepl(sha_pattern, run_key)) {
    stop("run_key must be a SHA-256 hex digest")
  }
  if (identical(status, "complete")) {
    required_hashes <- c(
      input_manifest_sha256, parameter_manifest_sha256,
      code_manifest_sha256, output_manifest_sha256
    )
    if (length(nested_analysis_key) != 1L || is.na(nested_analysis_key) ||
        !grepl(sha_pattern, nested_analysis_key) ||
        anyNA(required_hashes) || any(!grepl(sha_pattern, required_hashes))) {
      stop("Complete status requires nested key and four manifest hashes")
    }
  }
  now_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  data.table::fwrite(data.frame(
    status = status,
    run_key = run_key,
    nested_analysis_key = nested_analysis_key,
    input_manifest_sha256 = input_manifest_sha256,
    parameter_manifest_sha256 = parameter_manifest_sha256,
    code_manifest_sha256 = code_manifest_sha256,
    output_manifest_sha256 = output_manifest_sha256,
    started_at_utc = started_at_utc,
    finished_at_utc = if (identical(status, "running")) NA_character_ else now_utc,
    detail = as.character(detail),
    stringsAsFactors = FALSE
  ), path)
  invisible(path)
}

commit_run_directory <- function(tmp_dir, final_dir, run_key) {
  if (!dir.exists(tmp_dir)) stop("Temporary run directory does not exist")
  if (file.exists(final_dir) || dir.exists(final_dir)) {
    stop("Final run directory already exists; refusing to overwrite: ", final_dir)
  }
  if (!identical(basename(tmp_dir), paste0(run_key, ".tmp")) ||
      !identical(basename(final_dir), run_key)) {
    stop("Temporary/final directory names do not match run_key")
  }
  tmp_parent <- normalizePath(dirname(tmp_dir), winslash = "/", mustWork = TRUE)
  final_parent <- normalizePath(dirname(final_dir), winslash = "/", mustWork = TRUE)
  if (!identical(tmp_parent, final_parent)) {
    stop("Temporary and final run directories must share one parent directory")
  }
  if (!file.rename(tmp_dir, final_dir)) {
    stop("Atomic run-directory rename failed")
  }
  if (dir.exists(tmp_dir) || !dir.exists(final_dir)) {
    stop("Run-directory rename did not produce the expected final state")
  }
  invisible(final_dir)
}
