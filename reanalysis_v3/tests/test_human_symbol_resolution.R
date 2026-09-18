#!/usr/bin/env Rscript

# Regression tests for outcome-blind resolution of current human SYMBOLs to
# legacy expression row labels. All decision tests are synthetic except the
# ATP23/XRCC6BP1 annotation identity check against frozen org.Hs.eg.db.

options(stringsAsFactors = FALSE, warn = 1)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript")
this_file <- gsub("\\\\", "/", sub("^--file=", "", script_arg[[1L]]))
test_dir <- dirname(this_file)
v3_root <- dirname(test_dir)
work_root <- dirname(v3_root)
utils_file <- file.path(v3_root, "scripts", "utils.R")
if (!file.exists(utils_file)) stop("Cannot locate utils.R relative to test file")
Sys.setenv(CRC_WORK_ROOT = work_root)
if (!nzchar(Sys.getenv("CRC_SOURCE_ROOT"))) {
  Sys.setenv(CRC_SOURCE_ROOT = work_root)
}
source(utils_file)

expect_error <- function(expr, pattern) {
  captured <- NULL
  tryCatch(force(expr), error = function(e) captured <<- conditionMessage(e))
  if (is.null(captured)) stop("Expected an error but expression returned normally")
  if (!grepl(pattern, captured, ignore.case = TRUE)) {
    stop("Error did not match /", pattern, "/: ", captured)
  }
  invisible(captured)
}

make_expr <- function(row_ids) {
  values <- seq_len(length(row_ids) * 3L)
  out <- matrix(values, nrow = length(row_ids), byrow = TRUE)
  rownames(out) <- row_ids
  colnames(out) <- c("S1", "S2", "S3")
  out
}

make_crosswalk <- function(
    gene, entrez, aliases, eligible,
    reverse_unique = eligible,
    symbol_consistent = eligible,
    entrez_consistent = eligible) {
  stopifnot(
    length(aliases) == length(eligible),
    length(reverse_unique) == length(aliases),
    length(symbol_consistent) == length(aliases),
    length(entrez_consistent) == length(aliases)
  )
  data.frame(
    required_gene = rep(gene, length(aliases)),
    required_entrez_id = rep(entrez, length(aliases)),
    alias = aliases,
    reverse_symbols = ifelse(symbol_consistent, gene, paste0(gene, ";OTHER")),
    reverse_entrez_ids = ifelse(entrez_consistent, entrez, paste0(entrez, ";999")),
    reverse_symbol_count = ifelse(reverse_unique, 1L, 2L),
    reverse_entrez_count = ifelse(reverse_unique, 1L, 2L),
    reverse_unique = reverse_unique,
    symbol_consistent = symbol_consistent,
    entrez_consistent = entrez_consistent,
    eligible_alias = eligible,
    annotation_package = "synthetic.annotation",
    annotation_version = "1.0",
    stringsAsFactors = FALSE
  )
}

atp_crosswalk <- build_human_symbol_crosswalk(c("ATP23", "TP53"))
atp_alias <- atp_crosswalk[
  atp_crosswalk$required_gene == "ATP23" &
    atp_crosswalk$alias == "XRCC6BP1",
  , drop = FALSE
]
stopifnot(
  nrow(atp_alias) == 1L,
  identical(atp_alias$required_entrez_id, "91419"),
  isTRUE(atp_alias$reverse_unique),
  isTRUE(atp_alias$symbol_consistent),
  isTRUE(atp_alias$entrez_consistent),
  isTRUE(atp_alias$eligible_alias)
)

# Exact current SYMBOL must win even when an eligible historical alias row is
# also present in the expression matrix.
exact_expr <- make_expr(c("XRCC6BP1", "ATP23"))
exact <- resolve_expression_gene_rows(
  exact_expr, "ATP23", atp_crosswalk, cohort = "synthetic exact"
)
stopifnot(
  identical(rownames(exact$expr), "ATP23"),
  identical(as.numeric(exact$expr["ATP23", ]), as.numeric(exact_expr["ATP23", ])),
  identical(exact$audit$source_row_id, "ATP23"),
  identical(exact$audit$resolution_method, "exact_SYMBOL")
)

# When the exact current SYMBOL is absent, ATP23 must resolve to the sole
# reverse-unique, Entrez-consistent XRCC6BP1 row.
alias_expr <- make_expr(c("XRCC6BP1", "UNRELATED"))
alias <- resolve_expression_gene_rows(
  alias_expr, "ATP23", atp_crosswalk, cohort = "synthetic alias"
)
stopifnot(
  identical(rownames(alias$expr), "ATP23"),
  identical(alias$audit$source_row_id, "XRCC6BP1"),
  identical(
    alias$audit$resolution_method,
    "unique_reverse_unique_ENTREZ_alias"
  ),
  identical(alias$audit$required_entrez_id, "91419"),
  isTRUE(alias$audit$reverse_unique),
  isTRUE(alias$audit$symbol_consistent),
  isTRUE(alias$audit$entrez_consistent),
  identical(alias$audit$reverse_symbols, "ATP23"),
  identical(alias$audit$reverse_entrez_ids, "91419"),
  identical(alias$audit$reverse_symbol_count, 1L),
  identical(alias$audit$reverse_entrez_count, 1L),
  identical(
    as.numeric(alias$expr["ATP23", ]),
    as.numeric(alias_expr["XRCC6BP1", ])
  )
)

expect_error(
  resolve_expression_gene_rows(
    make_expr("UNRELATED"), "ATP23", atp_crosswalk,
    cohort = "synthetic missing"
  ),
  "no exact SYMBOL or annotated historical alias"
)

multiple_crosswalk <- make_crosswalk(
  "GENE1", "101", c("OLD1", "OLD2"), c(TRUE, TRUE)
)
expect_error(
  resolve_expression_gene_rows(
    make_expr(c("OLD1", "OLD2")), "GENE1", multiple_crosswalk,
    cohort = "synthetic multiple"
  ),
  "exactly one eligible historical alias"
)

ambiguous_crosswalk <- make_crosswalk(
  "GENE2", "202", "AMBIG", FALSE,
  reverse_unique = FALSE,
  symbol_consistent = FALSE,
  entrez_consistent = FALSE
)
expect_error(
  resolve_expression_gene_rows(
    make_expr("AMBIG"), "GENE2", ambiguous_crosswalk,
    cohort = "synthetic ambiguous"
  ),
  "ambiguous or ENTREZ-inconsistent"
)

# Both the remapped matrix and audit table must preserve the caller-supplied
# required-gene order, not annotation order or expression-row order.
ordered <- resolve_expression_gene_rows(
  make_expr(c("XRCC6BP1", "TP53")),
  c("TP53", "ATP23"),
  atp_crosswalk,
  cohort = "synthetic order"
)
stopifnot(
  identical(rownames(ordered$expr), c("TP53", "ATP23")),
  identical(ordered$audit$required_gene, c("TP53", "ATP23")),
  identical(ordered$audit$source_row_id, c("TP53", "XRCC6BP1"))
)

cat("Human SYMBOL/alias resolution tests: PASS\n")
