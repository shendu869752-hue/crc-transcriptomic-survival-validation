#!/usr/bin/env Rscript

# Focused regression test for one-to-one GEO probe mapping. Uses real GPL570
# probe IDs but a synthetic expression matrix, so it is fast and outcome-free.

options(stringsAsFactors = FALSE)
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

probe_ids <- c(
  "220077_at",     # CCDC134, one-to-one
  "235434_at",     # CCDC134, one-to-one
  "218746_at",     # TAPBPL, one-to-one
  "218747_s_at",   # TAPBPL, one-to-one
  "1552283_s_at", # ZDHHC11/ZDHHC11B, one-to-many
  "AFFX-BioB-3_at" # no human SYMBOL
)
probe_mat <- rbind(
  c(2, 3, 4),
  c(4, 3, 2), # ties 220077_at on mean; lexicographic PROBEID must win
  c(1, 1, 1),
  c(4, 4, 4), # higher-mean TAPBPL probe
  c(9, 9, 9),
  c(8, 8, 8)
)
rownames(probe_mat) <- probe_ids
colnames(probe_mat) <- c("S1", "S2", "S3")

highest <- map_probes(
  probe_mat,
  method = "unique_highest_mean",
  required_genes = c("CCDC134", "TAPBPL")
)
averaged <- map_probes(
  probe_mat,
  method = "unique_mean",
  required_genes = c("CCDC134", "TAPBPL")
)

stopifnot(
  setequal(rownames(highest$expr), c("CCDC134", "TAPBPL")),
  setequal(rownames(averaged$expr), c("CCDC134", "TAPBPL")),
  highest$mapping$actual_probe[highest$mapping$SYMBOL == "CCDC134"] == "220077_at",
  highest$mapping$actual_probe[highest$mapping$SYMBOL == "TAPBPL"] == "218747_s_at",
  highest$mapping$n_probes_gene[highest$mapping$SYMBOL == "CCDC134"] == 2L,
  highest$mapping$n_probes_gene[highest$mapping$SYMBOL == "TAPBPL"] == 2L,
  identical(as.numeric(averaged$expr["CCDC134", ]), c(3, 3, 3)),
  identical(as.numeric(averaged$expr["TAPBPL", ]), c(2.5, 2.5, 2.5))
)

multi <- highest$probe_audit[highest$probe_audit$PROBEID == "1552283_s_at", ]
no_symbol <- highest$probe_audit[highest$probe_audit$PROBEID == "AFFX-BioB-3_at", ]
ccdc <- highest$probe_audit[
  !is.na(highest$probe_audit$SYMBOL) &
    highest$probe_audit$SYMBOL == "CCDC134",
  ,
  drop = FALSE
]
stopifnot(
  nrow(multi) == 1L,
  multi$n_symbols_per_probe == 2L,
  !multi$eligible,
  !multi$used_in_expression,
  multi$selection_reason == "ineligible_multiple_SYMBOLs",
  nrow(no_symbol) == 1L,
  no_symbol$n_symbols_per_probe == 0L,
  !no_symbol$eligible,
  no_symbol$selection_reason == "ineligible_no_SYMBOL",
  nrow(ccdc) == 2L,
  all(ccdc$n_probes_gene == 2L),
  sum(ccdc$used_in_expression) == 1L,
  ccdc$selection_reason[ccdc$PROBEID == "235434_at"] ==
    "eligible_not_highest_mean",
  all(averaged$probe_audit$used_in_expression[averaged$probe_audit$eligible])
)

missing_error <- tryCatch(
  {
    map_probes(
      probe_mat,
      method = "unique_highest_mean",
      required_genes = "ZDHHC11"
    )
    NULL
  },
  error = function(e) conditionMessage(e)
)
stopifnot(
  is.character(missing_error),
  grepl("Required signature genes missing", missing_error, fixed = TRUE),
  grepl("ZDHHC11", missing_error, fixed = TRUE)
)

cat("Unique probe-mapping tests: PASS\n")
