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
launcher_file <- file.path(work_root, "run_v3_local.ps1")

if (!all(file.exists(c(utils_file, core_file, launcher_file)))) {
  stop("Cannot locate run-isolation source files relative to test file")
}

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

make_spec <- function(required, optional = character()) {
  paths <- c(required, optional)
  data.frame(
    artifact_id = paths,
    relative_path = paths,
    required = c(rep(TRUE, length(required)), rep(FALSE, length(optional))),
    stringsAsFactors = FALSE
  )
}

run_test("run key changes for input, code, and critical parameters", {
  td1 <- tempfile("crc9-key-a-")
  td2 <- tempfile("crc9-key-b-")
  dir.create(td1)
  dir.create(td2)
  on.exit(unlink(c(td1, td2), recursive = TRUE, force = TRUE), add = TRUE)
  for (td in c(td1, td2)) {
    writeLines("training-input", file.path(td, "training.csv"), useBytes = TRUE)
    writeLines("external-input", file.path(td, "external.csv"), useBytes = TRUE)
    writeLines("analysis-code", file.path(td, "core.R"), useBytes = TRUE)
    writeLines("launcher-code", file.path(td, "run.ps1"), useBytes = TRUE)
  }
  input_a <- c(training = file.path(td1, "training.csv"),
               external = file.path(td1, "external.csv"))
  code_a <- c(core = file.path(td1, "core.R"),
              launcher = file.path(td1, "run.ps1"))
  input_b <- c(training = file.path(td2, "training.csv"),
               external = file.path(td2, "external.csv"))
  code_b <- c(core = file.path(td2, "core.R"),
              launcher = file.path(td2, "run.ps1"))
  params <- list(
    repeats = 1L, bootstrap = list(reps = 100L, seed = 11L),
    mapping = list(primary = "unique_mean", sensitivity = "unique_highest_mean"),
    force_nested = FALSE
  )
  key <- make_run_key(input_a, code_a, params, return_manifest = TRUE)
  assert_true(grepl("^[0-9a-f]{64}$", key$key), "run key is not SHA-256")
  assert_true(identical(
    key$key, make_run_key(input_b, code_b, params)
  ), "absolute relocation changed a content-addressed run key")
  assert_true("launcher" %in% names(key$manifest$code),
              "launcher is absent from the code manifest")
  assert_true("org.Hs.eg.db" %in% names(key$manifest$packages),
              "org.Hs.eg.db version is absent from the run identity")

  changed <- params
  changed$bootstrap$reps <- 101L
  assert_true(!identical(key$key, make_run_key(input_a, code_a, changed)),
              "bootstrap change did not change run key")
  changed <- params
  changed$mapping$primary <- "unique_highest_mean"
  assert_true(!identical(key$key, make_run_key(input_a, code_a, changed)),
              "mapping change did not change run key")
  changed <- params
  changed$force_nested <- TRUE
  assert_true(!identical(key$key, make_run_key(input_a, code_a, changed)),
              "execution parameter change did not change run key")

  writeLines("changed-input", input_a[["external"]], useBytes = TRUE)
  assert_true(!identical(key$key, make_run_key(input_a, code_a, params)),
              "input change did not change run key")
  writeLines("external-input", input_a[["external"]], useBytes = TRUE)
  writeLines("changed-launcher", code_a[["launcher"]], useBytes = TRUE)
  assert_true(!identical(key$key, make_run_key(input_a, code_a, params)),
              "launcher change did not change run key")
})

run_test("parameter manifest preserves all list names and round-trips exactly", {
  td <- tempfile("crc9-parameter-manifest-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE, force = TRUE), add = TRUE)
  path <- file.path(td, "run_parameters.txt")
  parameters <- list(
    output_contract_version = "crc9-artifact-contract-test",
    pipeline_version = "synthetic-pipeline-v1",
    nested_pipeline = list(
      outer_repeats = 1L,
      lambda_ratio_grid = c(0.5, 0.25),
      flags = c(training_only = TRUE, fail_closed = TRUE)
    ),
    probe_mapping = list(
      primary = "unique_mean",
      sensitivity = "unique_highest_mean"
    )
  )
  write_dput_manifest(parameters, path)
  restored <- dget(path)
  assert_true(identical(restored, parameters),
              "parameter manifest did not round-trip exactly")
  assert_true(identical(names(restored), names(parameters)),
              "top-level parameter names were lost")
  assert_true(identical(
    names(restored$nested_pipeline), names(parameters$nested_pipeline)
  ), "nested pipeline parameter names were lost")
  assert_true(identical(
    names(restored$nested_pipeline$flags),
    names(parameters$nested_pipeline$flags)
  ), "named-vector parameter names were lost")
})

run_test("code snapshot is complete and byte-identical to hashed sources", {
  td <- tempfile("crc9-code-snapshot-")
  source_dir <- file.path(td, "source")
  run_root <- file.path(td, "run")
  dir.create(source_dir, recursive = TRUE)
  dir.create(run_root)
  on.exit(unlink(td, recursive = TRUE, force = TRUE), add = TRUE)
  sources <- c(
    core_analysis = file.path(source_dir, "core.R"),
    local_launcher = file.path(source_dir, "run.ps1")
  )
  writeLines(c("x <- 1L", "stopifnot(x == 1L)"), sources[[1L]], useBytes = TRUE)
  writeLines("Write-Output 'synthetic launcher'", sources[[2L]], useBytes = TRUE)

  snapshot <- snapshot_code_files(sources, run_root)
  assert_true(nrow(snapshot) == length(sources),
              "code snapshot omitted one or more source files")
  assert_true(identical(snapshot$name, names(sources)),
              "code snapshot identity order changed")
  copied <- file.path(run_root, snapshot$relative_path)
  assert_true(all(file.exists(copied)), "a code snapshot file is missing")
  assert_true(identical(
    unname(vapply(sources, sha256_file, character(1))),
    unname(vapply(copied, sha256_file, character(1)))
  ), "snapshot hashes differ from source hashes")
  assert_true(identical(
    unname(file.info(copied)$size), as.numeric(snapshot$bytes)
  ), "snapshot byte counts differ from recorded values")
  expect_error(snapshot_code_files(sources, run_root), "already exists")
})

run_test("registry rejects missing required and unresolved optional artifacts", {
  td <- tempfile("crc9-registry-missing-")
  dir.create(td)
  dir.create(file.path(td, "results"))
  on.exit(unlink(td, recursive = TRUE, force = TRUE), add = TRUE)
  registry <- new_artifact_registry(
    td,
    make_spec("results/required.csv", "results/optional.csv")
  )
  expect_error(
    finalize_artifact_registry(
      registry, file.path(td, "results", "output_manifest.csv"),
      control_relative_paths = "results/output_manifest.csv"
    ),
    "Required artifacts"
  )
  writeLines("required", file.path(td, "results", "required.csv"), useBytes = TRUE)
  expect_error(
    finalize_artifact_registry(
      registry, file.path(td, "results", "output_manifest.csv"),
      control_relative_paths = "results/output_manifest.csv"
    ),
    "Optional artifacts"
  )
})

run_test("registry rejects files outside the allow-list", {
  td <- tempfile("crc9-registry-rogue-")
  dir.create(td)
  dir.create(file.path(td, "results"))
  on.exit(unlink(td, recursive = TRUE, force = TRUE), add = TRUE)
  registry <- new_artifact_registry(td, make_spec("results/required.csv"))
  writeLines("required", file.path(td, "results", "required.csv"), useBytes = TRUE)
  writeLines("stale", file.path(td, "results", "stale.csv"), useBytes = TRUE)
  expect_error(
    finalize_artifact_registry(
      registry, file.path(td, "results", "output_manifest.csv"),
      control_relative_paths = "results/output_manifest.csv"
    ),
    "Unregistered files"
  )
})

run_test("failed run keeps tmp status and never creates final directory", {
  runs <- tempfile("crc9-runs-failed-")
  dir.create(runs)
  on.exit(unlink(runs, recursive = TRUE, force = TRUE), add = TRUE)
  key <- paste(rep("a", 64L), collapse = "")
  tmp <- file.path(runs, paste0(key, ".tmp"))
  final <- file.path(runs, key)
  dir.create(tmp)
  dir.create(file.path(tmp, "results"))
  status_path <- file.path(tmp, "results", "run_status.csv")
  write_run_status_file(status_path, "running", key)
  registry <- new_artifact_registry(tmp, make_spec("results/required.csv"))
  failure <- tryCatch({
    finalize_artifact_registry(
      registry, file.path(tmp, "results", "output_manifest.csv"),
      control_relative_paths = c(
        "results/run_status.csv", "results/output_manifest.csv"
      )
    )
    NULL
  }, error = function(e) e)
  assert_true(inherits(failure, "error"), "missing required artifact did not fail")
  write_run_status_file(
    status_path, "failed_or_interrupted", key,
    detail = conditionMessage(failure)
  )
  status <- data.table::fread(status_path, data.table = FALSE)
  assert_true(identical(status$status, "failed_or_interrupted"),
              "failed status was not retained")
  assert_true(dir.exists(tmp), "failed tmp directory was removed")
  assert_true(!dir.exists(final), "failed run created a final directory")
})

run_test("successful run atomically renames and preserves the four-hash chain", {
  runs <- tempfile("crc9-runs-success-")
  dir.create(runs)
  on.exit(unlink(runs, recursive = TRUE, force = TRUE), add = TRUE)
  key <- paste(rep("b", 64L), collapse = "")
  nested_key <- paste(rep("c", 64L), collapse = "")
  tmp <- file.path(runs, paste0(key, ".tmp"))
  final <- file.path(runs, key)
  dir.create(tmp)
  dir.create(file.path(tmp, "results"))
  dir.create(file.path(tmp, "logs"))
  required <- c(
    "results/input_manifest.csv", "results/code_manifest.csv",
    "results/run_parameters.txt", "results/result.csv", "logs/run.log"
  )
  optional <- "results/optional.csv"
  registry <- new_artifact_registry(tmp, make_spec(required, optional))
  writeLines("input", file.path(tmp, required[1]), useBytes = TRUE)
  writeLines("code", file.path(tmp, required[2]), useBytes = TRUE)
  writeLines("parameters", file.path(tmp, required[3]), useBytes = TRUE)
  writeLines("result", file.path(tmp, required[4]), useBytes = TRUE)
  writeLines("log", file.path(tmp, required[5]), useBytes = TRUE)
  register_not_generated(registry, optional, "eligibility_not_met")

  input_hash <- sha256_file(file.path(tmp, required[1]))
  code_hash <- sha256_file(file.path(tmp, required[2]))
  parameter_hash <- sha256_file(file.path(tmp, required[3]))
  status_path <- file.path(tmp, "results", "run_status.csv")
  manifest_path <- file.path(tmp, "results", "output_manifest.csv")
  write_run_status_file(
    status_path, "running", key, nested_analysis_key = nested_key,
    input_manifest_sha256 = input_hash,
    parameter_manifest_sha256 = parameter_hash,
    code_manifest_sha256 = code_hash
  )
  manifest <- finalize_artifact_registry(
    registry, manifest_path,
    control_relative_paths = c(
      "results/run_status.csv", "results/output_manifest.csv"
    )
  )
  output_hash <- sha256_file(manifest_path)
  write_run_status_file(
    status_path, "complete", key, nested_analysis_key = nested_key,
    input_manifest_sha256 = input_hash,
    parameter_manifest_sha256 = parameter_hash,
    code_manifest_sha256 = code_hash,
    output_manifest_sha256 = output_hash
  )
  commit_run_directory(tmp, final, key)

  assert_true(!dir.exists(tmp) && dir.exists(final),
              "atomic rename did not produce only the final directory")
  status <- data.table::fread(
    file.path(final, "results", "run_status.csv"), data.table = FALSE
  )
  assert_true(identical(status$status, "complete"), "final status is not complete")
  assert_true(identical(status$input_manifest_sha256, input_hash),
              "input manifest hash chain is broken")
  assert_true(identical(status$parameter_manifest_sha256, parameter_hash),
              "parameter manifest hash chain is broken")
  assert_true(identical(status$code_manifest_sha256, code_hash),
              "code manifest hash chain is broken")
  assert_true(identical(
    status$output_manifest_sha256,
    sha256_file(file.path(final, "results", "output_manifest.csv"))
  ), "output manifest hash chain is broken")
  assert_true(all(manifest$status %in% c("generated", "not_generated")),
              "manifest contains unresolved rows")
  assert_true(
    manifest$status[manifest$artifact_id == optional] == "not_generated" &&
      nzchar(manifest$reason[manifest$artifact_id == optional]),
    "optional artifact lacks an explicit not_generated reason"
  )
  generated <- manifest[manifest$status == "generated", , drop = FALSE]
  recomputed <- vapply(
    file.path(final, generated$relative_path), sha256_file, character(1)
  )
  assert_true(identical(unname(recomputed), generated$sha256),
              "a generated artifact hash does not verify")
})

run_test("core and launcher retain the fail-closed isolation contract", {
  core <- readLines(core_file, warn = FALSE)
  launcher <- readLines(launcher_file, warn = FALSE)
  key_line <- grep("current_run_key <- run_identity\\$key", core)
  first_write <- grep("fwrite\\(|writeLines\\(|saveRDS\\(", core)[1L]
  assert_true(length(key_line) == 1L && key_line < first_write,
              "core writes an artifact before computing run_key")
  assert_true(any(grepl("local_launcher =", core, fixed = TRUE)),
              "launcher is not included in code_files")
  assert_true(any(grepl("snapshot_code_files(code_files", core, fixed = TRUE)),
              "core does not preserve an executable code snapshot")
  assert_true(any(grepl("write_dput_manifest(run_parameters", core, fixed = TRUE)),
              "core does not use the round-trip parameter serializer")
  assert_true(any(grepl('file.path("code_snapshot", basename(code_files))',
                        core, fixed = TRUE)),
              "code snapshots are absent from the artifact allow-list")
  assert_true(any(grepl(
    '"TCGA_COAD_signature_gene_mapping.csv"', core, fixed = TRUE
  )), "TCGA gene-identity audit is absent from the artifact contract")
  assert_true(any(grepl(
    "exact_SYMBOL_then_unique_reverse_unique_ENTREZ_alias", core,
    fixed = TRUE
  )), "TCGA gene-identity resolution rule is absent from run parameters")
  assert_true(any(grepl("commit_run_directory", core, fixed = TRUE)),
              "core lacks the atomic commit call")
  assert_true(!any(grepl("CRC_RUN_OUTPUT_ROOT", c(core, launcher), fixed = TRUE)),
              "legacy timestamp/PID output root is still active")
  complete_line <- grep('write_run_status\\(.*"complete"|"complete",', core)
  commit_line <- grep("commit_run_directory", core)
  assert_true(length(complete_line) >= 1L && tail(complete_line, 1L) < tail(commit_line, 1L),
              "complete status is not immediately upstream of commit")
})

report <- do.call(rbind, results)
cat("\nSummary:\n")
print(table(report$status))
if (any(report$status != "PASS")) {
  quit(save = "no", status = 1L)
}
quit(save = "no", status = 0L)
