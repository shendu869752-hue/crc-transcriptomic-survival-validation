$ErrorActionPreference = "Stop"
$repoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$testRoot = Join-Path $repoRoot "reanalysis_v3/tests"
$rscriptCommand = Get-Command Rscript.exe -ErrorAction SilentlyContinue
if ($null -eq $rscriptCommand) {
    $rscriptCommand = Get-Command Rscript -ErrorAction Stop
}
$rscript = $rscriptCommand.Source

$tests = @(
    "run_nested_cv_tests.R",
    "test_unique_probe_mapping.R",
    "test_human_symbol_resolution.R",
    "test_run_isolation.R",
    "test_reporting_contract.R",
    "test_highest_mean_nested.R",
    "test_highest_mean_integration.R",
    "test_validate_completed_run.R"
)

foreach ($test in $tests) {
    Write-Output "Running $test"
    Push-Location $repoRoot
    try {
        # Keep the Rscript --file argument relative. On Windows, R 4.6.1 can
        # mis-handle non-ASCII components in an absolute repository path when
        # the process locale falls back from C.UTF-8.
        & $rscript (Join-Path "reanalysis_v3/tests" $test)
        if ($LASTEXITCODE -ne 0) { throw "$test failed with exit code $LASTEXITCODE" }
    }
    finally {
        Pop-Location
    }
}

Write-Output "All portable regression tests passed."
