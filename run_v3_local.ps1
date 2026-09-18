param(
    [ValidateRange(1, 100)]
    [int]$Repeats = 10,
    [ValidateRange(100, 100000)]
    [int]$BootstrapReps = 2000,
    [switch]$ForceNested,
    [switch]$ForceDownload
)

$ErrorActionPreference = "Stop"
$repoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$v3Root = Join-Path $repoRoot "reanalysis_v3"
$rscriptCommand = Get-Command Rscript.exe -ErrorAction SilentlyContinue
if ($null -eq $rscriptCommand) {
    $rscriptCommand = Get-Command Rscript -ErrorAction Stop
}
$rscript = $rscriptCommand.Source

$externalDir = Join-Path $repoRoot "data/external_GEO"
$rawDir = Join-Path $v3Root "raw_inputs"
$preparedDir = Join-Path $v3Root "prepared_inputs"
$provenanceDir = Join-Path $v3Root "provenance"
$runsRoot = Join-Path $v3Root "runs"

& $rscript (Join-Path $v3Root "scripts/00_download_external_geo.R") `
    "--output-dir=$externalDir" "--force=$($ForceDownload.IsPresent.ToString().ToLowerInvariant())"
if ($LASTEXITCODE -ne 0) { throw "External GEO download failed" }

& $rscript (Join-Path $v3Root "scripts/00_acquire_prepare_inputs.R") `
    "--action=write-a-level" "--source-root=$repoRoot" "--raw-dir=$rawDir" `
    "--output-dir=$preparedDir" "--provenance-dir=$provenanceDir" `
    "--download=true" "--force=$($ForceDownload.IsPresent.ToString().ToLowerInvariant())"
if ($LASTEXITCODE -ne 0) { throw "Input reconstruction failed" }

$env:CRC_SOURCE_ROOT = $repoRoot
$env:CRC_WORK_ROOT = $repoRoot
$env:CRC_NESTED_REPEATS = [string]$Repeats
$env:CRC_BOOTSTRAP_REPS = [string]$BootstrapReps
$env:CRC_FORCE_NESTED = if ($ForceNested) { "1" } else { "0" }
$env:CRC_PREPARED_INPUT_DIR = $preparedDir
$env:CRC_RUNS_ROOT = $runsRoot

& $rscript (Join-Path $v3Root "scripts/01_core_reanalysis_v3.R")
if ($LASTEXITCODE -ne 0) { throw "CRC v3 analysis failed with exit code $LASTEXITCODE" }

Write-Output "Analysis completed. Validate the generated run before accepting results."

