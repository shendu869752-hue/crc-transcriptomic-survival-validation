# Regression and certification tests

Run the portable suite from the repository root:

```powershell
.\run_tests.ps1
```

The suite covers strict repeated nested cross-validation, unique probe mapping, human SYMBOL/alias resolution, atomic run isolation, reporting contracts, the training-partition-only highest-mean sensitivity, integration checks, and completed-run validation.

To validate a completed formal run separately:

```powershell
Rscript reanalysis_v3/tests/validate_completed_run.R reanalysis_v3/runs/<run_key>
```

The certified manuscript run was:

```text
b2e747538453cf33b4321109514bdc2f06e38e688e14d1ee434c1677b4f6657c
```

Test-created temporary directories are excluded from the public release.
