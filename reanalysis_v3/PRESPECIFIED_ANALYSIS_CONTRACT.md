# CRC public-data survival analysis: v3 locked analysis contract

Locked on 2026-09-17 before accepting any v3 real-data result.

## Scientific objective

Quantify optimism, feature-selection stability, and cross-cohort transportability of a data-derived colorectal-cancer prognostic score. The study does not claim a locked deployable clinical model, external calibration, a clinical cutoff, or individual-patient prediction.

## Development and validation roles

- GSE39582 is the sole development cohort.
- TCGA-COAD, GSE14333, GSE17536, and GSE17537 do not define the development feature space or tune the development algorithm.
- External effects are cohort-standardized prognostic associations because score preprocessing is cohort-adaptive; they are not absolute-risk calibration estimates.

## Primary GPL570 mapping

`unique_mean` is primary. Retain probes mapping to exactly one non-empty gene SYMBOL in the frozen `hgu133plus2.db` annotation, then calculate the arithmetic mean of all eligible log2-scale probes for each gene. This transformation is annotation-only and fixed across samples, folds, and cohorts.

## Prespecified mapping sensitivity

`unique_highest_mean` is a complete-pipeline sensitivity, not a replacement primary analysis.

- In every inner or outer training partition, select the eligible probe with the highest mean expression for each gene; resolve exact ties lexicographically by PROBEID.
- Apply only those training-selected probe identities to the paired validation partition.
- Select the external-use probe identities using the complete GSE39582 development cohort, then apply those identities unchanged to external GPL570 cohorts.
- Report probe-selection frequencies and mapping availability. Never reselect probes using an external cohort.
- The historical first-SYMBOL/per-sample-maximum rule is excluded from scientific sensitivity analyses.

## Internal validation

- Repeated strict nested cross-validation: 10 repeats, 5 outer folds, and 5 inner folds.
- Feature filtering, univariable screening, fallbacks, scaling, penalized tuning, coefficient thresholding, and Cox refitting are learned only from the relevant training partition.
- The relative lambda grid is fixed in advance; the one-standard-error rule selects the largest admissible penalty.
- Every patient receives exactly one out-of-fold prediction per repeat. Failed folds stop the run and remain auditable.
- Primary internal performance reporting is the distribution of repeat-level patient OOF Harrell C-indices. Fold quantiles are not called confidence intervals.

## External evaluation and uncertainty

- Primary external effect: Cox hazard ratio per one within-cohort SD of the continuous score, with 95% confidence interval.
- Discrimination: Harrell C-index with explicitly conditional patient-level percentile-bootstrap intervals for precomputed scores.
- Time-dependent AUC is reported only when the prespecified event and at-risk thresholds are met; non-estimable results and reasons are retained.
- Dichotomized Kaplan-Meier plots are descriptive visualizations and do not establish a clinical cutoff.

## Meta-analysis

- Primary OS synthesis: REML random effects with Knapp-Hartung/t inference.
- Report the 95% prediction interval, Q, Q-test P value, tau-squared, I-squared, degrees of freedom, and exact cohort inputs.
- Normal/Wald inference is a sensitivity analysis.
- Leave-one-out analyses retain the same primary Knapp-Hartung/t inference, with Wald results shown separately.

## Diagnostics and clinical adjustment

- Assess proportional hazards for score terms and preserve all results. Any prespecified response to a material violation must be reported rather than selected by significance.
- Assess score functional form using an explicitly named Cox model comparison; describe it as exploratory where event counts are limited.
- Clinical adjustment requires an explicit cohort-eligibility table. Non-eligible cohorts receive a recorded reason rather than being silently skipped.
- Within-stage-stratum C-indices are apparent descriptive quantities and are labeled accordingly.

## Integrity rules

- Do not reuse v2 numerical results in v3.
- Do not change methods in response to favorable or unfavorable v3 findings.
- Preserve null, inconsistent, and non-estimable findings.
- Each accepted run must be content-keyed, isolated, hash-manifested, complete, and reproducible from recorded inputs, code, parameters, R, and package versions.
- Pilot runs are smoke tests only and cannot supply manuscript numbers.

