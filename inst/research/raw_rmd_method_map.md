# trajSurvR Method Mapping from `raw/*.Rmd`

This note documents how the original notebooks were consolidated into package functions.

## Source Rmd inventory and extracted methods

- `raw/1.Input_visits.Rmd`
  - Visit timeline processing and longitudinal subject-level ordering.
  - In package: used as data-shaping conventions (`SampleID`, `UIDx`, visit-level keys).

- `raw/1.Differentially_Expressed_Proteins.Rmd`
  - QC filters (`SampleQC == PASS`, `AssayQC == PASS`), low-read filtering.
  - Mixed-effect differential analysis with subject random intercept.
  - In package: `qc_olink_data()`, `fit_differential_mixed_models()`.

- `raw/4.Correlation_Study.Rmd`
  - Correlation matrices (Pearson/Spearman), complete-case protein subset handling.
  - In package: `compute_correlation_study()`.

- `raw/3.Protein_Trajectories_All.Rmd` and `raw/1.Protein_Trajectory_NEFL.Rmd`
  - GAM/GAMM trajectory fitting against `YrSinceOs`, confidence-band significance runs.
  - In package: `fit_trajectory_gam()`, `run_trajectory_batch()`.

- `raw/1.Cox_model_baseline.Rmd`
  - Baseline survival dataset construction from first visit.
  - Per-protein adjusted Cox and multivariable LASSO Cox.
  - In package: `build_baseline_survival_data()`, `fit_baseline_cox()`, `fit_lasso_cox()`.

- `raw/2.Time-depentd-cox.Rmd`
  - Counting-process interval construction and time-dependent Cox models.
  - In package: `make_counting_process_data()`, `fit_time_dependent_cox()`.

- `raw/2.Binary_Logisics_regression_local.Rmd`,
  `raw/3.Binary_Logisics_regression_data_driven.Rmd`,
  `raw/4.Panel_Evaluation.Rmd`
  - Logistic event prediction, CV AUC, panel comparison, forward greedy selection.
  - In package: `make_event_outcome()`, `cv_logistic_auc()`, `forward_greedy_selection()`,
    `fit_event_prediction_lasso()`.

- `raw/2.Phenoconversion_timing_estimation_LOOCV.Rmd`
  - LOOCV timing prediction and error metrics.
  - In package: `fit_onset_time_regression()`.

## Feature selection strategies implemented

- LASSO:
  - Logistic: `fit_event_prediction_lasso()`
  - Cox: `fit_lasso_cox()`
  - Generic API: `select_features_lasso()`

- Forward greedy search:
  - `forward_greedy_selection()`

## Covered method groups requested for trajSurvR

1. QC
2. Differential
3. Correlation study
4. Trajectories
5. Time-to-event analysis (multiple models)
6. Event prediction (logistic regression)
7. Time-of-onset prediction (regression with baseline)
8. Feature selection (LASSO + forward greedy)
