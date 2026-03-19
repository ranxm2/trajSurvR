# trajSurvR

**Longitudinal Trajectory and Survival Analysis for Proteomic Data**

`trajSurvR` wraps a complete longitudinal proteomics analysis workflow into
reusable R functions, covering quality control through trajectory modelling,
survival analysis, event prediction, and feature selection.

## Installation

```r
install.packages("remotes")
remotes::install_github("ranxm2/trajSurvR", build_vignettes = TRUE)
```

## Quick Start

```r
library(trajSurvR)

# Generate a synthetic dataset (180 subjects × 4 visits)
dat <- simulate_trajsurv_data(n_subjects = 180, n_visits = 4, seed = 2026)
head(dat)
```

Feature columns in demos use generic names (`protein_001`, `delta_001`, …).
Use `standardize_feature_names()` to map legacy marker names.

---

## Function Reference

### Data Simulation & Utilities

| Function | Description |
|---|---|
| `simulate_trajsurv_data()` | Generate a synthetic longitudinal cohort (subjects, visits, proteins, outcomes) |
| `generate_example_data()` | Low-level generator underlying `simulate_trajsurv_data()` |
| `load_example_data()` | Load the bundled example CSV from the package |
| `standardize_feature_names()` | Rename legacy marker columns to generic IDs (`protein_001`, `delta_001`, …) |

---

### Quality Control

| Function | Description |
|---|---|
| `qc_olink_data()` | Filter by sample/assay QC flags, collapse duplicates, widen to protein matrix, apply missing-rate filter |
| `summarize_qc()` | Extract stage-wise sample and protein counts from a QC result |

---

### Differential Expression

| Function | Description |
|---|---|
| `fit_differential_mixed_models()` | Per-protein linear mixed model (random subject intercept); supports constant and time-varying covariates; returns BH-adjusted FDR |

---

### Correlation Study

| Function | Description |
|---|---|
| `compute_correlation_study()` | Pairwise Pearson/Spearman/Kendall correlation matrix and long-format pair table |

---

### Trajectory Modelling

| Function | Description |
|---|---|
| `fit_trajectory_gam()` | GAM/GAMM with smooth time term and subject random intercept for a single protein delta |
| `run_trajectory_batch()` | Iterate `fit_trajectory_gam()` over multiple proteins; returns summary table and combined predictions |

---

### Time-to-Event Analysis

| Function | Description |
|---|---|
| `build_baseline_survival_data()` | Build per-subject baseline dataset with event status and survival time |
| `fit_baseline_cox()` | Univariate baseline Cox model per protein |
| `fit_lasso_cox()` | Penalised LASSO Cox via `cv.glmnet` |
| `make_counting_process_data()` | Convert visit table to counting-process interval format for time-dependent covariates |
| `fit_time_dependent_cox()` | Univariate time-dependent Cox model per protein |

---

### Event Prediction

| Function | Description |
|---|---|
| `make_event_outcome()` | Create binary `CASE` label for onset-within-horizon prediction |
| `cv_logistic_auc()` | Cross-validated AUROC for a logistic regression model |
| `fit_event_prediction_lasso()` | LASSO logistic regression for event prediction |

---

### Onset-Time Prediction

| Function | Description |
|---|---|
| `fit_onset_time_regression()` | Linear regression with LOOCV to predict continuous time-to-onset; returns MAE, RMSE, correlation |

---

### Feature Selection

| Function | Description |
|---|---|
| `select_features_lasso()` | General-purpose LASSO selector for `binomial`, `gaussian`, or `cox` families |
| `forward_greedy_selection()` | Forward greedy selection maximising any user-supplied scoring function (default: CV AUC) |

---

### Plotting

| Function | Description |
|---|---|
| `plot_differential_volcano()` | Volcano plot of effect sizes vs −log10(p) |
| `plot_correlation_heatmap()` | Heatmap of pairwise correlation matrix |
| `plot_trajectory_fit()` | Observed data + GAM smooth curve with 95% CI ribbon |
| `plot_cox_forest()` | Forest-style dot plot of Cox hazard ratios |
| `plot_greedy_history()` | AUC gain at each step of forward greedy selection |
| `plot_onset_predictions()` | Observed vs LOOCV-predicted onset time scatter plot |

---

## Vignette

A full end-to-end workflow demo (PDF) is available in
[`docs/trajSurvR_demo.pdf`](docs/trajSurvR_demo.pdf), or run:

```r
vignette("trajSurvR_demo", package = "trajSurvR")
```

## License

MIT © Ximing Ran
