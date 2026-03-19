# trajSurvR

`trajSurvR` wraps the analysis logic from `raw/*.Rmd` into reusable R functions for:

- QC
- Differential analysis
- Correlation study
- Trajectory modeling
- Time-to-event analysis (baseline Cox, LASSO Cox, time-dependent Cox)
- Event prediction (logistic regression)
- Time-of-onset prediction (regression with baseline)
- Feature selection (LASSO and forward greedy search)

## Installation

Install from GitHub:

```r
install.packages("remotes")
remotes::install_github("ranxm2/trajSurvR", build_vignettes = TRUE)
```

## Quick start

```r
library(trajSurvR)
dat <- simulate_trajsurv_data(n_subjects = 180, n_visits = 4, seed = 2026, generic_feature_names = TRUE)
head(dat)
```

Feature names in demos are intentionally generic (`protein_001`, `delta_001`, etc.).
If you start from legacy marker names, use `standardize_feature_names()`.

See `vignettes/trajSurvR_demo.Rmd` for a full workflow demo.

