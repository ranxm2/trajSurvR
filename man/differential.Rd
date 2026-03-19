\name{fit_differential_mixed_models}
\alias{fit_differential_mixed_models}
\title{Differential protein analysis via linear mixed-effects models}
\description{
Fits one linear mixed model per protein using \pkg{lme4} (or \pkg{lmerTest}
when available for Satterthwaite p-values) with a random intercept per
subject. Returns BH-adjusted FDR and fold-change direction.
}
\usage{
fit_differential_mixed_models(
  analysis_df,
  proteins,
  id_col = "UIDx",
  group_col = "Group",
  constant_covariates = c("Sex", "GenoGroup"),
  time_varying_covariates = c("CollAge"),
  case_level = "Affected",
  ref_level = "Control",
  min_rows = 20
)
}
\arguments{
  \item{analysis_df}{Data frame with subject IDs, group, covariates, and
    protein columns.}
  \item{proteins}{Character vector of protein column names to test.}
  \item{id_col}{Subject identifier column.}
  \item{group_col}{Group/condition column.}
  \item{constant_covariates}{Subject-level covariates (e.g. sex, genotype)
    that do not change between visits.}
  \item{time_varying_covariates}{Visit-level covariates that may change between
    visits (e.g. age, weight).}
  \item{case_level}{Factor level to compare against the reference.}
  \item{ref_level}{Reference factor level.}
  \item{min_rows}{Minimum complete rows required to fit a model.}
}
\value{
A tibble with columns \code{Protein}, \code{N}, \code{coef}, \code{se},
\code{t}, \code{p_value}, \code{fdr}, \code{direction}, and \code{significant},
sorted by ascending p-value.
}
\examples{
\dontrun{
dat <- simulate_trajsurv_data()
prots <- grep("^protein_", colnames(dat), value = TRUE)
res <- fit_differential_mixed_models(
  dat, proteins = prots,
  case_level = "Convert", ref_level = "Control"
)
head(res)
}
}
\keyword{models}
