\name{fit_trajectory_gam}
\alias{fit_trajectory_gam}
\alias{run_trajectory_batch}
\title{GAM/GAMM trajectory modelling}
\description{
\code{fit_trajectory_gam} fits a GAM or GAMM (with subject random intercept)
for a single protein delta over time. \code{run_trajectory_batch} iterates
over multiple proteins and returns a summary table plus combined predictions.
}
\usage{
fit_trajectory_gam(
  delta_df,
  time_col = "YrSinceOs",
  delta_col = "delta",
  id_col = "UIDx",
  covariates = character(0),
  k = 10,
  min_obs = 10,
  n_grid = 300
)

run_trajectory_batch(
  delta_wide,
  visit_design,
  sample_id_col = "SampleID",
  proteins = NULL,
  id_col = "UIDx",
  time_col = "YrSinceOs",
  constant_covariates = character(0),
  time_varying_covariates = character(0),
  min_obs = 10
)
}
\arguments{
  \item{delta_df}{Data frame with time, delta response, and subject ID.}
  \item{time_col}{Time-to-onset column.}
  \item{delta_col}{Delta (response) column.}
  \item{id_col}{Subject ID column.}
  \item{covariates}{Additional adjustment covariates.}
  \item{k}{Basis dimension for the smooth term.}
  \item{min_obs}{Minimum observations required to attempt fitting.}
  \item{n_grid}{Number of grid points for the prediction curve.}
  \item{delta_wide}{Wide delta matrix with sample ID and protein columns.}
  \item{visit_design}{Visit-level design data with sample ID, time, and
    subject ID.}
  \item{sample_id_col}{Sample identifier column.}
  \item{proteins}{Protein columns to fit. Defaults to all non-ID columns of
    \code{delta_wide}.}
  \item{constant_covariates}{Subject-level covariates (do not change between
    visits).}
  \item{time_varying_covariates}{Visit-level covariates.}
}
\value{
\code{fit_trajectory_gam} returns a list with elements \code{model},
\code{gam}, \code{data}, \code{predictions}, \code{p_value}, and
\code{intervals}, or \code{NULL} if fitting failed. \code{run_trajectory_batch}
returns a list with \code{models}, \code{summary}, and \code{predictions}.
}
\keyword{models}
