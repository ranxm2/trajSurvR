\name{build_baseline_survival_data}
\alias{build_baseline_survival_data}
\alias{fit_baseline_cox}
\alias{fit_lasso_cox}
\alias{make_counting_process_data}
\alias{fit_time_dependent_cox}
\title{Time-to-event analysis functions}
\description{
A family of functions for Cox proportional-hazards modelling:

\itemize{
  \item \code{build_baseline_survival_data} constructs a per-subject baseline
    dataset with event status and survival time.
  \item \code{fit_baseline_cox} fits univariate baseline Cox models per protein.
  \item \code{fit_lasso_cox} fits a penalised Cox model via \code{cv.glmnet}.
  \item \code{make_counting_process_data} converts a visit table to
    counting-process (interval) format for time-dependent covariates.
  \item \code{fit_time_dependent_cox} fits univariate time-dependent Cox models
    per protein.
}
}
\usage{
build_baseline_survival_data(
  visit_info,
  protein_wide,
  sample_id_col = "SampleID",
  id_col = "UIDx",
  group_col = "Group",
  converter_label = "Convert",
  pre_converter_label = "Potential Convert",
  onset_col = "YrSinceOs",
  censor_col = "YrSinceCen",
  constant_covariates = c("Sex", "GenoGroup"),
  time_varying_covariates = c("CollAge"),
  visit_order_col = "CollNum"
)

fit_baseline_cox(
  dat,
  proteins,
  constant_covariates = c("Sex", "GenoGroup"),
  time_varying_covariates = c("CollAge"),
  time_col = "time_yr",
  event_col = "event",
  min_rows = 20,
  min_events = 5
)

fit_lasso_cox(
  dat,
  proteins,
  constant_covariates = c("Sex", "GenoGroup"),
  time_varying_covariates = c("CollAge"),
  time_col = "time_yr",
  event_col = "event",
  nfolds = 5
)

make_counting_process_data(
  visit_df,
  id_col = "UIDx",
  visit_time_col = "time_since_baseline",
  stop_time_col = "stop_time",
  event_col = "event",
  keep_visit_cols = character(0)
)

fit_time_dependent_cox(
  td_data,
  proteins,
  constant_covariates = c("Sex", "GenoGroup"),
  time_varying_covariates = c("CollAge"),
  id_col = "UIDx",
  tstart_col = "tstart",
  tstop_col = "tstop",
  event_col = "event",
  min_subjects = 10,
  min_events = 5
)
}
\arguments{
  \item{visit_info}{Visit-level data frame.}
  \item{protein_wide}{Wide protein matrix keyed by sample ID.}
  \item{sample_id_col}{Sample identifier column.}
  \item{id_col}{Subject identifier column.}
  \item{group_col}{Group/condition column.}
  \item{converter_label}{Label for converters/cases.}
  \item{pre_converter_label}{Label for pre-converters or censored subjects.}
  \item{onset_col}{Years-since-onset column.}
  \item{censor_col}{Years-since-censoring column.}
  \item{constant_covariates}{Subject-level covariates.}
  \item{time_varying_covariates}{Visit-level covariates.}
  \item{visit_order_col}{Column used to order visits for baseline selection.}
  \item{dat}{Survival data frame with event/time and predictors.}
  \item{proteins}{Character vector of protein column names.}
  \item{time_col}{Survival time column.}
  \item{event_col}{Event indicator column (0/1).}
  \item{min_rows}{Minimum complete rows for model fitting.}
  \item{min_events}{Minimum events for model fitting.}
  \item{nfolds}{Number of cross-validation folds for LASSO.}
  \item{visit_df}{Visit table for counting-process conversion.}
  \item{visit_time_col}{Time-since-baseline column.}
  \item{stop_time_col}{Subject-level follow-up stop time.}
  \item{keep_visit_cols}{Additional visit-level columns to carry into each
    interval row.}
  \item{td_data}{Counting-process data from \code{make_counting_process_data}.}
  \item{tstart_col}{Interval start column.}
  \item{tstop_col}{Interval stop column.}
  \item{min_subjects}{Minimum subjects required.}
}
\value{
\code{build_baseline_survival_data} returns a data frame with columns
\code{event} and \code{time_yr} plus joined proteins and covariates.

\code{fit_baseline_cox} and \code{fit_time_dependent_cox} return tibbles
with per-protein Cox statistics (\code{coef}, \code{HR}, \code{p_value},
\code{fdr}, \code{direction}).

\code{fit_lasso_cox} returns a list with \code{cv_fit}, \code{lambda_min},
\code{lambda_1se}, \code{selected_min}, \code{selected_1se}.

\code{make_counting_process_data} returns an interval-level data frame with
columns \code{tstart}, \code{tstop}, and \code{event}.
}
\keyword{survival}
