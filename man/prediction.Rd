\name{make_event_outcome}
\alias{make_event_outcome}
\alias{cv_logistic_auc}
\alias{fit_event_prediction_lasso}
\alias{fit_onset_time_regression}
\alias{select_features_lasso}
\alias{forward_greedy_selection}
\title{Event prediction and feature selection}
\description{
Functions for binary event prediction and onset-time estimation:

\itemize{
  \item \code{make_event_outcome} creates a binary \code{CASE} column based on
    a time horizon.
  \item \code{cv_logistic_auc} estimates cross-validated AUROC for a logistic
    regression model.
  \item \code{fit_event_prediction_lasso} fits a LASSO logistic model for event
    prediction.
  \item \code{fit_onset_time_regression} fits a linear regression with LOOCV
    for onset-time prediction.
  \item \code{select_features_lasso} is a general-purpose LASSO feature selector
    that works with \code{glmnet}.
  \item \code{forward_greedy_selection} performs forward greedy feature selection
    maximising a user-supplied evaluator (defaults to CV AUC).
}
}
\usage{
make_event_outcome(
  visit_df,
  horizon_years = 2,
  group_col = "Group",
  converter_label = "Convert",
  control_label = c("Control", "Potential Convert"),
  onset_col = "YrSinceOs"
)

cv_logistic_auc(
  data,
  feature_cols,
  outcome_col = "CASE",
  fold_col = NULL,
  n_folds = 5,
  seed = 2026
)

fit_event_prediction_lasso(
  data,
  outcome_col = "CASE",
  feature_cols,
  nfolds = 5
)

fit_onset_time_regression(
  data,
  outcome_col = "Y",
  feature_cols,
  constant_covariates = c("Sex", "GenoGroup"),
  time_varying_covariates = c("CollAge"),
  id_col = "UIDx"
)

select_features_lasso(
  x, y,
  family = c("binomial", "gaussian", "cox"),
  nfolds = 5,
  alpha = 1
)

forward_greedy_selection(
  data,
  candidate_features,
  base_features = character(0),
  evaluator = function(dat, feats) cv_logistic_auc(dat, feats),
  min_improvement = 1e-4,
  max_steps = 20
)
}
\arguments{
  \item{visit_df}{Visit-level data frame.}
  \item{horizon_years}{Prediction horizon in years.}
  \item{group_col}{Group column.}
  \item{converter_label}{Label for converters.}
  \item{control_label}{Label(s) for controls/non-converters.}
  \item{onset_col}{Years-since-onset column.}
  \item{data}{Modelling data frame.}
  \item{feature_cols}{Predictor column names.}
  \item{outcome_col}{Outcome column (binary for logistic, continuous for
    regression).}
  \item{fold_col}{Optional precomputed fold membership column.}
  \item{n_folds}{Number of cross-validation folds.}
  \item{seed}{Random seed.}
  \item{nfolds}{Number of LASSO cross-validation folds.}
  \item{constant_covariates}{Subject-level covariates always included.}
  \item{time_varying_covariates}{Visit-level covariates always included.}
  \item{id_col}{Subject ID for LOOCV.}
  \item{x}{Predictor matrix for \code{select_features_lasso}.}
  \item{y}{Response vector or \code{Surv} object.}
  \item{family}{glmnet family: \code{"binomial"}, \code{"gaussian"}, or
    \code{"cox"}.}
  \item{alpha}{Elastic-net mixing parameter (1 = LASSO).}
  \item{candidate_features}{Candidate variables for forward selection.}
  \item{base_features}{Variables always included in every model during
    selection.}
  \item{evaluator}{Function accepting \code{(data, features)} returning a
    scalar score.}
  \item{min_improvement}{Minimum score gain required to continue selection.}
  \item{max_steps}{Maximum number of selection steps.}
}
\value{
\code{make_event_outcome} returns the input data frame with an added
\code{CASE} column.

\code{cv_logistic_auc} returns a scalar AUROC.

\code{fit_event_prediction_lasso} and \code{select_features_lasso} return a
list with \code{cv_fit}, \code{selected_min}, and \code{selected_1se}.

\code{fit_onset_time_regression} returns a list with \code{model},
\code{loocv_predictions}, and \code{metrics}.

\code{forward_greedy_selection} returns a list with \code{selected},
\code{history}, and \code{final_score}.
}
\keyword{models}
