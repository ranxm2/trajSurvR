\name{plot_differential_volcano}
\alias{plot_differential_volcano}
\alias{plot_correlation_heatmap}
\alias{plot_trajectory_fit}
\alias{plot_cox_forest}
\alias{plot_greedy_history}
\alias{plot_onset_predictions}
\title{Plotting helpers for trajSurvR results}
\description{
Convenience functions that produce \pkg{ggplot2} graphics from the outputs
of the main analysis functions.

\itemize{
  \item \code{plot_differential_volcano} — volcano plot of effect sizes vs
    -log10(p).
  \item \code{plot_correlation_heatmap} — heatmap of pairwise correlations.
  \item \code{plot_trajectory_fit} — observed data plus fitted GAM curve with
    95\% confidence ribbon.
  \item \code{plot_cox_forest} — forest-style dot plot of hazard ratios.
  \item \code{plot_greedy_history} — AUC trajectory during forward greedy
    selection.
  \item \code{plot_onset_predictions} — scatter of observed vs LOOCV-predicted
    onset times.
}
}
\usage{
plot_differential_volcano(
  diff_result,
  p_col = "p_value",
  coef_col = "coef",
  fdr_col = "fdr"
)

plot_correlation_heatmap(corr_result)

plot_trajectory_fit(trajectory_fit, time_col = "YrSinceOs")

plot_cox_forest(
  cox_result,
  hr_col = "HR",
  protein_col = "Protein",
  top_n = 20
)

plot_greedy_history(greedy_result)

plot_onset_predictions(onset_fit)
}
\arguments{
  \item{diff_result}{Output of \code{fit_differential_mixed_models()}.}
  \item{p_col}{P-value column name.}
  \item{coef_col}{Coefficient column name.}
  \item{fdr_col}{FDR column name for colour mapping.}
  \item{corr_result}{Output of \code{compute_correlation_study()}, or a
    numeric correlation matrix.}
  \item{trajectory_fit}{Output of \code{fit_trajectory_gam()}.}
  \item{time_col}{Time column in the predictions data frame.}
  \item{cox_result}{Output of \code{fit_baseline_cox()} or
    \code{fit_time_dependent_cox()}.}
  \item{hr_col}{Hazard-ratio column name.}
  \item{protein_col}{Feature label column name.}
  \item{top_n}{Maximum number of features to display.}
  \item{greedy_result}{Output of \code{forward_greedy_selection()}.}
  \item{onset_fit}{Output of \code{fit_onset_time_regression()}.}
}
\value{
All functions return a \code{ggplot} object.
}
\keyword{hplot}
