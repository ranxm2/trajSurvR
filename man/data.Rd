\name{simulate_trajsurv_data}
\alias{simulate_trajsurv_data}
\alias{generate_example_data}
\alias{load_example_data}
\alias{standardize_feature_names}
\title{Simulate and load example datasets}
\description{
Functions for generating synthetic longitudinal proteomic datasets
and for standardising feature column names to generic identifiers.
}
\usage{
simulate_trajsurv_data(
  n_subjects = 180,
  n_visits = 4,
  seed = 2026,
  generic_feature_names = TRUE
)

generate_example_data(
  n_subjects = 180,
  n_visits = 4,
  seed = 2026,
  generic_feature_names = TRUE
)

load_example_data()

standardize_feature_names(
  data,
  mapping = c(
    NPX_NEFL = "protein_001",
    NPX_GPNMB = "protein_002",
    NPX_CXCL13 = "protein_003",
    delta_NEFL = "delta_001",
    delta_GPNMB = "delta_002",
    delta_CXCL13 = "delta_003"
  )
)
}
\arguments{
  \item{n_subjects}{Number of synthetic subjects (default 180).}
  \item{n_visits}{Visits per subject (default 4).}
  \item{seed}{Random seed for reproducibility.}
  \item{generic_feature_names}{If \code{TRUE}, feature columns are renamed to
    generic identifiers such as \code{protein_001}, \code{delta_001}.}
  \item{data}{Input data frame whose columns are to be renamed.}
  \item{mapping}{Named character vector: names are old column names, values are
    the replacement generic names.}
}
\value{
\code{simulate_trajsurv_data} and \code{generate_example_data} return a
longitudinal visit-level data frame. \code{load_example_data} reads a CSV
bundled in the package. \code{standardize_feature_names} returns the input
data frame with selected columns renamed.
}
\examples{
dat <- simulate_trajsurv_data(n_subjects = 40, n_visits = 3, seed = 1)
head(dat)
}
\keyword{datasets}
