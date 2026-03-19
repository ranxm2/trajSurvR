\name{compute_correlation_study}
\alias{compute_correlation_study}
\title{Correlation study for protein features}
\description{
Computes a pairwise correlation matrix and returns both the matrix and a
long-format pair table.
}
\usage{
compute_correlation_study(mat_wide, proteins = NULL, method = "pearson")
}
\arguments{
  \item{mat_wide}{Wide data frame or matrix with one row per sample.}
  \item{proteins}{Optional character vector of column names to include.
    Defaults to all numeric columns.}
  \item{method}{Correlation method: \code{"pearson"}, \code{"spearman"}, or
    \code{"kendall"}.}
}
\value{
A list with elements \code{correlation_matrix} (square numeric matrix) and
\code{pair_table} (tibble with columns \code{Protein1}, \code{Protein2},
\code{Correlation}).
}
\keyword{models}
