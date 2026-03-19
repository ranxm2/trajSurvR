\name{qc_olink_data}
\alias{qc_olink_data}
\alias{summarize_qc}
\title{Quality control for Olink-like long-format data}
\description{
\code{qc_olink_data} filters samples and assays by QC flags, collapses
duplicate measurements, widens the data to a protein matrix, and applies a
missing-rate filter. \code{summarize_qc} extracts the summary table from the
result.
}
\usage{
qc_olink_data(
  df,
  sample_qc_col = "SampleQC",
  assay_qc_col = "AssayQC",
  pass_value = "PASS",
  sample_id_col = "SampleID",
  assay_col = "Assay",
  value_col = "NPX",
  missing_rate_threshold = 0.40
)

summarize_qc(qc_result)
}
\arguments{
  \item{df}{Long-format data frame with one row per sample-assay measurement.}
  \item{sample_qc_col}{Column containing sample QC status.}
  \item{assay_qc_col}{Column containing assay QC status.}
  \item{pass_value}{Value indicating a passing QC flag.}
  \item{sample_id_col}{Sample identifier column.}
  \item{assay_col}{Assay/protein name column.}
  \item{value_col}{Numeric measurement column (e.g. NPX).}
  \item{missing_rate_threshold}{Maximum allowed per-protein missing rate.}
  \item{qc_result}{Output from \code{qc_olink_data()}.}
}
\value{
\code{qc_olink_data} returns a list with elements \code{long_qc},
\code{wide_matrix}, \code{protein_missing_rate}, \code{keep_proteins}, and
\code{summary}. \code{summarize_qc} returns the summary data frame.
}
\keyword{models}
