#' Quality control and matrix generation for Olink-like data
#'
#' @param df Long-format data with sample ID, assay name, and value.
#' @param sample_qc_col Column name for sample QC status.
#' @param assay_qc_col Column name for assay QC status.
#' @param pass_value QC value that indicates pass.
#' @param sample_id_col Sample ID column.
#' @param assay_col Assay/protein column.
#' @param value_col Numeric value column (e.g. NPX).
#' @param missing_rate_threshold Maximum allowed missing-rate per protein.
#'
#' @return A list with QC-filtered long data, wide matrix, retained proteins, and summary.
qc_olink_data <- function(
  df,
  sample_qc_col = "SampleQC",
  assay_qc_col = "AssayQC",
  pass_value = "PASS",
  sample_id_col = "SampleID",
  assay_col = "Assay",
  value_col = "NPX",
  missing_rate_threshold = 0.40
) {
  stopifnot(is.data.frame(df))

  required_cols <- c(sample_id_col, assay_col, value_col, sample_qc_col, assay_qc_col)
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  n_raw_samples <- length(unique(df[[sample_id_col]]))
  n_raw_assays <- length(unique(df[[assay_col]]))

  df_qc <- df[df[[sample_qc_col]] == pass_value & df[[assay_qc_col]] == pass_value, , drop = FALSE]
  df_qc <- dplyr::as_tibble(df_qc)

  n_qc_samples <- length(unique(df_qc[[sample_id_col]]))
  n_qc_assays <- length(unique(df_qc[[assay_col]]))

  # Average duplicated sample-assay measurements.
  df_collapsed <- df_qc |>
    dplyr::group_by(.data[[sample_id_col]], .data[[assay_col]]) |>
    dplyr::summarise(value = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop")

  wide <- df_collapsed |>
    tidyr::pivot_wider(
      names_from = dplyr::all_of(assay_col),
      values_from = dplyr::all_of("value")
    ) |>
    dplyr::as_tibble()

  protein_cols <- setdiff(colnames(wide), sample_id_col)
  if (length(protein_cols) == 0) {
    stop("No protein columns after QC and widening.")
  }

  miss_rate <- vapply(
    protein_cols,
    function(p) mean(is.na(wide[[p]])),
    numeric(1)
  )
  keep_proteins <- protein_cols[miss_rate <= missing_rate_threshold]

  wide_keep <- wide[, c(sample_id_col, keep_proteins), drop = FALSE]
  summary_tbl <- dplyr::tibble(
    Stage = c("Raw", "After QC", "After missing-rate filter"),
    N_samples = c(n_raw_samples, n_qc_samples, nrow(wide_keep)),
    N_proteins = c(n_raw_assays, n_qc_assays, length(keep_proteins))
  )

  list(
    long_qc = df_qc,
    wide_matrix = wide_keep,
    protein_missing_rate = dplyr::tibble(Protein = protein_cols, MissingRate = miss_rate),
    keep_proteins = keep_proteins,
    summary = summary_tbl
  )
}


#' Return a compact QC summary table
#'
#' @param qc_result Output from \code{qc_olink_data()}.
#'
#' @return A data frame with stage-wise counts.
summarize_qc <- function(qc_result) {
  if (is.null(qc_result$summary)) {
    stop("Input must be result of qc_olink_data().")
  }
  qc_result$summary
}
