#' Correlation study for selected proteins
#'
#' @param mat_wide Wide matrix/data frame with one row per sample.
#' @param proteins Optional vector of proteins. Default uses all numeric columns.
#' @param method Correlation method ("pearson", "spearman", "kendall").
#'
#' @return A list with correlation matrix and long-form pair table.
compute_correlation_study <- function(mat_wide, proteins = NULL, method = "pearson") {
  stopifnot(is.data.frame(mat_wide) || is.matrix(mat_wide))
  mat_wide <- as.data.frame(mat_wide)

  if (is.null(proteins)) {
    is_num <- vapply(mat_wide, is.numeric, logical(1))
    proteins <- names(is_num)[is_num]
  }
  proteins <- intersect(proteins, colnames(mat_wide))
  if (length(proteins) < 2) stop("Need at least 2 proteins for correlation study.")

  x <- mat_wide[, proteins, drop = FALSE]
  cor_mat <- stats::cor(x, method = method, use = "pairwise.complete.obs")

  cor_long <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
  colnames(cor_long) <- c("Protein1", "Protein2", "Correlation")
  cor_long <- dplyr::as_tibble(cor_long) |>
    dplyr::filter(Protein1 != Protein2)

  list(correlation_matrix = cor_mat, pair_table = cor_long)
}
