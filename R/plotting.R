#' Plot differential-analysis volcano
#'
#' @param diff_result Output of \code{fit_differential_mixed_models()}.
#' @param p_col P-value column name.
#' @param coef_col Coefficient column name.
#' @param fdr_col FDR column name; used to colour significant points.
#'
#' @return A \code{ggplot} object.
#' @export
plot_differential_volcano <- function(
  diff_result,
  p_col = "p_value",
  coef_col = "coef",
  fdr_col = "fdr"
) {
  req <- c(p_col, coef_col)
  miss <- setdiff(req, colnames(diff_result))
  if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

  df <- dplyr::as_tibble(diff_result) |>
    dplyr::filter(!is.na(.data[[p_col]]), !is.na(.data[[coef_col]])) |>
    dplyr::mutate(
      neglog10p = -log10(.data[[p_col]]),
      sig = if (fdr_col %in% colnames(diff_result)) .data[[fdr_col]] < 0.05 else .data[[p_col]] < 0.05
    )

  ggplot2::ggplot(df, ggplot2::aes(x = .data[[coef_col]], y = .data$neglog10p, color = .data$sig)) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::scale_color_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey50")) +
    ggplot2::labs(x = "Effect size", y = "-log10(p)", color = "Significant") +
    ggplot2::theme_bw()
}


#' Plot correlation matrix heatmap
#'
#' @param corr_result Result from \code{compute_correlation_study()}, or a
#'   numeric correlation matrix directly.
#'
#' @return A \code{ggplot} object.
#' @export
plot_correlation_heatmap <- function(corr_result) {
  cor_mat <- if (is.list(corr_result) && !is.null(corr_result$correlation_matrix)) {
    corr_result$correlation_matrix
  } else {
    corr_result
  }
  cor_long <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
  colnames(cor_long) <- c("Protein1", "Protein2", "Correlation")

  ggplot2::ggplot(cor_long, ggplot2::aes(x = .data$Protein1, y = .data$Protein2, fill = .data$Correlation)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1)) +
    ggplot2::labs(x = NULL, y = NULL, fill = "r") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}


#' Plot one trajectory fit
#'
#' @param trajectory_fit Result from \code{fit_trajectory_gam()}.
#' @param time_col Time column name in the predictions data frame.
#'
#' @return A \code{ggplot} object.
#' @export
plot_trajectory_fit <- function(trajectory_fit, time_col = "YrSinceOs") {
  if (is.null(trajectory_fit) || is.null(trajectory_fit$predictions) || is.null(trajectory_fit$data)) {
    stop("trajectory_fit must be the output of fit_trajectory_gam().")
  }

  dat  <- trajectory_fit$data
  pred <- trajectory_fit$predictions

  ggplot2::ggplot() +
    ggplot2::geom_point(data = dat, ggplot2::aes(x = .data[[time_col]], y = .data$delta), alpha = 0.5, size = 1.2) +
    ggplot2::geom_ribbon(
      data = pred,
      ggplot2::aes(x = .data[[time_col]], ymin = .data$lower, ymax = .data$upper),
      alpha = 0.2, fill = "grey40"
    ) +
    ggplot2::geom_line(data = pred, ggplot2::aes(x = .data[[time_col]], y = .data$fit), linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::labs(x = time_col, y = "Delta", title = "Trajectory fit") +
    ggplot2::theme_bw()
}


#' Forest-style dot plot of Cox model hazard ratios
#'
#' @param cox_result Output from \code{fit_baseline_cox()} or
#'   \code{fit_time_dependent_cox()}.
#' @param hr_col Hazard-ratio column name.
#' @param protein_col Feature/protein label column name.
#' @param top_n Maximum number of features to display.
#'
#' @return A \code{ggplot} object.
#' @export
plot_cox_forest <- function(
  cox_result,
  hr_col      = "HR",
  protein_col = "Protein",
  top_n       = 20
) {
  req  <- c(hr_col, protein_col)
  miss <- setdiff(req, colnames(cox_result))
  if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

  df <- dplyr::as_tibble(cox_result) |>
    dplyr::filter(!is.na(.data[[hr_col]]), is.finite(.data[[hr_col]]))
  if ("p_value" %in% colnames(df)) {
    df <- dplyr::arrange(df, .data$p_value)
  }
  df <- utils::head(df, top_n)
  df[[protein_col]] <- factor(df[[protein_col]], levels = rev(df[[protein_col]]))

  ggplot2::ggplot(df, ggplot2::aes(x = .data[[hr_col]], y = .data[[protein_col]])) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "Hazard ratio (log scale)", y = NULL, title = "Cox effect summary") +
    ggplot2::theme_bw()
}


#' Plot forward-greedy feature-selection history
#'
#' @param greedy_result Output from \code{forward_greedy_selection()}.
#'
#' @return A \code{ggplot} object.
#' @export
plot_greedy_history <- function(greedy_result) {
  if (is.null(greedy_result$history) || nrow(greedy_result$history) == 0) {
    stop("greedy_result has empty history.")
  }

  df <- greedy_result$history
  ggplot2::ggplot(df, ggplot2::aes(x = .data$Step, y = .data$Score)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::geom_text(ggplot2::aes(label = .data$Added), nudge_y = 0.01, size = 3, check_overlap = TRUE) +
    ggplot2::labs(x = "Step", y = "Score (AUC)", title = "Forward greedy selection path") +
    ggplot2::theme_bw()
}


#' Plot observed vs predicted onset time
#'
#' @param onset_fit Output from \code{fit_onset_time_regression()}.
#'
#' @return A \code{ggplot} object.
#' @export
plot_onset_predictions <- function(onset_fit) {
  if (is.null(onset_fit$loocv_predictions)) {
    stop("onset_fit must contain a loocv_predictions element.")
  }

  df <- onset_fit$loocv_predictions |>
    dplyr::filter(is.finite(.data$y_true), is.finite(.data$y_pred))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$y_pred, y = .data$y_true)) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::labs(x = "Predicted onset time", y = "Observed onset time",
                  title = "Onset prediction performance") +
    ggplot2::theme_bw()
}
