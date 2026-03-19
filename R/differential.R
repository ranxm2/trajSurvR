#' Differential protein analysis via linear mixed-effects models
#'
#' Per-protein linear mixed model with a random intercept per subject.
#' Uses \pkg{lmerTest} for Satterthwaite p-values when installed; falls
#' back to \pkg{lme4} with Wald approximation otherwise.
#'
#' @param analysis_df Data frame containing subject ID, group, covariates,
#'   and protein columns.
#' @param proteins Character vector of protein column names to test.
#' @param id_col Subject ID column.
#' @param group_col Group column (e.g. \code{"Group"}).
#' @param constant_covariates Subject-level covariates that do not change
#'   between visits (e.g. sex, genotype).
#' @param time_varying_covariates Visit-level covariates that may change
#'   between visits (e.g. age, weight, height).
#' @param case_level Factor level to compare against the reference.
#' @param ref_level Reference factor level.
#' @param min_rows Minimum complete rows required to fit a model.
#'
#' @return A \code{tibble} with one row per protein: coefficient, standard
#'   error, t-statistic, p-value, BH-adjusted FDR, and direction.
#' @export
fit_differential_mixed_models <- function(
  analysis_df,
  proteins,
  id_col = "UIDx",
  group_col = "Group",
  constant_covariates = c("Sex", "GenoGroup"),
  time_varying_covariates = c("CollAge"),
  case_level = "Affected",
  ref_level = "Control",
  min_rows = 20
) {
  stopifnot(is.data.frame(analysis_df))
  stopifnot(is.character(proteins), length(proteins) > 0)

  use_lmerTest <- requireNamespace("lmerTest", quietly = TRUE)
  model_covariates <- unique(c(constant_covariates, time_varying_covariates))
  req <- c(id_col, group_col, model_covariates)
  miss <- setdiff(req, colnames(analysis_df))
  if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))

  dat0 <- dplyr::as_tibble(analysis_df)
  dat0[[group_col]] <- factor(dat0[[group_col]])
  if (all(c(ref_level, case_level) %in% levels(dat0[[group_col]]))) {
    dat0[[group_col]] <- stats::relevel(dat0[[group_col]], ref = ref_level)
  }

  out <- vector("list", length(proteins))

  for (i in seq_along(proteins)) {
    p <- proteins[i]
    if (!p %in% colnames(dat0)) {
      out[[i]] <- dplyr::tibble(Protein = p, N = NA_integer_, coef = NA_real_,
        se = NA_real_, t = NA_real_, p_value = NA_real_, Status = "Protein missing")
      next
    }

    dat <- dat0 |>
      dplyr::select(dplyr::all_of(c(id_col, group_col, model_covariates, p))) |>
      tidyr::drop_na()

    if (nrow(dat) < min_rows || dplyr::n_distinct(dat[[id_col]]) < 8) {
      out[[i]] <- dplyr::tibble(Protein = p, N = nrow(dat), coef = NA_real_,
        se = NA_real_, t = NA_real_, p_value = NA_real_, Status = "Insufficient data")
      next
    }

    rhs <- paste(c(group_col, model_covariates), collapse = " + ")
    fml <- stats::as.formula(
      paste0("`", p, "` ~ ", rhs, " + (1 | ", id_col, ")")
    )

    fit <- tryCatch({
      if (use_lmerTest) {
        lmerTest::lmer(fml, data = dat, REML = FALSE)
      } else {
        lme4::lmer(fml, data = dat, REML = FALSE)
      }
    }, error = function(e) NULL)

    if (is.null(fit)) {
      out[[i]] <- dplyr::tibble(Protein = p, N = nrow(dat), coef = NA_real_,
        se = NA_real_, t = NA_real_, p_value = NA_real_, Status = "Model error")
      next
    }

    cf <- as.data.frame(summary(fit)$coefficients)
    rn <- rownames(cf)
    target_term <- paste0(group_col, case_level)
    if (!target_term %in% rn) {
      idx <- grep(paste0("^", group_col), rn)
      if (length(idx) == 0) {
        out[[i]] <- dplyr::tibble(Protein = p, N = nrow(dat), coef = NA_real_,
          se = NA_real_, t = NA_real_, p_value = NA_real_, Status = "No group term")
        next
      }
      target_term <- rn[idx[1]]
    }

    p_col <- if ("Pr(>|t|)" %in% colnames(cf)) "Pr(>|t|)" else "Pr(>|z|)"
    out[[i]] <- dplyr::tibble(
      Protein = p,
      N = nrow(dat),
      coef = cf[target_term, "Estimate"],
      se = cf[target_term, "Std. Error"],
      t = cf[target_term, "t value"],
      p_value = if (p_col %in% colnames(cf)) cf[target_term, p_col] else NA_real_,
      Status = "Success"
    )
  }

  dplyr::bind_rows(out) |>
    dplyr::mutate(
      fdr = ifelse(is.na(.data$p_value), NA_real_,
                   stats::p.adjust(.data$p_value, method = "BH")),
      direction = dplyr::case_when(
        is.na(.data$coef) ~ NA_character_,
        .data$coef > 0 ~ "UP",
        .data$coef < 0 ~ "DOWN",
        TRUE ~ "NO"
      ),
      significant = dplyr::if_else(
        !is.na(.data$fdr) & .data$fdr < 0.05, "Significant", "Not Significant"
      )
    ) |>
    dplyr::arrange(.data$p_value)
}
