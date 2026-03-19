#' Differential protein analysis via mixed-effects models
#'
#' Mirrors the mixed-model strategy used in the raw differential-expression Rmd:
#' per-protein model with repeated measurements by participant.
#'
#' @param analysis_df Data with subject ID, outcome group, covariates, and proteins.
#' @param proteins Character vector of protein columns to test.
#' @param id_col Subject ID column.
#' @param group_col Group column (e.g. "Group", with Control as reference).
#' @param constant_covariates Subject-level covariates that do not change over visits
#'   (e.g. sex, genotype).
#' @param time_varying_covariates Visit-level covariates that can change by visit
#'   (e.g. age, weight, height).
#' @param case_level Level to compare against reference.
#' @param ref_level Reference level for group factor.
#' @param min_rows Minimum complete rows required for a model.
#'
#' @return Data frame of coefficients, p-values, and FDR.
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
      out[[i]] <- dplyr::tibble(
        Protein = p, N = NA_integer_, coef = NA_real_, se = NA_real_,
        t = NA_real_, p_value = NA_real_, Status = "Protein missing"
      )
      next
    }

    dat <- dat0 |>
      dplyr::select(dplyr::all_of(c(id_col, group_col, model_covariates, p))) |>
      tidyr::drop_na()
    if (nrow(dat) < min_rows || dplyr::n_distinct(dat[[id_col]]) < 8) {
      out[[i]] <- dplyr::tibble(
        Protein = p, N = nrow(dat), coef = NA_real_, se = NA_real_,
        t = NA_real_, p_value = NA_real_, Status = "Insufficient data"
      )
      next
    }

    rhs_terms <- c(group_col, model_covariates)
    fml <- stats::as.formula(paste0("`", p, "` ~ ", paste(rhs_terms, collapse = " + "), " + (1|", id_col, ")"))

    fit <- tryCatch(
      lmerTest::lmer(fml, data = dat, REML = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      out[[i]] <- dplyr::tibble(
        Protein = p, N = nrow(dat), coef = NA_real_, se = NA_real_,
        t = NA_real_, p_value = NA_real_, Status = "Model error"
      )
      next
    }

    cf <- as.data.frame(summary(fit)$coefficients)
    rn <- rownames(cf)
    target_term <- paste0(group_col, case_level)
    if (!target_term %in% rn) {
      # fallback to first non-intercept group term
      target_idx <- grep(paste0("^", group_col), rn)
      if (length(target_idx) == 0) {
        out[[i]] <- dplyr::tibble(
          Protein = p, N = nrow(dat), coef = NA_real_, se = NA_real_,
          t = NA_real_, p_value = NA_real_, Status = "No group term"
        )
        next
      }
      target_term <- rn[target_idx[1]]
    }

    out[[i]] <- dplyr::tibble(
      Protein = p,
      N = nrow(dat),
      coef = cf[target_term, "Estimate"],
      se = cf[target_term, "Std. Error"],
      t = cf[target_term, "t value"],
      p_value = cf[target_term, "Pr(>|t|)"],
      Status = "Success"
    )
  }

  res <- dplyr::bind_rows(out) |>
    dplyr::mutate(
      fdr = ifelse(is.na(p_value), NA_real_, stats::p.adjust(p_value, method = "BH")),
      direction = dplyr::case_when(
        is.na(coef) ~ NA_character_,
        coef > 0 ~ "UP",
        coef < 0 ~ "DOWN",
        TRUE ~ "NO"
      ),
      significant = dplyr::case_when(
        !is.na(fdr) & fdr < 0.05 ~ "Significant",
        TRUE ~ "Not Significant"
      )
    ) |>
    dplyr::arrange(p_value)

  res
}
