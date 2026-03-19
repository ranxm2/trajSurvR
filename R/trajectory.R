#' Fit trajectory model for one protein delta
#'
#' Implements the GAM/GAMM trajectory approach used in the raw Rmd:
#' \code{delta ~ s(YrSinceOs)} with random intercept by subject where possible.
#'
#' @param delta_df Data frame with time, delta, and subject ID.
#' @param time_col Time-to-onset column.
#' @param delta_col Delta response column.
#' @param id_col Subject ID column.
#' @param covariates Additional covariates to adjust in trajectory model.
#' @param k Basis dimension for smooth term.
#' @param min_obs Minimum observations required.
#' @param n_grid Number of grid points for smooth prediction.
#'
#' @return List with model, prediction grid, p-value, and significance intervals.
fit_trajectory_gam <- function(
  delta_df,
  time_col = "YrSinceOs",
  delta_col = "delta",
  id_col = "UIDx",
  covariates = character(0),
  k = 10,
  min_obs = 10,
  n_grid = 300
) {
  req <- c(time_col, delta_col, id_col, covariates)
  miss <- setdiff(req, colnames(delta_df))
  if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))

  dat <- dplyr::as_tibble(delta_df) |>
    dplyr::select(dplyr::all_of(req)) |>
    dplyr::filter(!is.na(.data[[time_col]]), !is.na(.data[[delta_col]]))
  if (nrow(dat) < min_obs) return(NULL)

  dat[[id_col]] <- as.factor(dat[[id_col]])
  rhs_terms <- c(paste0("s(", time_col, ", k = ", k, ")"), covariates)
  formula_str <- paste0(delta_col, " ~ ", paste(rhs_terms, collapse = " + "))
  fml <- stats::as.formula(formula_str)

  mod <- tryCatch(
    mgcv::gamm(fml, random = stats::setNames(list(~1), id_col), data = dat, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(mod)) {
    mod <- tryCatch(list(gam = mgcv::gam(fml, data = dat, method = "REML")), error = function(e) NULL)
  }
  if (is.null(mod)) return(NULL)

  gam_obj <- if (inherits(mod, "gam")) mod else mod$gam
  st <- tryCatch(summary(gam_obj)$s.table, error = function(e) NULL)
  pval <- NA_real_
  if (!is.null(st)) {
    smooth_name <- paste0("s(", time_col, ")")
    if (smooth_name %in% rownames(st)) {
      pval <- st[smooth_name, "p-value"]
    } else if (nrow(st) > 0) {
      pval <- st[1, "p-value"]
    }
  }

  grid <- data.frame(
    YrSinceOs = seq(min(dat[[time_col]], na.rm = TRUE), max(dat[[time_col]], na.rm = TRUE), length.out = n_grid)
  )
  colnames(grid)[1] <- time_col
  if (length(covariates) > 0) {
    for (cv in covariates) {
      if (is.numeric(dat[[cv]])) {
        grid[[cv]] <- stats::median(dat[[cv]], na.rm = TRUE)
      } else {
        lev <- names(sort(table(dat[[cv]]), decreasing = TRUE))
        fill_val <- if (length(lev) > 0) lev[1] else as.character(dat[[cv]][1])
        if (is.factor(dat[[cv]])) {
          grid[[cv]] <- factor(fill_val, levels = levels(dat[[cv]]))
        } else {
          grid[[cv]] <- fill_val
        }
      }
    }
  }
  pred <- stats::predict(gam_obj, newdata = grid, se.fit = TRUE)
  grid$fit <- as.numeric(pred$fit)
  grid$se <- as.numeric(pred$se.fit)
  grid$lower <- grid$fit - 1.96 * grid$se
  grid$upper <- grid$fit + 1.96 * grid$se
  grid$sig <- dplyr::case_when(
    grid$lower > 0 ~ 1L,
    grid$upper < 0 ~ -1L,
    TRUE ~ 0L
  )

  r <- rle(grid$sig)
  ends <- cumsum(r$lengths)
  starts <- c(1, utils::head(ends, -1) + 1)
  nzi <- which(r$values != 0)
  intervals <- dplyr::tibble(
    start = grid[[time_col]][starts[nzi]],
    end = grid[[time_col]][ends[nzi]],
    direction = ifelse(r$values[nzi] > 0, "increase", "decrease")
  )

  list(
    model = mod,
    gam = gam_obj,
    data = dat,
    predictions = dplyr::as_tibble(grid),
    p_value = pval,
    intervals = intervals
  )
}


#' Run trajectory analysis for multiple proteins
#'
#' @param delta_wide Wide delta matrix with sample ID + protein columns.
#' @param visit_design Visit-level design data with sample ID, time, and subject ID.
#' @param sample_id_col Sample ID column name.
#' @param proteins Optional vector of proteins. Default: all non-sample columns.
#' @param id_col Subject ID column in visit data.
#' @param time_col Time column in visit data.
#' @param constant_covariates Subject-level covariates from visit data.
#' @param time_varying_covariates Visit-level covariates from visit data.
#' @param min_obs Minimum observations to fit one protein.
#'
#' @return List with per-protein models, summary table, and combined predictions.
run_trajectory_batch <- function(
  delta_wide,
  visit_design,
  sample_id_col = "SampleID",
  proteins = NULL,
  id_col = "UIDx",
  time_col = "YrSinceOs",
  constant_covariates = character(0),
  time_varying_covariates = character(0),
  min_obs = 10
) {
  covariates <- unique(c(constant_covariates, time_varying_covariates))
  req_visit <- c(sample_id_col, id_col, time_col, covariates)
  miss_visit <- setdiff(req_visit, colnames(visit_design))
  if (length(miss_visit) > 0) {
    stop("visit_design missing required columns: ", paste(miss_visit, collapse = ", "))
  }
  if (!sample_id_col %in% colnames(delta_wide)) {
    stop("delta_wide missing sample ID column: ", sample_id_col)
  }

  if (is.null(proteins)) proteins <- setdiff(colnames(delta_wide), sample_id_col)
  proteins <- intersect(proteins, colnames(delta_wide))
  if (length(proteins) == 0) stop("No proteins available in delta_wide.")

  fits <- list()
  pred_all <- list()
  sum_rows <- list()

  for (p in proteins) {
    dat <- delta_wide |>
      dplyr::select(dplyr::all_of(sample_id_col), dplyr::all_of(p)) |>
      dplyr::rename(delta = dplyr::all_of(p)) |>
      dplyr::left_join(visit_design[, req_visit, drop = FALSE], by = sample_id_col)

    fit <- fit_trajectory_gam(
      dat,
      time_col = time_col,
      delta_col = "delta",
      id_col = id_col,
      covariates = covariates,
      min_obs = min_obs
    )
    if (is.null(fit)) next
    fits[[p]] <- fit
    pred_all[[p]] <- dplyr::mutate(fit$predictions, Protein = p)
    sum_rows[[p]] <- dplyr::tibble(
      Protein = p,
      N_obs = nrow(fit$data),
      N_subjects = dplyr::n_distinct(fit$data[[id_col]]),
      p_value = fit$p_value,
      fdr = NA_real_,
      first_sig_time = ifelse(nrow(fit$intervals) > 0, min(fit$intervals$start), NA_real_)
    )
  }

  summary_df <- dplyr::bind_rows(sum_rows)
  if (nrow(summary_df) > 0) {
    summary_df$fdr <- stats::p.adjust(summary_df$p_value, method = "BH")
    summary_df <- dplyr::arrange(summary_df, p_value)
  }

  list(
    models = fits,
    summary = summary_df,
    predictions = dplyr::bind_rows(pred_all)
  )
}
