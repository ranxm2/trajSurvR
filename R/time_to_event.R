#' Build baseline survival dataset from visit and protein data
#'
#' @param visit_info Visit data.
#' @param protein_wide Protein matrix keyed by sample ID.
#' @param sample_id_col Sample ID column.
#' @param id_col Subject ID column.
#' @param group_col Group column.
#' @param converter_label Converter group label.
#' @param pre_converter_label Pre-converter/censored group label.
#' @param onset_col Years since onset column.
#' @param censor_col Years since censor column.
#' @param constant_covariates Subject-level covariates retained from baseline visit.
#' @param time_varying_covariates Visit-level covariates retained from baseline visit.
#' @param visit_order_col Visit ordering column used to define baseline visit.
#'
#' @return Data frame with \code{event} and \code{time_yr} plus joined proteins.
build_baseline_survival_data <- function(
  visit_info,
  protein_wide,
  sample_id_col = "SampleID",
  id_col = "UIDx",
  group_col = "Group",
  converter_label = "Convert",
  pre_converter_label = "Potential Convert",
  onset_col = "YrSinceOs",
  censor_col = "YrSinceCen",
  constant_covariates = c("Sex", "GenoGroup"),
  time_varying_covariates = c("CollAge"),
  visit_order_col = "CollNum"
) {
  covariates <- unique(c(constant_covariates, time_varying_covariates))
  req <- c(sample_id_col, id_col, group_col, onset_col, censor_col, visit_order_col, covariates)
  miss <- setdiff(req, colnames(visit_info))
  if (length(miss) > 0) stop("visit_info missing columns: ", paste(miss, collapse = ", "))
  if (!sample_id_col %in% colnames(protein_wide)) {
    stop("protein_wide missing sample column: ", sample_id_col)
  }

  baseline <- visit_info |>
    dplyr::filter(.data[[group_col]] %in% c(converter_label, pre_converter_label)) |>
    dplyr::arrange(.data[[id_col]], .data[[visit_order_col]]) |>
    dplyr::group_by(.data[[id_col]]) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      event = ifelse(.data[[group_col]] == converter_label, 1L, 0L),
      time_yr = dplyr::case_when(
        .data[[group_col]] == converter_label ~ -as.numeric(.data[[onset_col]]),
        .data[[group_col]] == pre_converter_label ~ -as.numeric(.data[[censor_col]]),
        TRUE ~ NA_real_
      )
    ) |>
    dplyr::filter(is.finite(time_yr), !is.na(time_yr), time_yr > 0)

  join_cols <- setdiff(colnames(protein_wide), colnames(baseline))
  protein_join <- protein_wide[, unique(c(sample_id_col, join_cols)), drop = FALSE]
  baseline |>
    dplyr::left_join(protein_join, by = sample_id_col)
}


#' Univariate baseline Cox models per protein
#'
#' @param dat Survival dataset with event/time and proteins.
#' @param proteins Protein names.
#' @param constant_covariates Subject-level covariates to include.
#' @param time_varying_covariates Visit-level covariates to include.
#' @param time_col Survival time column.
#' @param event_col Event indicator column.
#' @param min_rows Minimum complete rows.
#' @param min_events Minimum events.
#'
#' @return Data frame of per-protein Cox statistics.
fit_baseline_cox <- function(
  dat,
  proteins,
  constant_covariates = c("Sex", "GenoGroup"),
  time_varying_covariates = c("CollAge"),
  time_col = "time_yr",
  event_col = "event",
  min_rows = 20,
  min_events = 5
) {
  covariates <- unique(c(constant_covariates, time_varying_covariates))
  req <- c(time_col, event_col, covariates)
  miss <- setdiff(req, colnames(dat))
  if (length(miss) > 0) stop("Missing columns in dat: ", paste(miss, collapse = ", "))

  out <- vector("list", length(proteins))
  for (i in seq_along(proteins)) {
    p <- proteins[i]
    if (!p %in% colnames(dat)) {
      out[[i]] <- dplyr::tibble(Protein = p, N = NA_integer_, N_events = NA_integer_, coef = NA_real_, HR = NA_real_, p_value = NA_real_, Status = "Missing protein")
      next
    }

    tmp <- dat |>
      dplyr::select(dplyr::all_of(c(time_col, event_col, covariates, p))) |>
      tidyr::drop_na()
    n_events <- sum(tmp[[event_col]] == 1)
    if (nrow(tmp) < min_rows || n_events < min_events) {
      out[[i]] <- dplyr::tibble(Protein = p, N = nrow(tmp), N_events = n_events, coef = NA_real_, HR = NA_real_, p_value = NA_real_, Status = "Insufficient data")
      next
    }

    fml <- stats::as.formula(
      paste0("survival::Surv(", time_col, ", ", event_col, ") ~ `", p, "` + ", paste(covariates, collapse = " + "))
    )
    fit <- tryCatch(survival::coxph(fml, data = tmp), error = function(e) NULL)
    if (is.null(fit)) {
      out[[i]] <- dplyr::tibble(Protein = p, N = nrow(tmp), N_events = n_events, coef = NA_real_, HR = NA_real_, p_value = NA_real_, Status = "Model error")
      next
    }

    coef_tab <- summary(fit)$coefficients
    out[[i]] <- dplyr::tibble(
      Protein = p,
      N = nrow(tmp),
      N_events = n_events,
      coef = coef_tab[1, "coef"],
      HR = coef_tab[1, "exp(coef)"],
      p_value = coef_tab[1, "Pr(>|z|)"],
      Status = "Success"
    )
  }

  dplyr::bind_rows(out) |>
    dplyr::mutate(
      fdr = ifelse(is.na(p_value), NA_real_, stats::p.adjust(p_value, method = "BH")),
      direction = dplyr::case_when(coef > 0 ~ "Risk", coef < 0 ~ "Protective", TRUE ~ NA_character_)
    ) |>
    dplyr::arrange(p_value)
}


#' LASSO Cox model for multivariable time-to-event analysis
#'
#' @param dat Data frame containing event/time and predictors.
#' @param proteins Protein columns to penalize.
#' @param constant_covariates Subject-level covariates.
#' @param time_varying_covariates Visit-level covariates.
#' @param time_col Time column.
#' @param event_col Event column.
#' @param nfolds CV folds for \code{cv.glmnet}.
#'
#' @return List with cv fit and selected features at lambda.min/lambda.1se.
fit_lasso_cox <- function(
  dat,
  proteins,
  constant_covariates = c("Sex", "GenoGroup"),
  time_varying_covariates = c("CollAge"),
  time_col = "time_yr",
  event_col = "event",
  nfolds = 5
) {
  covariates <- unique(c(constant_covariates, time_varying_covariates))
  cols <- unique(c(time_col, event_col, covariates, proteins))
  tmp <- dat |>
    dplyr::select(dplyr::all_of(cols)) |>
    tidyr::drop_na()
  if (nrow(tmp) < 30) stop("Too few complete rows for LASSO Cox.")

  x <- stats::model.matrix(
    stats::as.formula(paste0("~ ", paste(c(covariates, proteins), collapse = " + "))),
    data = tmp
  )
  x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
  y <- survival::Surv(tmp[[time_col]], tmp[[event_col]])

  cvfit <- glmnet::cv.glmnet(x = x, y = y, family = "cox", alpha = 1, nfolds = nfolds)
  coef_min <- as.matrix(glmnet::coef.glmnet(cvfit, s = "lambda.min"))
  coef_1se <- as.matrix(glmnet::coef.glmnet(cvfit, s = "lambda.1se"))

  list(
    cv_fit = cvfit,
    lambda_min = cvfit$lambda.min,
    lambda_1se = cvfit$lambda.1se,
    selected_min = rownames(coef_min)[coef_min[, 1] != 0],
    selected_1se = rownames(coef_1se)[coef_1se[, 1] != 0]
  )
}


#' Construct counting-process intervals from visit table
#'
#' @param visit_df Visit data including visit time and stop time.
#' @param id_col Subject ID.
#' @param visit_time_col Time since baseline at each visit.
#' @param stop_time_col Subject-level stop time.
#' @param event_col Subject-level event.
#' @param keep_visit_cols Visit-level columns to carry onto each interval start.
#'
#' @return Interval-level data with \code{tstart}, \code{tstop}, and \code{event}.
make_counting_process_data <- function(
  visit_df,
  id_col = "UIDx",
  visit_time_col = "time_since_baseline",
  stop_time_col = "stop_time",
  event_col = "event",
  keep_visit_cols = character(0)
) {
  req <- c(id_col, visit_time_col, stop_time_col, event_col, keep_visit_cols)
  miss <- setdiff(req, colnames(visit_df))
  if (length(miss) > 0) stop("visit_df missing: ", paste(miss, collapse = ", "))

  split_subj <- split(visit_df, visit_df[[id_col]])
  intervals <- lapply(split_subj, function(subj) {
    subj <- subj[order(subj[[visit_time_col]]), , drop = FALSE]
    stop_time <- unique(subj[[stop_time_col]])[1]
    has_event <- unique(subj[[event_col]])[1]
    times <- unique(subj[[visit_time_col]])
    times <- times[times < stop_time]
    if (!is.finite(stop_time) || length(times) == 0) return(NULL)

    tstart <- times
    tstop <- c(times[-1], stop_time)
    out <- data.frame(
      UIDx = subj[[id_col]][1],
      tstart = tstart,
      tstop = tstop,
      event = 0L
    )
    if (isTRUE(has_event == 1)) out$event[nrow(out)] <- 1L
    if (length(keep_visit_cols) > 0) {
      keep_map <- subj[, c(visit_time_col, keep_visit_cols), drop = FALSE]
      colnames(keep_map)[1] <- "tstart"
      out <- merge(out, keep_map, by = "tstart", all.x = TRUE, sort = FALSE)
      out <- out[order(out$tstart), , drop = FALSE]
      out[[id_col]] <- subj[[id_col]][1]
    }
    out <- out[out$tstart < out$tstop, , drop = FALSE]
    out
  })

  dplyr::bind_rows(intervals)
}


#' Time-dependent Cox model per protein (counting-process data)
#'
#' @param td_data Data with \code{tstart}, \code{tstop}, \code{event}, proteins, covariates.
#' @param proteins Protein names.
#' @param constant_covariates Subject-level covariates.
#' @param time_varying_covariates Visit-level covariates.
#' @param id_col Subject ID column.
#' @param tstart_col Interval start column.
#' @param tstop_col Interval stop column.
#' @param event_col Event indicator column.
#' @param min_subjects Minimum subjects to fit.
#' @param min_events Minimum events to fit.
#'
#' @return Data frame with per-protein time-dependent Cox statistics.
fit_time_dependent_cox <- function(
  td_data,
  proteins,
  constant_covariates = c("Sex", "GenoGroup"),
  time_varying_covariates = c("CollAge"),
  id_col = "UIDx",
  tstart_col = "tstart",
  tstop_col = "tstop",
  event_col = "event",
  min_subjects = 10,
  min_events = 5
) {
  covariates <- unique(c(constant_covariates, time_varying_covariates))
  req <- c(id_col, tstart_col, tstop_col, event_col, covariates)
  miss <- setdiff(req, colnames(td_data))
  if (length(miss) > 0) stop("td_data missing columns: ", paste(miss, collapse = ", "))

  out <- vector("list", length(proteins))
  for (i in seq_along(proteins)) {
    p <- proteins[i]
    if (!p %in% colnames(td_data)) {
      out[[i]] <- dplyr::tibble(Protein = p, N_subjects = NA_integer_, N_events = NA_integer_, coef = NA_real_, HR = NA_real_, p_value = NA_real_, Status = "Missing protein")
      next
    }
    dat <- td_data |>
      dplyr::select(dplyr::all_of(c(id_col, tstart_col, tstop_col, event_col, covariates, p))) |>
      tidyr::drop_na()
    n_subj <- dplyr::n_distinct(dat[[id_col]])
    n_evt <- dat |>
      dplyr::group_by(.data[[id_col]]) |>
      dplyr::summarise(evt = max(.data[[event_col]]), .groups = "drop") |>
      dplyr::summarise(n = sum(evt), .groups = "drop") |>
      dplyr::pull(n)

    if (n_subj < min_subjects || n_evt < min_events) {
      out[[i]] <- dplyr::tibble(Protein = p, N_subjects = n_subj, N_events = n_evt, coef = NA_real_, HR = NA_real_, p_value = NA_real_, Status = "Insufficient data")
      next
    }

    fml <- stats::as.formula(
      paste0("survival::Surv(", tstart_col, ", ", tstop_col, ", ", event_col, ") ~ `", p, "` + ", paste(covariates, collapse = " + "))
    )
    fit <- tryCatch(survival::coxph(fml, data = dat), error = function(e) NULL)
    if (is.null(fit)) {
      out[[i]] <- dplyr::tibble(Protein = p, N_subjects = n_subj, N_events = n_evt, coef = NA_real_, HR = NA_real_, p_value = NA_real_, Status = "Model error")
      next
    }
    co <- summary(fit)$coefficients
    out[[i]] <- dplyr::tibble(
      Protein = p,
      N_subjects = n_subj,
      N_events = n_evt,
      coef = co[1, "coef"],
      HR = co[1, "exp(coef)"],
      p_value = co[1, "Pr(>|z|)"],
      Status = "Success"
    )
  }

  dplyr::bind_rows(out) |>
    dplyr::mutate(fdr = ifelse(is.na(p_value), NA_real_, stats::p.adjust(p_value, method = "BH"))) |>
    dplyr::arrange(p_value)
}
