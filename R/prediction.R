#' Create binary event label for onset-within-horizon prediction
#'
#' @param visit_df Visit table.
#' @param horizon_years Prediction horizon in years.
#' @param group_col Group column.
#' @param converter_label Converter label.
#' @param control_label Control/non-converter label.
#' @param onset_col Years since onset.
#'
#' @return Input data with \code{CASE} column.
make_event_outcome <- function(
  visit_df,
  horizon_years = 2,
  group_col = "Group",
  converter_label = "Convert",
  control_label = c("Control", "Potential Convert"),
  onset_col = "YrSinceOs"
) {
  dat <- dplyr::as_tibble(visit_df)
  if (!all(c(group_col, onset_col) %in% colnames(dat))) {
    stop("visit_df must contain group and onset columns.")
  }
  dat |>
    dplyr::mutate(
      onset_in_years = -as.numeric(.data[[onset_col]]),
      CASE = dplyr::case_when(
        .data[[group_col]] == converter_label & onset_in_years <= horizon_years ~ 1L,
        .data[[group_col]] %in% control_label ~ 0L,
        TRUE ~ NA_integer_
      )
    )
}


#' Cross-validated AUROC for logistic regression
#'
#' @param data Data frame.
#' @param feature_cols Predictors in logistic model.
#' @param outcome_col Binary outcome column.
#' @param fold_col Optional precomputed fold column.
#' @param n_folds Number of folds if \code{fold_col} missing.
#' @param seed Random seed.
#'
#' @return Scalar AUROC.
cv_logistic_auc <- function(
  data,
  feature_cols,
  outcome_col = "CASE",
  fold_col = NULL,
  n_folds = 5,
  seed = 2026
) {
  req <- c(feature_cols, outcome_col)
  miss <- setdiff(req, colnames(data))
  if (length(miss) > 0) stop("Missing data columns: ", paste(miss, collapse = ", "))

  dat <- data |>
    dplyr::select(dplyr::all_of(req), dplyr::any_of(fold_col)) |>
    tidyr::drop_na()
  if (nrow(dat) < 30) return(NA_real_)

  if (is.null(fold_col) || !fold_col %in% colnames(dat)) {
    set.seed(seed)
    dat$Fold <- sample(rep(seq_len(n_folds), length.out = nrow(dat)))
    fold_col <- "Fold"
  }

  preds <- numeric(nrow(dat))
  truth <- dat[[outcome_col]]
  for (k in sort(unique(dat[[fold_col]]))) {
    idx_test <- which(dat[[fold_col]] == k)
    idx_train <- which(dat[[fold_col]] != k)
    train_df <- dat[idx_train, c(feature_cols, outcome_col), drop = FALSE]
    test_df <- dat[idx_test, feature_cols, drop = FALSE]

    fml <- stats::as.formula(paste0(outcome_col, " ~ ", paste(feature_cols, collapse = " + ")))
    fit <- tryCatch(stats::glm(fml, data = train_df, family = "binomial"), error = function(e) NULL)
    if (is.null(fit)) {
      preds[idx_test] <- NA_real_
      next
    }
    preds[idx_test] <- as.numeric(stats::predict(fit, newdata = test_df, type = "response"))
  }

  ok <- is.finite(preds) & !is.na(truth)
  if (sum(ok) < 10 || length(unique(truth[ok])) < 2) return(NA_real_)
  roc_obj <- pROC::roc(response = truth[ok], predictor = preds[ok], quiet = TRUE)
  as.numeric(pROC::auc(roc_obj))
}


#' Forward greedy feature selection
#'
#' @param data Modeling dataset.
#' @param candidate_features Candidate variables.
#' @param base_features Always-included variables.
#' @param evaluator Function accepting \code{(data, features)} and returning score.
#' @param min_improvement Minimal score gain to continue.
#' @param max_steps Maximum selection steps.
#'
#' @return A list with selected features and step history.
forward_greedy_selection <- function(
  data,
  candidate_features,
  base_features = character(0),
  evaluator = function(dat, feats) cv_logistic_auc(dat, feats),
  min_improvement = 1e-4,
  max_steps = 20
) {
  selected <- character(0)
  history <- dplyr::tibble(Step = integer(), Added = character(), Score = numeric(), Delta = numeric())
  current_score <- evaluator(data, c(base_features, selected))
  if (!is.finite(current_score)) current_score <- -Inf

  for (step in seq_len(min(max_steps, length(candidate_features)))) {
    remaining <- setdiff(candidate_features, selected)
    if (length(remaining) == 0) break

    scores <- vapply(remaining, function(cand) {
      evaluator(data, c(base_features, selected, cand))
    }, numeric(1))

    finite_idx <- which(is.finite(scores))
    if (length(finite_idx) == 0) break
    best_idx <- finite_idx[which.max(scores[finite_idx])]
    best_cand <- remaining[best_idx]
    best_score <- scores[best_idx]
    delta <- best_score - current_score
    if (!is.finite(delta) || delta < min_improvement) break

    selected <- c(selected, best_cand)
    history <- dplyr::bind_rows(history, dplyr::tibble(
      Step = step, Added = best_cand, Score = best_score, Delta = delta
    ))
    current_score <- best_score
  }

  list(selected = selected, history = history, final_score = current_score)
}


#' LASSO feature selection helper (binomial/gaussian/cox)
#'
#' @param x Predictor matrix or data frame.
#' @param y Response vector (or \code{Surv} object for cox).
#' @param family glmnet family.
#' @param nfolds Number of CV folds.
#' @param alpha Elastic-net mixing (1 = LASSO).
#'
#' @return List with cv fit and selected features.
select_features_lasso <- function(x, y, family = c("binomial", "gaussian", "cox"), nfolds = 5, alpha = 1) {
  family <- match.arg(family)
  x <- as.matrix(x)
  fit <- glmnet::cv.glmnet(x = x, y = y, family = family, alpha = alpha, nfolds = nfolds)
  coef_min <- as.matrix(glmnet::coef.glmnet(fit, s = "lambda.min"))
  coef_1se <- as.matrix(glmnet::coef.glmnet(fit, s = "lambda.1se"))
  list(
    cv_fit = fit,
    selected_min = rownames(coef_min)[coef_min[, 1] != 0],
    selected_1se = rownames(coef_1se)[coef_1se[, 1] != 0]
  )
}


#' Fit logistic LASSO for event prediction
#'
#' @param data Data frame.
#' @param outcome_col Binary outcome column.
#' @param feature_cols Predictor columns.
#' @param nfolds CV folds.
#'
#' @return LASSO feature-selection result for binomial model.
fit_event_prediction_lasso <- function(data, outcome_col = "CASE", feature_cols, nfolds = 5) {
  req <- c(outcome_col, feature_cols)
  miss <- setdiff(req, colnames(data))
  if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))
  tmp <- data |>
    dplyr::select(dplyr::all_of(req)) |>
    tidyr::drop_na()
  x <- as.matrix(tmp[, feature_cols, drop = FALSE])
  y <- tmp[[outcome_col]]
  select_features_lasso(x, y, family = "binomial", nfolds = nfolds, alpha = 1)
}


#' Onset-time regression with baseline covariates
#'
#' @param data Data frame.
#' @param outcome_col Continuous onset-time column.
#' @param feature_cols Protein features.
#' @param constant_covariates Subject-level covariates to always include.
#' @param time_varying_covariates Visit-level covariates to always include.
#' @param id_col Subject ID for LOOCV.
#'
#' @return List with fitted model, LOOCV predictions, and metrics.
fit_onset_time_regression <- function(
  data,
  outcome_col = "Y",
  feature_cols,
  constant_covariates = c("Sex", "GenoGroup"),
  time_varying_covariates = c("CollAge"),
  id_col = "UIDx"
) {
  baseline_cols <- unique(c(constant_covariates, time_varying_covariates))
  req <- unique(c(outcome_col, id_col, feature_cols, baseline_cols))
  miss <- setdiff(req, colnames(data))
  if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

  dat <- data |>
    dplyr::select(dplyr::all_of(req)) |>
    tidyr::drop_na()
  fml <- stats::as.formula(
    paste0(outcome_col, " ~ ", paste(c(feature_cols, baseline_cols), collapse = " + "))
  )
  fit <- stats::lm(fml, data = dat)

  ids <- unique(dat[[id_col]])
  pred_df <- vector("list", length(ids))
  for (i in seq_along(ids)) {
    uid <- ids[i]
    tr <- dat[dat[[id_col]] != uid, , drop = FALSE]
    te <- dat[dat[[id_col]] == uid, , drop = FALSE]
    fit_i <- tryCatch(stats::lm(fml, data = tr), error = function(e) NULL)
    pred <- if (is.null(fit_i)) rep(NA_real_, nrow(te)) else as.numeric(stats::predict(fit_i, newdata = te))
    pred_df[[i]] <- dplyr::tibble(UIDx = uid, y_true = te[[outcome_col]], y_pred = pred)
  }
  cv <- dplyr::bind_rows(pred_df)
  ok <- is.finite(cv$y_true) & is.finite(cv$y_pred)
  resid <- cv$y_true[ok] - cv$y_pred[ok]
  metrics <- dplyr::tibble(
    MAE = mean(abs(resid)),
    RMSE = sqrt(mean(resid^2)),
    Correlation = suppressWarnings(stats::cor(cv$y_true[ok], cv$y_pred[ok])),
    Within1yr = mean(abs(resid) < 1) * 100
  )

  list(model = fit, loocv_predictions = cv, metrics = metrics)
}
