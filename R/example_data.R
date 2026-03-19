#' Load bundled example dataset
#'
#' @return A data frame for trajectory/survival/prediction demos.
load_example_data <- function() {
  f <- system.file("extdata", "trajSurvR_example.csv", package = "trajSurvR")
  if (f == "") stop("Bundled example dataset not found.")
  utils::read.csv(f, stringsAsFactors = FALSE)
}


#' Rename feature columns to generic names
#'
#' Converts legacy marker names from original notebooks into generic names
#' for broader reuse in package demos and templates.
#'
#' @param data Input data frame.
#' @param mapping Named character vector where names are old columns and
#'   values are new generic names.
#'
#' @return Data frame with renamed columns (only columns present are renamed).
standardize_feature_names <- function(
  data,
  mapping = c(
    "NPX_NEFL" = "protein_001",
    "NPX_GPNMB" = "protein_002",
    "NPX_CXCL13" = "protein_003",
    "delta_NEFL" = "delta_001",
    "delta_GPNMB" = "delta_002",
    "delta_CXCL13" = "delta_003"
  )
) {
  stopifnot(is.data.frame(data))
  cols <- intersect(names(mapping), colnames(data))
  if (length(cols) == 0) return(data)
  stats::setNames(data, ifelse(colnames(data) %in% cols, mapping[colnames(data)], colnames(data)))
}


#' Generate larger synthetic example dataset
#'
#' @param n_subjects Number of synthetic subjects.
#' @param n_visits Number of visits per subject.
#' @param seed Random seed.
#' @param generic_feature_names If TRUE, return generic feature names
#'   (`protein_001`, `delta_001`, ...). If FALSE, return legacy names.
#'
#' @return Synthetic longitudinal dataset suitable for package demos.
generate_example_data <- function(
  n_subjects = 180,
  n_visits = 4,
  seed = 2026,
  generic_feature_names = TRUE
) {
  set.seed(seed)
  if (n_subjects < 12) stop("n_subjects should be >= 12.")
  if (n_visits < 2) stop("n_visits should be >= 2.")

  groups <- rep(c("Control", "Potential Convert", "Convert"), length.out = n_subjects)
  uid <- sprintf("U%04d", seq_len(n_subjects))
  sex <- sample(c("Male", "Female"), n_subjects, replace = TRUE)
  geno <- sample(c("Other Genotype", "SOD1 A4V", "SOD1 nonA4V", "C9orf72"), n_subjects, replace = TRUE)
  base_age <- round(rnorm(n_subjects, mean = 52, sd = 8), 1)

  out <- vector("list", n_subjects * n_visits)
  idx <- 1L
  for (i in seq_len(n_subjects)) {
    grp <- groups[i]
    onset_year <- if (grp == "Convert") runif(1, 1.5, 6) else NA_real_
    censor_year <- if (grp != "Convert") runif(1, 3, 8) else NA_real_
    stop_time <- if (grp == "Convert") onset_year else censor_year
    event <- if (grp == "Convert") 1L else 0L

    visit_times <- sort(runif(n_visits, min = 0, max = max(1.2, stop_time - 0.2)))
    for (v in seq_len(n_visits)) {
      t <- visit_times[v]
      coll_age <- base_age[i] + t
      weight <- round(58 + 0.25 * coll_age + ifelse(sex[i] == "Male", 6, 0) + rnorm(1, 0, 2), 2)
      height <- round(162 + ifelse(sex[i] == "Male", 10, 0) + rnorm(1, 0, 2), 2)

      progression <- if (grp == "Convert") t / max(onset_year, 1e-6) else t / max(censor_year, 1e-6)
      progression <- pmin(pmax(progression, 0), 1.5)
      signal <- if (grp == "Convert") 0.45 * progression else 0.08 * progression

      delta1 <- round(rnorm(1, mean = signal, sd = 0.07), 3)
      delta2 <- round(rnorm(1, mean = signal * 0.75, sd = 0.07), 3)
      delta3 <- round(rnorm(1, mean = signal * 0.6, sd = 0.07), 3)
      p1 <- round(6.2 + delta1 * 2.4 + rnorm(1, 0, 0.12), 3)
      p2 <- round(5.6 + delta2 * 2.2 + rnorm(1, 0, 0.12), 3)
      p3 <- round(5.0 + delta3 * 2.1 + rnorm(1, 0, 0.12), 3)

      yrs_since_os <- if (grp == "Convert") -(onset_year - t) else NA_real_
      yrs_since_cen <- if (grp != "Convert") -(censor_year - t) else NA_real_
      y_time <- if (grp == "Convert") yrs_since_os else yrs_since_cen
      case <- if (grp == "Convert" && !is.na(yrs_since_os) && (-yrs_since_os) <= 2) 1L else 0L

      out[[idx]] <- data.frame(
        SampleID = sprintf("S_%s_%02d", uid[i], v),
        UIDx = uid[i],
        CollNum = v,
        Group = grp,
        Sex = sex[i],
        GenoGroup = geno[i],
        CollAge = coll_age,
        Weight = weight,
        Height = height,
        YrSinceOs = yrs_since_os,
        YrSinceCen = yrs_since_cen,
        time_since_baseline = t,
        stop_time = stop_time,
        event = event,
        delta_NEFL = delta1,
        delta_GPNMB = delta2,
        delta_CXCL13 = delta3,
        NPX_NEFL = p1,
        NPX_GPNMB = p2,
        NPX_CXCL13 = p3,
        Y = y_time,
        CASE = case,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  dat <- dplyr::bind_rows(out)
  if (generic_feature_names) {
    dat <- standardize_feature_names(dat)
  }
  dat
}


#' Simulate longitudinal trajectory-survival dataset
#'
#' High-level simulation API for demos and testing. This is the preferred
#' function for generating synthetic input data for \code{trajSurvR}.
#'
#' @param n_subjects Number of synthetic subjects.
#' @param n_visits Number of visits per subject.
#' @param seed Random seed.
#' @param generic_feature_names If TRUE, output uses generic feature names
#'   (`protein_001`, `delta_001`, ...).
#'
#' @return Simulated visit-level data frame with outcomes, covariates, and features.
simulate_trajsurv_data <- function(
  n_subjects = 180,
  n_visits = 4,
  seed = 2026,
  generic_feature_names = TRUE
) {
  generate_example_data(
    n_subjects = n_subjects,
    n_visits = n_visits,
    seed = seed,
    generic_feature_names = generic_feature_names
  )
}
