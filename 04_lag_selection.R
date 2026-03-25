###############################################################################
#  04_lag_selection.R  –  Lag-Order Selection for Individual VARX* Models
#
#  For each unit i the individual model is a VARX*(p_i, q_i):
#
#    x_{it} = c_i + A_{i1} x_{i,t-1} + ... + A_{ip} x_{i,t-p}
#                 + B_{i0} x*_{it}    + ... + B_{iq} x*_{i,t-q}  + u_{it}
#
#  We select (p_i, q_i) by evaluating AIC and BIC over a grid of candidate
#  orders, holding q constant while searching p, and vice versa, then
#  doing a joint grid search.
#
#  Contents
#  --------
#  1. select_lag_single()   – AIC/BIC grid search for one unit
#  2. select_lags_all()     – loop over all units and return a summary table
###############################################################################


# ───────────────────────────────────────────────────────────────────────────────
# 1.  Lag Selection for a Single Unit
# ───────────────────────────────────────────────────────────────────────────────

#' Grid-search over (p, q) lag orders for a single unit's VARX* model.
#'
#' For every candidate pair (p, q) ∈ {1,...,max_p} × {1,...,max_q} we
#' estimate the individual model by OLS and record AIC and BIC.
#'
#' @param domestic  T × k_i  matrix of domestic variables
#' @param star      T × k*   matrix of star variables
#' @param max_p     Integer; maximum domestic lag order to test (default 4)
#' @param max_q     Integer; maximum foreign  lag order to test (default 4)
#' @param ic        Character; "aic", "bic", or "both" (default "both")
#' @return  A list with:
#'   \item{grid}{Data frame of all (p, q, aic, bic) evaluations}
#'   \item{best_aic}{Named vector c(p=..., q=...) minimising AIC}
#'   \item{best_bic}{Named vector c(p=..., q=...) minimising BIC}
select_lag_single <- function(domestic, star, max_p = 4, max_q = 4,
                               ic = "both", deterministic = "intercept") {

  domestic <- as.matrix(domestic)
  star     <- as.matrix(star)

  det_info <- build_deterministic_columns(1, deterministic)
  n_det <- det_info$n_det

  # Pre-allocate results grid
  grid <- expand.grid(p = 1:max_p, q = 1:max_q)
  grid$aic <- NA_real_
  grid$bic <- NA_real_

  for (row in seq_len(nrow(grid))) {
    p <- grid$p[row]
    q <- grid$q[row]

    # Use prepare_country_data() from 02_data_preparation.R
    cd <- tryCatch(
      prepare_country_data(domestic, star, p_lag = p, q_lag = q),
      error = function(e) NULL
    )
    if (is.null(cd)) next   # skip infeasible lag combos

    # Build deterministic columns and prepend
    det <- build_deterministic_columns(nrow(cd$Y), deterministic)
    X_full <- if (!is.null(det$D_mat)) cbind(det$D_mat, cd$X) else cd$X

    # Estimate by OLS
    fit <- tryCatch(
      ols_estimate(cd$Y, X_full, intercept = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit)) next

    # Count parameters per equation:
    #   n_det + p * k_dom + (1 + q) * k_star   [contemporaneous star included]
    k_dom  <- cd$k_dom
    k_star <- cd$k_star
    n_params_per_eq <- n_det + p * k_dom + (1 + q) * k_star
    n_params_total  <- n_params_per_eq * k_dom   # system-wide

    ic_val <- information_criteria(fit$residuals, n_params_total)
    grid$aic[row] <- ic_val$aic
    grid$bic[row] <- ic_val$bic
  }

  # Remove rows that failed

  grid <- grid[complete.cases(grid), ]

  if (nrow(grid) == 0) {
    warning("No feasible lag combination found. Returning defaults (1,1).")
    return(list(
      grid     = grid,
      best_aic = c(p = 1, q = 1),
      best_bic = c(p = 1, q = 1)
    ))
  }

  best_aic_row <- grid[which.min(grid$aic), ]
  best_bic_row <- grid[which.min(grid$bic), ]

  return(list(
    grid     = grid,
    best_aic = c(p = best_aic_row$p, q = best_aic_row$q),
    best_bic = c(p = best_bic_row$p, q = best_bic_row$q)
  ))
}


# ───────────────────────────────────────────────────────────────────────────────
# 2.  Lag Selection Across All Units
# ───────────────────────────────────────────────────────────────────────────────

#' Select lag orders for all units and return a tidy summary.
#'
#' @param data_list   Named list of T × k_i matrices (domestic)
#' @param star_list   Named list of T × k*  matrices (star variables)
#' @param max_p       Integer; max domestic lag
#' @param max_q       Integer; max foreign lag
#' @param criterion   "aic" or "bic" – which criterion to use for the
#'                    recommended lag order (default "bic")
#' @return  A list with:
#'   \item{summary_table}{Data frame: unit, best_p, best_q, aic, bic}
#'   \item{p_vec}{Named integer vector of chosen p lags per unit}
#'   \item{q_vec}{Named integer vector of chosen q lags per unit}
#'   \item{details}{Full per-unit grid results}
select_lags_all <- function(data_list, star_list,
                             max_p = 4, max_q = 4,
                             criterion = "bic",
                             deterministic = "intercept") {

  print_banner("Lag-Order Selection (AIC & BIC)")

  unit_names <- names(data_list)
  N <- length(unit_names)

  summary_df <- data.frame(
    unit   = character(N),
    p_aic  = integer(N),
    q_aic  = integer(N),
    p_bic  = integer(N),
    q_bic  = integer(N),
    stringsAsFactors = FALSE
  )

  details <- list()

  for (i in seq_len(N)) {
    u <- unit_names[i]
    message(sprintf("  Searching lags for unit '%s' ...", u))

    res <- select_lag_single(
      domestic      = data_list[[u]],
      star          = star_list[[u]],
      max_p         = max_p,
      max_q         = max_q,
      deterministic = deterministic
    )

    summary_df$unit[i]  <- u
    summary_df$p_aic[i] <- res$best_aic["p"]
    summary_df$q_aic[i] <- res$best_aic["q"]
    summary_df$p_bic[i] <- res$best_bic["p"]
    summary_df$q_bic[i] <- res$best_bic["q"]

    details[[u]] <- res
  }

  # Print summary
  cat("\n")
  print(summary_df, row.names = FALSE)
  cat("\n")

  # Build the recommended lag vectors
  if (criterion == "aic") {
    p_vec <- setNames(summary_df$p_aic, summary_df$unit)
    q_vec <- setNames(summary_df$q_aic, summary_df$unit)
  } else {
    p_vec <- setNames(summary_df$p_bic, summary_df$unit)
    q_vec <- setNames(summary_df$q_bic, summary_df$unit)
  }

  message(sprintf(
    "[GVAR] Recommended lags (criterion = %s): p = {%s}, q = {%s}",
    toupper(criterion),
    paste(p_vec, collapse = ","),
    paste(q_vec, collapse = ",")
  ))

  return(list(
    summary_table = summary_df,
    p_vec         = p_vec,
    q_vec         = q_vec,
    details       = details
  ))
}
