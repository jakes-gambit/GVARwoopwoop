###############################################################################
#  07_impulse_response.R  –  Generalised Impulse Response Functions (GIRF)
#
#  In a GVAR context the standard approach is to compute Generalised Impulse
#  Responses (Pesaran & Shin, 1998) which do not require a Cholesky ordering.
#
#  The GIRF to a one-s.e. shock to variable j is:
#
#    GIRF(h, j) = F_comp^h  *  Sigma * e_j  /  sqrt(sigma_jj)
#
#  where e_j is a k-vector with 1 in position j and 0 elsewhere.
#
#  Confidence bands are obtained by bootstrap (resample residuals).
#
#  Contents
#  --------
#  1. girf_analytic()      – compute GIRF from the companion form
#  2. bootstrap_girf()     – bootstrap confidence bands
#  3. plot_irf()           – ggplot2 visualisation of IRFs + bands
#  4. compute_all_irfs()   – convenience wrapper for all shock–response pairs
###############################################################################


# ───────────────────────────────────────────────────────────────────────────────
# 1.  Analytic GIRF
# ───────────────────────────────────────────────────────────────────────────────

#' Compute Generalised IRFs for a one-standard-error shock to variable j.
#'
#' @param gvar_model  Fitted GVAR model (from estimate_gvar)
#' @param shock_var   Integer index (or character name) of the shocked variable
#' @param horizon     Integer; number of periods for the IRF (default 20)
#' @return            Matrix of dimension (horizon+1) × k_total where row h
#'                    gives the impulse response at horizon h (h=0,...,horizon)
girf_analytic <- function(gvar_model, shock_var, horizon = 20) {

  k_total <- gvar_model$k_total
  comp    <- gvar_model$companion
  Sigma   <- gvar_model$Sigma
  kp      <- nrow(comp)

  # Resolve shock_var to an integer index
  if (is.character(shock_var)) {
    shock_idx <- which(gvar_model$var_names == shock_var)
    if (length(shock_idx) == 0)
      stop("Variable '", shock_var, "' not found in the global model.")
  } else {
    shock_idx <- shock_var
  }

  # Standard deviation of the shocked variable
  sigma_jj <- sqrt(Sigma[shock_idx, shock_idx])

  # Selection vector e_j
  e_j <- rep(0, k_total)
  e_j[shock_idx] <- 1

  # Initial impact:  delta_0 = Sigma %*% e_j / sigma_jj
  delta_0 <- as.numeric(Sigma %*% e_j) / sigma_jj

  # Extend delta_0 to companion dimension
  delta_0_comp <- c(delta_0, rep(0, kp - k_total))

  # Iterate forward using the companion matrix
  irf_mat <- matrix(0, nrow = horizon + 1, ncol = k_total)
  irf_mat[1, ] <- delta_0   # h = 0

  state <- delta_0_comp
  for (h in seq_len(horizon)) {
    state <- comp %*% state
    irf_mat[h + 1, ] <- state[1:k_total]
  }

  colnames(irf_mat) <- gvar_model$var_names
  rownames(irf_mat) <- paste0("h=", 0:horizon)

  return(irf_mat)
}


# ───────────────────────────────────────────────────────────────────────────────
# 2.  Bootstrap Confidence Bands
# ───────────────────────────────────────────────────────────────────────────────

#' Bootstrap GIRFs to obtain confidence intervals.
#'
#' Procedure:
#'   1. Resample residuals (with replacement) from the estimated model
#'   2. Simulate new data under the null (keeping coefficients fixed)
#'   3. Re-estimate the GVAR on the simulated data
#'   4. Compute GIRFs from the re-estimated model
#'   5. Repeat B times and extract percentile-based bands
#'
#' For speed, a simplified approach is used:
#'   - We keep the *same* coefficient matrices and only resample the
#'     residual covariance (non-parametric).  This gives bands that
#'     reflect sampling uncertainty in Sigma while treating the point
#'     estimates as fixed.
#'
#' @param gvar_model   Fitted GVAR model
#' @param shock_var    Index or name of shocked variable
#' @param horizon      IRF horizon
#' @param n_boot       Number of bootstrap replications (default 500)
#' @param ci_level     Confidence level (default 0.95 for 95% bands)
#' @param seed         Random seed for reproducibility
#' @return  A list:
#'   \item{median}{(horizon+1) × k  median IRF}
#'   \item{lower}{(horizon+1) × k  lower band}
#'   \item{upper}{(horizon+1) × k  upper band}
#'   \item{point}{(horizon+1) × k  analytic (point-estimate) IRF}
bootstrap_girf <- function(gvar_model, shock_var, horizon = 20,
                            n_boot = 500, ci_level = 0.95, seed = 42) {

  set.seed(seed)

  k_total  <- gvar_model$k_total
  comp     <- gvar_model$companion
  resid    <- gvar_model$residuals   # T_eff × k_total
  TT       <- nrow(resid)

  # Point estimate
  point_irf <- girf_analytic(gvar_model, shock_var, horizon)

  # Storage for bootstrap draws
  boot_array <- array(NA, dim = c(horizon + 1, k_total, n_boot))

  alpha <- 1 - ci_level

  for (b in seq_len(n_boot)) {
    # Resample residuals
    boot_idx <- sample(seq_len(TT), TT, replace = TRUE)
    boot_resid <- resid[boot_idx, ]

    # Recompute covariance from resampled residuals
    Sigma_b <- crossprod(boot_resid) / TT

    # Build a temporary model with this Sigma
    gvar_b <- gvar_model
    gvar_b$Sigma <- Sigma_b

    # Compute GIRF
    boot_array[, , b] <- girf_analytic(gvar_b, shock_var, horizon)
  }

  # Extract percentiles
  lower <- apply(boot_array, c(1, 2), quantile, probs = alpha / 2, na.rm = TRUE)
  upper <- apply(boot_array, c(1, 2), quantile, probs = 1 - alpha / 2, na.rm = TRUE)
  med   <- apply(boot_array, c(1, 2), median, na.rm = TRUE)

  colnames(lower) <- colnames(upper) <- colnames(med) <- gvar_model$var_names
  rownames(lower) <- rownames(upper) <- rownames(med) <- rownames(point_irf)

  return(list(
    point  = point_irf,
    median = med,
    lower  = lower,
    upper  = upper,
    ci_level = ci_level,
    n_boot   = n_boot
  ))
}


# ───────────────────────────────────────────────────────────────────────────────
# 3.  Plotting
# ───────────────────────────────────────────────────────────────────────────────

#' Plot IRF for selected response variables with confidence bands.
#'
#' @param irf_result   Output from bootstrap_girf()
#' @param response_vars  Integer indices or character names of variables to plot
#' @param shock_label   Character; label for the shock (for the title)
#' @return              A ggplot object (also printed)
plot_irf <- function(irf_result, response_vars = NULL, shock_label = "Shock") {

  point <- irf_result$point
  lower <- irf_result$lower
  upper <- irf_result$upper
  ci    <- irf_result$ci_level

  k <- ncol(point)
  H <- nrow(point) - 1
  horizons <- 0:H

  # If no response vars specified, plot all (up to 12)
  if (is.null(response_vars)) {
    response_vars <- seq_len(min(k, 12))
  }
  if (is.character(response_vars)) {
    response_vars <- match(response_vars, colnames(point))
  }

  # Build a long-format data frame for ggplot
  plot_df <- data.frame()
  for (j in response_vars) {
    vname <- colnames(point)[j]
    df_j <- data.frame(
      horizon  = horizons,
      point    = point[, j],
      lower    = lower[, j],
      upper    = upper[, j],
      variable = vname,
      stringsAsFactors = FALSE
    )
    plot_df <- rbind(plot_df, df_j)
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = horizon)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         fill = "steelblue", alpha = 0.25) +
    ggplot2::geom_line(ggplot2::aes(y = point), colour = "steelblue",
                       linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    ggplot2::facet_wrap(~ variable, scales = "free_y") +
    ggplot2::labs(
      title    = paste0("Generalised IRF – ", shock_label),
      subtitle = paste0(ci * 100, "% bootstrap confidence bands (",
                        irf_result$n_boot, " replications)"),
      x = "Horizon", y = "Response"
    ) +
    ggplot2::theme_minimal(base_size = 11)

  print(p)
  return(invisible(p))
}


# ───────────────────────────────────────────────────────────────────────────────
# 4.  Compute IRFs for All Shocks (or a Subset)
# ───────────────────────────────────────────────────────────────────────────────

#' Convenience wrapper: compute bootstrap GIRFs for multiple shock variables.
#'
#' @param gvar_model    Fitted GVAR model
#' @param shock_vars    Integer or character vector of variables to shock.
#'                      Default: all variables.
#' @param horizon       IRF horizon
#' @param n_boot        Bootstrap replications
#' @param ci_level      Confidence level
#' @param plot          Logical; if TRUE, produce plots for each shock
#' @return  Named list of bootstrap_girf() results, one per shock variable
compute_all_irfs <- function(gvar_model, shock_vars = NULL,
                              horizon = 20, n_boot = 500,
                              ci_level = 0.95, plot = TRUE) {

  print_banner("Generalised Impulse Response Functions")

  if (is.null(shock_vars)) {
    shock_vars <- seq_len(gvar_model$k_total)
  }
  if (is.character(shock_vars)) {
    shock_idx <- match(shock_vars, gvar_model$var_names)
  } else {
    shock_idx <- shock_vars
    shock_vars <- gvar_model$var_names[shock_idx]
  }

  irf_list <- setNames(vector("list", length(shock_vars)), shock_vars)

  for (i in seq_along(shock_idx)) {
    sv <- shock_idx[i]
    sname <- shock_vars[i]
    message(sprintf("  Computing GIRF for shock to '%s' ...", sname))

    irf_list[[sname]] <- bootstrap_girf(
      gvar_model, shock_var = sv,
      horizon = horizon, n_boot = n_boot, ci_level = ci_level
    )

    if (plot) {
      plot_irf(irf_list[[sname]], shock_label = sname)
    }
  }

  return(irf_list)
}
