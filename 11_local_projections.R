###############################################################################
#  11_local_projections.R  –  Local Projections (Jordà 2005) for the GVAR
#
#  Local Projections (LPs) estimate impulse responses directly by running a
#  separate regression at each horizon h:
#
#    y_{i,t+h} = alpha_h + beta_h * shock_t + gamma_h' * X_{t} + eps_{t+h,h}
#
#  where:
#    - shock_t   is the identified shock variable at time t
#    - X_t       is a vector of control variables (lags of all system variables)
#    - beta_h    gives the impulse response of variable i at horizon h
#
#  Advantages over VAR-based GIRFs:
#    - Robust to model misspecification (no restrictions beyond horizon h)
#    - Handles potential nonlinearities at each horizon independently
#    - Confidence bands are based on OLS standard errors (no bootstrap needed,
#      though bootstrap is also supported here)
#
#  Contents
#  --------
#  1. lp_single()          – LP regression for one response variable, all horizons
#  2. lp_gvar()            – LP for all variables in the GVAR system
#  3. bootstrap_lp()       – Bootstrap confidence bands for LP-IRFs
#  4. plot_lp_irf()        – Plot LP impulse responses with bands
#  5. compare_lp_girf()    – Side-by-side comparison of LP vs VAR-GIRF
###############################################################################


# ───────────────────────────────────────────────────────────────────────────────
# 1.  Local Projection for a Single Response Variable
# ───────────────────────────────────────────────────────────────────────────────

#' Estimate Local Projection impulse responses for one response variable.
#'
#' For each horizon h = 0, 1, ..., H runs the OLS regression:
#'
#'   y_{t+h} = alpha_h + beta_h * shock_t + gamma_h' * Z_t + eps_{t+h,h}
#'
#' where Z_t includes lags 1..n_lags of all control variables (and the shock).
#'
#' @param y          T-vector or T×1 matrix; the response variable.
#' @param shock      T-vector; the shock series (e.g. unit-shock residual or
#'                   level of a variable).
#' @param controls   T × m matrix of control variables (lags will be created
#'                   internally); pass NULL to use no controls beyond the shock.
#' @param n_lags     Number of lags of controls to include (default 1).
#' @param horizon    Maximum forecast horizon H (default 12).
#' @param ci_level   Confidence level for Newey-West HAC standard errors
#'                   (default 0.95).
#' @param nw_lags    Number of lags for Newey-West covariance (default:
#'                   floor(horizon^(2/3)), the data-dependent rule).
#' @return  A data frame with columns:
#'   h, irf, se, lower, upper
lp_single <- function(y, shock, controls = NULL, n_lags = 1,
                      horizon = 12, ci_level = 0.95, nw_lags = NULL) {

  y     <- as.numeric(y)
  shock <- as.numeric(shock)
  TT    <- length(y)

  if (!is.null(controls)) controls <- as.matrix(controls)

  z_alpha <- qnorm(1 - (1 - ci_level) / 2)

  results <- data.frame(
    h     = 0:horizon,
    irf   = NA_real_,
    se    = NA_real_,
    lower = NA_real_,
    upper = NA_real_
  )

  for (h in 0:horizon) {

    # Response: y shifted h periods ahead
    y_h <- c(y[(h + 1):TT], rep(NA, h))

    # Number of NW lags (Andrews 1991 rule: floor(h^(2/3)) or user-supplied)
    nw <- if (!is.null(nw_lags)) nw_lags else max(1L, floor(h^(2/3)))

    # Build regressor matrix: intercept | shock | lags of [shock, controls]
    # We always include n_lags lags of the shock to absorb autocorrelation.
    X_list <- list(shock)
    if (!is.null(controls)) {
      X_list <- c(X_list, lapply(seq_len(ncol(controls)),
                                  function(j) controls[, j]))
    }

    # Stack contemporaneous regressors, then add lags
    X_contemp <- do.call(cbind, X_list)

    # Build lagged block
    lag_blocks <- list()
    for (l in seq_len(n_lags)) {
      shifted <- rbind(matrix(NA, l, ncol(X_contemp)),
                       X_contemp[1:(TT - l), , drop = FALSE])
      lag_blocks[[l]] <- shifted
    }
    X_full <- cbind(1, X_contemp, do.call(cbind, lag_blocks))

    # Trim rows with NAs (from lags and from the h-step lead of y)
    valid <- complete.cases(y_h, X_full)
    y_reg <- y_h[valid]
    X_reg <- X_full[valid, , drop = FALSE]
    n_obs <- sum(valid)

    if (n_obs < ncol(X_reg) + 2) {
      next   # not enough obs for this horizon
    }

    # OLS
    XtX_inv <- tryCatch(solve(crossprod(X_reg)), error = function(e) NULL)
    if (is.null(XtX_inv)) next

    b    <- XtX_inv %*% crossprod(X_reg, y_reg)
    e    <- y_reg - X_reg %*% b
    irf_h <- b[2]   # coefficient on the shock (column 2 after intercept)

    # Newey-West HAC covariance (Bartlett kernel)
    S <- .nw_cov(X_reg, e, nw)
    V <- XtX_inv %*% S %*% XtX_inv
    se_h <- sqrt(max(V[2, 2], 0))

    results$irf[h + 1]   <- irf_h
    results$se[h + 1]    <- se_h
    results$lower[h + 1] <- irf_h - z_alpha * se_h
    results$upper[h + 1] <- irf_h + z_alpha * se_h
  }

  return(results)
}


# Internal: Newey-West HAC covariance of X'e (Bartlett kernel)
.nw_cov <- function(X, e, n_lags) {
  n  <- nrow(X)
  p  <- ncol(X)
  Xe <- X * as.vector(e)   # n × p: elementwise

  # Lag-0 term
  S <- crossprod(Xe) / n

  for (l in seq_len(n_lags)) {
    w     <- 1 - l / (n_lags + 1)   # Bartlett weight
    Gamma <- crossprod(Xe[(l + 1):n, , drop = FALSE],
                       Xe[1:(n - l), , drop = FALSE]) / n
    S <- S + w * (Gamma + t(Gamma))
  }

  return(S)
}


# ───────────────────────────────────────────────────────────────────────────────
# 2.  Local Projections for the Full GVAR System
# ───────────────────────────────────────────────────────────────────────────────

#' Estimate Local Projection IRFs for all variables in the GVAR system.
#'
#' Identifies the shock as a one-standard-deviation innovation to
#' \code{shock_var} using the residuals from the estimated GVAR (i.e. a
#' GIRF-style identification: the shock series is the GVAR residual for
#' \code{shock_var}, scaled by its standard deviation).
#'
#' @param gvar_model  Fitted GVAR model (output of \code{estimate_gvar}).
#' @param shock_var   Character; name of the shocked variable (in
#'                    "UNIT.varname" format, e.g. "DEU.gdp_logdiff").
#' @param data_list   Named list of T × k_i data matrices (same as used for
#'                    estimation).
#' @param n_lags      Number of control lags in the LP regression (default 1).
#' @param horizon     Maximum IRF horizon (default 12).
#' @param ci_level    Confidence level (default 0.95).
#' @param nw_lags     Newey-West lags (NULL = automatic).
#' @return  A named list of data frames (one per response variable), each with
#'   columns h, irf, se, lower, upper.  Also has attributes:
#'   \item{shock_var}{Name of the shocked variable}
#'   \item{var_names}{All response variable names}
lp_gvar <- function(gvar_model, shock_var, data_list,
                    n_lags = 1, horizon = 12, ci_level = 0.95,
                    nw_lags = NULL) {

  print_banner(sprintf("Local Projections  –  shock: %s", shock_var))

  var_names <- gvar_model$var_names
  k_total   <- gvar_model$k_total

  # ── Build global data matrix ──────────────────────────────────────────────
  unit_names <- names(data_list)
  TT <- nrow(data_list[[unit_names[1]]])

  x_global <- matrix(NA, nrow = TT, ncol = k_total)
  colnames(x_global) <- var_names
  for (u in unit_names) {
    mat_u <- as.matrix(data_list[[u]])
    for (v in colnames(mat_u)) {
      gname <- paste0(u, ".", v)
      if (gname %in% var_names) x_global[, gname] <- mat_u[, v]
    }
  }

  # ── Identify shock series ─────────────────────────────────────────────────
  shock_idx <- which(var_names == shock_var)
  if (length(shock_idx) == 0) {
    stop(sprintf("[LP] shock_var '%s' not found in var_names.", shock_var))
  }

  # Use GVAR residuals for the shock variable, scaled to unit s.d.
  resid_global <- gvar_model$residuals   # T_eff × k_total
  T_eff <- nrow(resid_global)

  raw_shock  <- resid_global[, shock_idx]
  shock_sd   <- sd(raw_shock, na.rm = TRUE)
  shock_series_eff <- raw_shock / shock_sd   # unit-s.d. shock

  # Align shock series with x_global (GVAR residuals are T_eff = T - p_global rows)
  p_global <- gvar_model$p_global
  t_start  <- p_global + 1   # first observation used in estimation
  shock_full <- rep(NA, TT)
  shock_full[t_start:TT] <- shock_series_eff

  # ── Controls: all GVAR variables (x_global) ──────────────────────────────
  controls <- x_global

  # ── Run LP for each response variable ────────────────────────────────────
  lp_results <- setNames(vector("list", k_total), var_names)

  for (j in seq_len(k_total)) {
    resp_name <- var_names[j]
    message(sprintf("  LP horizon 0-%d for response '%s' ...", horizon, resp_name))

    lp_j <- lp_single(
      y        = x_global[, j],
      shock    = shock_full,
      controls = controls,
      n_lags   = n_lags,
      horizon  = horizon,
      ci_level = ci_level,
      nw_lags  = nw_lags
    )
    lp_results[[resp_name]] <- lp_j
  }

  attr(lp_results, "shock_var") <- shock_var
  attr(lp_results, "var_names") <- var_names
  attr(lp_results, "ci_level")  <- ci_level

  message(sprintf("[LP] Done. %d response variables, horizons 0-%d.", k_total, horizon))
  return(lp_results)
}


# ───────────────────────────────────────────────────────────────────────────────
# 3.  Bootstrap Confidence Bands
# ───────────────────────────────────────────────────────────────────────────────

#' Bootstrap confidence bands for Local Projection IRFs.
#'
#' Resamples residuals (pairs bootstrap, also called the wild bootstrap) to
#' produce robust confidence intervals that are valid under heteroskedasticity.
#'
#' @param gvar_model  Fitted GVAR model.
#' @param shock_var   Shocked variable name.
#' @param data_list   Raw data list.
#' @param n_lags      LP control lags.
#' @param horizon     Max horizon.
#' @param n_boot      Number of bootstrap replications (default 500).
#' @param ci_level    Confidence level (default 0.95).
#' @param seed        Random seed.
#' @return  A named list of data frames (one per response variable) with columns
#'   h, irf, lower, upper.
bootstrap_lp <- function(gvar_model, shock_var, data_list,
                         n_lags = 1, horizon = 12,
                         n_boot = 500, ci_level = 0.95, seed = 42) {

  print_banner(sprintf("Bootstrap LP  –  shock: %s  (%d reps)", shock_var, n_boot))
  set.seed(seed)

  var_names <- gvar_model$var_names
  k_total   <- gvar_model$k_total
  unit_names <- names(data_list)
  TT <- nrow(data_list[[unit_names[1]]])

  # Build global data matrix (same as lp_gvar)
  x_global <- matrix(NA, nrow = TT, ncol = k_total)
  colnames(x_global) <- var_names
  for (u in unit_names) {
    mat_u <- as.matrix(data_list[[u]])
    for (v in colnames(mat_u)) {
      gname <- paste0(u, ".", v)
      if (gname %in% var_names) x_global[, gname] <- mat_u[, v]
    }
  }

  shock_idx  <- which(var_names == shock_var)
  resid_glb  <- gvar_model$residuals
  T_eff      <- nrow(resid_glb)
  p_global   <- gvar_model$p_global
  t_start    <- p_global + 1

  raw_shock  <- resid_glb[, shock_idx]
  shock_sd   <- sd(raw_shock, na.rm = TRUE)
  shock_full <- rep(NA, TT)
  shock_full[t_start:TT] <- raw_shock / shock_sd

  # Point estimates from the analytic LP
  lp_point <- lp_gvar(gvar_model, shock_var, data_list,
                       n_lags = n_lags, horizon = horizon,
                       ci_level = ci_level)

  # Storage: boot_draws[[resp_var]][ boot_rep, horizon ]
  boot_draws <- lapply(var_names, function(v) {
    matrix(NA, nrow = n_boot, ncol = horizon + 1)
  })
  names(boot_draws) <- var_names

  controls <- x_global

  for (b in seq_len(n_boot)) {
    # Wild bootstrap: multiply shock series by i.i.d. ±1 (Rademacher)
    eta <- sample(c(-1, 1), TT, replace = TRUE)
    shock_b <- shock_full * eta

    for (j in seq_len(k_total)) {
      lp_b <- lp_single(
        y        = x_global[, j],
        shock    = shock_b,
        controls = controls,
        n_lags   = n_lags,
        horizon  = horizon,
        ci_level = ci_level
      )
      boot_draws[[var_names[j]]][b, ] <- lp_b$irf
    }
  }

  alpha <- (1 - ci_level) / 2
  results <- setNames(vector("list", k_total), var_names)

  for (v in var_names) {
    draws_v <- boot_draws[[v]]
    point_v <- lp_point[[v]]$irf
    lo <- apply(draws_v, 2, quantile, probs = alpha,     na.rm = TRUE)
    hi <- apply(draws_v, 2, quantile, probs = 1 - alpha, na.rm = TRUE)
    results[[v]] <- data.frame(
      h     = 0:horizon,
      irf   = point_v,
      lower = lo,
      upper = hi
    )
  }

  attr(results, "shock_var") <- shock_var
  attr(results, "var_names") <- var_names
  attr(results, "ci_level")  <- ci_level

  message(sprintf("[Bootstrap LP] Done (%d reps).", n_boot))
  return(results)
}


# ───────────────────────────────────────────────────────────────────────────────
# 4.  Plot LP Impulse Responses
# ───────────────────────────────────────────────────────────────────────────────

#' Plot Local Projection impulse responses with confidence bands.
#'
#' @param lp_result   Output from \code{lp_gvar()} or \code{bootstrap_lp()}.
#' @param variables   Character vector of response variables to plot (NULL = all).
#' @param max_vars    Maximum number of panels when variables = NULL (default 9).
#' @param shock_label Character; label for the shock (used in plot title).
#' @param ncol_panels Integer; number of columns in the panel layout (default 3).
#' @return  Invisible NULL; plots produced as side effects.
plot_lp_irf <- function(lp_result, variables = NULL, max_vars = 9,
                        shock_label = NULL, ncol_panels = 3) {

  var_names <- attr(lp_result, "var_names")
  shock_var <- attr(lp_result, "shock_var")
  ci_level  <- attr(lp_result, "ci_level")

  if (is.null(shock_label)) shock_label <- shock_var

  if (is.null(variables)) {
    variables <- var_names[seq_len(min(max_vars, length(var_names)))]
  }
  n_plot <- length(variables)

  nrow_panels <- ceiling(n_plot / ncol_panels)
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(nrow_panels, ncol_panels),
      mar   = c(3, 3, 2.5, 1),
      oma   = c(0, 0, 3, 0))

  ci_pct <- round(ci_level * 100)

  for (v in variables) {
    df <- lp_result[[v]]
    if (is.null(df)) next

    ylim <- range(c(df$lower, df$upper), na.rm = TRUE, finite = TRUE)
    if (!all(is.finite(ylim))) ylim <- c(-1, 1)

    plot(df$h, df$irf,
         type = "l", lwd = 2, col = "navy",
         ylim = ylim,
         xlab = "Horizon", ylab = "Response",
         main = v, cex.main = 0.85)

    # Confidence band
    polygon(c(df$h, rev(df$h)),
            c(df$lower, rev(df$upper)),
            col = adjustcolor("steelblue", alpha.f = 0.25),
            border = NA)
    lines(df$h, df$lower, lty = 2, col = "steelblue", lwd = 1)
    lines(df$h, df$upper, lty = 2, col = "steelblue", lwd = 1)

    abline(h = 0, col = "grey50", lty = 3)
  }

  mtext(sprintf("Local Projection IRFs  –  Shock: %s  (%d%% CI)",
                shock_label, ci_pct),
        outer = TRUE, cex = 1.1, font = 2)

  invisible(NULL)
}


# ───────────────────────────────────────────────────────────────────────────────
# 5.  Compare LP vs VAR-based GIRF
# ───────────────────────────────────────────────────────────────────────────────

#' Side-by-side comparison of Local Projection and VAR-based GIRF.
#'
#' @param lp_result    Output from \code{lp_gvar()} or \code{bootstrap_lp()}.
#' @param girf_result  Output from \code{bootstrap_girf()}.
#' @param variables    Character vector of response variables to plot (NULL = first 4).
#' @param shock_label  Character; label for the plot title.
#' @return  Invisible NULL.
compare_lp_girf <- function(lp_result, girf_result, variables = NULL,
                             shock_label = NULL) {

  var_names <- attr(lp_result, "var_names")
  shock_var <- attr(lp_result, "shock_var")
  if (is.null(shock_label)) shock_label <- shock_var

  if (is.null(variables)) {
    variables <- var_names[seq_len(min(4, length(var_names)))]
  }
  n_plot <- length(variables)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(n_plot, 2),
      mar   = c(3, 3, 2.5, 1),
      oma   = c(0, 0, 3, 0))

  for (v in variables) {

    # ── LP panel ──────────────────────────────────────────────────────────
    df_lp <- lp_result[[v]]
    if (!is.null(df_lp)) {
      ylim <- range(c(df_lp$lower, df_lp$upper), na.rm = TRUE, finite = TRUE)
      if (!all(is.finite(ylim))) ylim <- c(-1, 1)
      plot(df_lp$h, df_lp$irf,
           type = "l", lwd = 2, col = "navy",
           ylim = ylim,
           xlab = "Horizon", ylab = "Response",
           main = paste0(v, " (LP)"), cex.main = 0.85)
      polygon(c(df_lp$h, rev(df_lp$h)),
              c(df_lp$lower, rev(df_lp$upper)),
              col = adjustcolor("steelblue", alpha.f = 0.25), border = NA)
      lines(df_lp$h, df_lp$lower, lty = 2, col = "steelblue", lwd = 1)
      lines(df_lp$h, df_lp$upper, lty = 2, col = "steelblue", lwd = 1)
      abline(h = 0, col = "grey50", lty = 3)
    }

    # ── GIRF panel ────────────────────────────────────────────────────────
    if (!is.null(girf_result) && v %in% names(girf_result$irf)) {
      irf_v <- girf_result$irf[[v]]
      lo_v  <- girf_result$lower[[v]]
      hi_v  <- girf_result$upper[[v]]
      hh    <- seq_along(irf_v) - 1

      ylim2 <- range(c(lo_v, hi_v), na.rm = TRUE, finite = TRUE)
      if (!all(is.finite(ylim2))) ylim2 <- c(-1, 1)
      plot(hh, irf_v,
           type = "l", lwd = 2, col = "darkred",
           ylim = ylim2,
           xlab = "Horizon", ylab = "Response",
           main = paste0(v, " (VAR GIRF)"), cex.main = 0.85)
      polygon(c(hh, rev(hh)), c(lo_v, rev(hi_v)),
              col = adjustcolor("salmon", alpha.f = 0.25), border = NA)
      lines(hh, lo_v, lty = 2, col = "salmon", lwd = 1)
      lines(hh, hi_v, lty = 2, col = "salmon", lwd = 1)
      abline(h = 0, col = "grey50", lty = 3)
    } else {
      plot.new()
      text(0.5, 0.5, paste0(v, " (VAR GIRF)\nNot available"),
           cex = 0.9, col = "grey40")
    }
  }

  mtext(sprintf("LP vs VAR GIRF  –  Shock: %s", shock_label),
        outer = TRUE, cex = 1.1, font = 2)

  invisible(NULL)
}
