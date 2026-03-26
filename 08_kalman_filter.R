###############################################################################
#  08_kalman_filter.R  –  Conditional Forecasting via the Kalman Filter
#
#  Conditional forecasting allows the user to impose known future paths on
#  a subset of variables (the "conditions") and produce model-consistent
#  forecasts for the remaining variables.
#
#  The approach follows Doan, Litterman & Sims (1984) / Waggoner & Zha (1999):
#
#    1. Write the GVAR in state-space form.
#    2. Treat the conditions as noisy observations of the state.
#    3. Run the Kalman filter / smoother to obtain the conditional distribution
#       of the unconstrained variables.
#
#  State-Space Representation
#  --------------------------
#  Transition:   α_{t+1}  =  T * α_t  +  c  +  R * η_t      η_t ~ N(0, Q)
#  Observation:  y_t       =  Z * α_t  +  d  +  ε_t          ε_t ~ N(0, H)
#
#  where:
#    α_t  = stacked state  [x_t; x_{t-1}; ... ; x_{t-p+1}]   (kp × 1)
#    T    = companion matrix
#    c    = companion-form intercept
#    R    = selection matrix for shocks  (kp × k)
#    Q    = Σ  (k × k residual covariance)
#    Z    = selection matrix picking out the conditioned variables
#    H    = measurement noise covariance (set very small for hard conditions,
#           or positive for soft conditions)
#
#  Contents
#  --------
#  1. build_state_space()             – set up T, c, R, Q from the GVAR
#  2. build_condition_matrices()      – set up Z, d, H from the conditions
#  3. kalman_filter()                 – forward Kalman filter
#  4. kalman_smoother()               – Rauch-Tung-Striebel backward smoother
#  5. conditional_forecast()          – master function
#  6. plot_conditional_forecast()     – visualisation
#  7. kalman_temporal_disaggregate()  – annual → quarterly via Kalman smoother
###############################################################################


# ───────────────────────────────────────────────────────────────────────────────
# 1.  Build State-Space Matrices from the GVAR
# ───────────────────────────────────────────────────────────────────────────────

#' Construct the state-space transition matrices from the global VAR.
#'
#' @param gvar_model  Fitted GVAR model (from estimate_gvar)
#' @return List with elements: TT (transition), cc (intercept), RR (shock
#'         selection), QQ (state noise covariance), k_total, kp
build_state_space <- function(gvar_model) {

  k_total  <- gvar_model$k_total
  p_global <- gvar_model$p_global
  kp       <- k_total * p_global

  # Transition matrix = companion matrix  (kp × kp)
  TT <- gvar_model$companion

  # Intercept in companion form  (kp × 1)
  cc <- c(gvar_model$F0, rep(0, kp - k_total))

  # Shock selection: shocks hit only the first k_total elements of the state
  RR <- rbind(diag(k_total), matrix(0, nrow = kp - k_total, ncol = k_total))

  # State noise covariance  (k_total × k_total)
  QQ <- gvar_model$Sigma

  return(list(
    TT      = TT,
    cc      = cc,
    RR      = RR,
    QQ      = QQ,
    k_total = k_total,
    kp      = kp
  ))
}


# ───────────────────────────────────────────────────────────────────────────────
# 2.  Build Observation (Condition) Matrices
# ───────────────────────────────────────────────────────────────────────────────

#' Set up the observation equation matrices for the conditioned variables.
#'
#' @param gvar_model  Fitted GVAR model
#' @param conditions  A data frame (or list) with columns:
#'   \describe{
#'     \item{variable}{Character; name of the conditioned variable
#'                     (must match gvar_model$var_names)}
#'     \item{horizon}{Integer; forecast horizon at which the condition applies
#'                    (1 = one step ahead, etc.)}
#'     \item{value}{Numeric; the imposed value}
#'     \item{hardness}{Numeric in (0, Inf]; small values → hard constraint,
#'                     large → soft.  Use 1e-8 for hard constraints.
#'                     Default: 1e-8 for all.}
#'   }
#' @param max_h       Integer; maximum forecast horizon
#' @return A list with one element per horizon h = 1,...,max_h.  Each element
#'         contains Z_h (m_h × kp), d_h (m_h × 1), H_h (m_h × m_h)
#'         where m_h is the number of conditions at horizon h.
#'         If no conditions at horizon h, Z_h = NULL.
build_condition_matrices <- function(gvar_model, conditions, max_h) {

  k_total <- gvar_model$k_total
  kp      <- k_total * gvar_model$p_global
  vnames  <- gvar_model$var_names

  # Ensure conditions is a data frame
  conditions <- as.data.frame(conditions, stringsAsFactors = FALSE)

  if (!"hardness" %in% names(conditions)) {
    conditions$hardness <- 1e-8   # hard constraint by default
  }

  obs_mats <- vector("list", max_h)

  for (h in seq_len(max_h)) {
    cond_h <- conditions[conditions$horizon == h, , drop = FALSE]

    if (nrow(cond_h) == 0) {
      obs_mats[[h]] <- list(Z = NULL, d = NULL, H = NULL)
      next
    }

    m_h <- nrow(cond_h)

    # Z matrix: each row selects one variable from the state vector
    Z_h <- matrix(0, nrow = m_h, ncol = kp)
    d_h <- numeric(m_h)
    H_h <- diag(m_h)

    for (r in seq_len(m_h)) {
      var_name <- cond_h$variable[r]
      var_idx  <- which(vnames == var_name)
      if (length(var_idx) == 0) {
        stop("Conditioned variable '", var_name, "' not found in the GVAR model.")
      }
      # Select the variable from the first k_total positions of the state
      Z_h[r, var_idx] <- 1
      d_h[r]          <- cond_h$value[r]
      H_h[r, r]       <- cond_h$hardness[r]
    }

    obs_mats[[h]] <- list(Z = Z_h, d = d_h, H = H_h)
  }

  return(obs_mats)
}


# ───────────────────────────────────────────────────────────────────────────────
# 3.  Kalman Filter (Forward Pass)
# ───────────────────────────────────────────────────────────────────────────────

#' Run the Kalman filter forward over the forecast horizon.
#'
#' At each step h:
#'   Predict:
#'     α_{h|h-1} = T * α_{h-1|h-1} + c
#'     P_{h|h-1} = T * P_{h-1|h-1} * T' + R * Q * R'
#'
#'   Update (if observations/conditions exist at h):
#'     v_h        = y_h - Z_h * α_{h|h-1} - d_h
#'     S_h        = Z_h * P_{h|h-1} * Z_h' + H_h
#'     K_h        = P_{h|h-1} * Z_h' * S_h^{-1}
#'     α_{h|h}    = α_{h|h-1} + K_h * v_h
#'     P_{h|h}    = (I - K_h * Z_h) * P_{h|h-1}
#'
#' @param ss         State-space matrices from build_state_space()
#' @param obs_mats   List of observation matrices from build_condition_matrices()
#' @param alpha_init Initial state mean  (kp × 1)
#' @param P_init     Initial state covariance (kp × kp)
#' @param max_h      Forecast horizon
#' @return List of filtered states, covariances, and log-likelihood
kalman_filter <- function(ss, obs_mats, alpha_init, P_init, max_h) {

  kp <- ss$kp
  TT <- ss$TT
  cc <- ss$cc
  RR <- ss$RR
  QQ <- ss$QQ

  # Pre-compute R * Q * R'
  RQRT <- RR %*% QQ %*% t(RR)

  # Storage
  alpha_pred <- alpha_filt <- vector("list", max_h)
  P_pred     <- P_filt     <- vector("list", max_h)
  log_lik <- 0

  alpha_prev <- alpha_init
  P_prev     <- P_init

  for (h in seq_len(max_h)) {

    # ---- Prediction step ----
    a_pred <- as.numeric(TT %*% alpha_prev + cc)
    P_pr   <- TT %*% P_prev %*% t(TT) + RQRT

    # Symmetrise (numerical safety)
    P_pr <- (P_pr + t(P_pr)) / 2

    alpha_pred[[h]] <- a_pred
    P_pred[[h]]     <- P_pr

    # ---- Update step (if conditions exist at this horizon) ----
    obs_h <- obs_mats[[h]]

    if (!is.null(obs_h$Z)) {
      Z_h <- obs_h$Z
      d_h <- obs_h$d   # the conditioned VALUES
      H_h <- obs_h$H

      # Innovation
      v_h <- d_h - Z_h %*% a_pred

      # Innovation covariance
      S_h <- Z_h %*% P_pr %*% t(Z_h) + H_h
      S_h <- (S_h + t(S_h)) / 2

      # Kalman gain
      S_inv <- tryCatch(solve(S_h), error = function(e) MASS::ginv(S_h))
      K_h <- P_pr %*% t(Z_h) %*% S_inv

      # Filtered state
      a_filt <- a_pred + as.numeric(K_h %*% v_h)
      P_f    <- (diag(kp) - K_h %*% Z_h) %*% P_pr
      P_f    <- (P_f + t(P_f)) / 2

      # Log-likelihood contribution
      m_h <- length(d_h)
      log_lik <- log_lik - 0.5 * (m_h * log(2 * pi) +
                                    log(det(S_h)) +
                                    as.numeric(t(v_h) %*% S_inv %*% v_h))

    } else {
      # No observation: filtered = predicted
      a_filt <- a_pred
      P_f    <- P_pr
    }

    alpha_filt[[h]] <- a_filt
    P_filt[[h]]     <- P_f

    # Advance
    alpha_prev <- a_filt
    P_prev     <- P_f
  }

  return(list(
    alpha_pred = alpha_pred,
    P_pred     = P_pred,
    alpha_filt = alpha_filt,
    P_filt     = P_filt,
    log_lik    = log_lik
  ))
}


# ───────────────────────────────────────────────────────────────────────────────
# 4.  Kalman Smoother (Backward Pass)
# ───────────────────────────────────────────────────────────────────────────────

#' Rauch-Tung-Striebel smoother for the state-space model.
#'
#' Starting from the last filtered state and working backwards:
#'
#'   J_h       = P_{h|h} * T' * P_{h+1|h}^{-1}
#'   α_{h|H}   = α_{h|h} + J_h * (α_{h+1|H} - α_{h+1|h})
#'   P_{h|H}   = P_{h|h} + J_h * (P_{h+1|H} - P_{h+1|h}) * J_h'
#'
#' @param kf_result   Output from kalman_filter()
#' @param ss          State-space matrices
#' @param max_h       Forecast horizon
#' @return List of smoothed states and covariances
kalman_smoother <- function(kf_result, ss, max_h) {

  kp <- ss$kp
  TT <- ss$TT

  alpha_smooth <- vector("list", max_h)
  P_smooth     <- vector("list", max_h)

  # Initialise at the last time step: smoothed = filtered

  alpha_smooth[[max_h]] <- kf_result$alpha_filt[[max_h]]
  P_smooth[[max_h]]     <- kf_result$P_filt[[max_h]]

  for (h in (max_h - 1):1) {
    P_filt_h     <- kf_result$P_filt[[h]]
    P_pred_h1    <- kf_result$P_pred[[h + 1]]
    alpha_filt_h <- kf_result$alpha_filt[[h]]
    alpha_pred_h1 <- kf_result$alpha_pred[[h + 1]]

    # Smoother gain
    P_pred_inv <- tryCatch(solve(P_pred_h1),
                           error = function(e) MASS::ginv(P_pred_h1))
    J_h <- P_filt_h %*% t(TT) %*% P_pred_inv

    # Smoothed state
    alpha_smooth[[h]] <- alpha_filt_h +
      as.numeric(J_h %*% (alpha_smooth[[h + 1]] - alpha_pred_h1))

    # Smoothed covariance
    P_smooth[[h]] <- P_filt_h +
      J_h %*% (P_smooth[[h + 1]] - P_pred_h1) %*% t(J_h)
    P_smooth[[h]] <- (P_smooth[[h]] + t(P_smooth[[h]])) / 2
  }

  return(list(
    alpha_smooth = alpha_smooth,
    P_smooth     = P_smooth
  ))
}


# ───────────────────────────────────────────────────────────────────────────────
# 5.  Master Conditional Forecasting Function
# ───────────────────────────────────────────────────────────────────────────────

#' Produce conditional forecasts from the GVAR model.
#'
#' The user specifies a set of "conditions" — known future paths for certain
#' variables — and the function returns model-consistent forecasts (with
#' uncertainty) for the remaining (free) variables.
#'
#' @param gvar_model   Fitted GVAR model
#' @param conditions   Data frame with columns: variable, horizon, value,
#'                     (optional) hardness.  See build_condition_matrices().
#' @param max_h        Maximum forecast horizon (default: max of conditions$horizon)
#' @param data_list    Named list of raw data matrices (for initial state)
#' @param ci_level     Confidence level for uncertainty bands (default 0.90)
#' @return A list:
#'   \item{forecast_mean}{(max_h) × k_total  matrix of conditional forecast means}
#'   \item{forecast_lower}{Lower confidence band}
#'   \item{forecast_upper}{Upper confidence band}
#'   \item{forecast_sd}{Standard deviations}
#'   \item{conditions}{The input conditions for reference}
#'   \item{var_names}{Variable names}
conditional_forecast <- function(gvar_model, conditions, max_h = NULL,
                                  data_list = NULL, ci_level = 0.90) {

  print_banner("Conditional Forecasting (Kalman Filter)")

  conditions <- as.data.frame(conditions, stringsAsFactors = FALSE)

  if (is.null(max_h)) max_h <- max(conditions$horizon)

  k_total  <- gvar_model$k_total
  p_global <- gvar_model$p_global
  kp       <- k_total * p_global
  vnames   <- gvar_model$var_names

  # --- Build state-space representation ---
  ss <- build_state_space(gvar_model)

  # --- Build condition matrices ---
  obs_mats <- build_condition_matrices(gvar_model, conditions, max_h)

  # --- Set initial state from the last p observations ---
  if (!is.null(data_list)) {
    TT_data <- nrow(data_list[[names(data_list)[1]]])

    # Build the initial state using gvar_model$var_names to select exactly the
    # variables that belong to the model, in the correct order.
    # var_names entries are "UNIT.varname"; some units (non-dominant) may have
    # fewer columns than data_list[[unit]] because global vars were stripped.
    alpha_init <- c()
    for (l in 0:(p_global - 1)) {
      t_idx <- TT_data - l
      x_t <- numeric(k_total)
      for (v in seq_len(k_total)) {
        vname  <- vnames[v]
        dot    <- regexpr("\\.", vname)[1]
        u      <- substr(vname, 1, dot - 1)
        varstr <- substr(vname, dot + 1, nchar(vname))
        if (u %in% names(data_list) && varstr %in% colnames(data_list[[u]])) {
          x_t[v] <- data_list[[u]][t_idx, varstr]
        }
      }
      alpha_init <- c(alpha_init, x_t)
    }
  } else {
    # Default: zero initial state (less accurate)
    alpha_init <- rep(0, kp)
    message("  [!] No data_list supplied; using zero initial state.")
  }

  # Initial state covariance: use a diffuse prior
  # (large diagonal reflecting uncertainty about the initial state)
  P_init <- diag(kp) * 1e4

  # --- Run the Kalman filter ---
  message("  Running Kalman filter ...")
  kf <- kalman_filter(ss, obs_mats, alpha_init, P_init, max_h)

  # --- Run the smoother ---
  message("  Running Kalman smoother ...")
  ks <- kalman_smoother(kf, ss, max_h)

  # --- Extract forecasts ---
  forecast_mean <- matrix(NA, nrow = max_h, ncol = k_total)
  forecast_sd   <- matrix(NA, nrow = max_h, ncol = k_total)

  for (h in seq_len(max_h)) {
    # Mean: first k_total elements of the smoothed state
    forecast_mean[h, ] <- ks$alpha_smooth[[h]][1:k_total]

    # Standard deviation: sqrt of diagonal of smoothed covariance
    P_h <- ks$P_smooth[[h]]
    forecast_sd[h, ] <- sqrt(pmax(diag(P_h)[1:k_total], 0))
  }

  colnames(forecast_mean) <- colnames(forecast_sd) <- vnames
  rownames(forecast_mean) <- rownames(forecast_sd) <- paste0("h=", 1:max_h)

  # Confidence bands
  z_val <- qnorm(1 - (1 - ci_level) / 2)
  forecast_lower <- forecast_mean - z_val * forecast_sd
  forecast_upper <- forecast_mean + z_val * forecast_sd

  message(sprintf(
    "  Conditional forecast produced for %d horizons, %d variables.",
    max_h, k_total))
  message(sprintf(
    "  Conditions imposed on %d variable-horizon pairs.",
    nrow(conditions)))

  return(list(
    forecast_mean  = forecast_mean,
    forecast_lower = forecast_lower,
    forecast_upper = forecast_upper,
    forecast_sd    = forecast_sd,
    conditions     = conditions,
    var_names      = vnames,
    ci_level       = ci_level,
    max_h          = max_h,
    kf_result      = kf,
    ks_result      = ks
  ))
}


# ───────────────────────────────────────────────────────────────────────────────
# 6.  Plot Conditional Forecasts
# ───────────────────────────────────────────────────────────────────────────────

#' Visualise conditional forecasts with uncertainty bands and conditions.
#'
#' @param cf_result     Output from conditional_forecast()
#' @param plot_vars     Character or integer vector of variables to plot
#' @param history       Optional: data frame or matrix of historical data to
#'                      prepend to the plot (same columns as forecast_mean)
#' @return  A ggplot2 object
plot_conditional_forecast <- function(cf_result, plot_vars = NULL,
                                      history = NULL) {

  fm    <- cf_result$forecast_mean
  fl    <- cf_result$forecast_lower
  fu    <- cf_result$forecast_upper
  conds <- cf_result$conditions
  vnames <- cf_result$var_names

  max_h <- nrow(fm)

  if (is.null(plot_vars)) {
    plot_vars <- seq_len(min(ncol(fm), 9))
  }
  if (is.character(plot_vars)) {
    plot_vars <- match(plot_vars, vnames)
  }

  plot_df <- data.frame()
  for (j in plot_vars) {
    vname <- vnames[j]

    df_j <- data.frame(
      horizon   = 1:max_h,
      mean      = fm[, j],
      lower     = fl[, j],
      upper     = fu[, j],
      variable  = vname,
      type      = "forecast",
      stringsAsFactors = FALSE
    )

    # Mark conditioned horizons
    cond_j <- conds[conds$variable == vname, ]
    if (nrow(cond_j) > 0) {
      cond_points <- data.frame(
        horizon  = cond_j$horizon,
        mean     = cond_j$value,
        lower    = NA,
        upper    = NA,
        variable = vname,
        type     = "condition",
        stringsAsFactors = FALSE
      )
      df_j <- rbind(df_j, cond_points)
    }

    plot_df <- rbind(plot_df, df_j)
  }

  fc_df   <- plot_df[plot_df$type == "forecast", ]
  cond_df <- plot_df[plot_df$type == "condition", ]

  p <- ggplot2::ggplot(fc_df, ggplot2::aes(x = horizon)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         fill = "darkorange", alpha = 0.2) +
    ggplot2::geom_line(ggplot2::aes(y = mean), colour = "darkorange",
                       linewidth = 0.8) +
    ggplot2::facet_wrap(~ variable, scales = "free_y") +
    ggplot2::labs(
      title    = "Conditional Forecast (Kalman Filter)",
      subtitle = paste0(cf_result$ci_level * 100, "% confidence bands"),
      x = "Forecast Horizon", y = "Value"
    ) +
    ggplot2::theme_minimal(base_size = 11)

  # Add condition points
  if (nrow(cond_df) > 0) {
    p <- p + ggplot2::geom_point(
      data = cond_df,
      ggplot2::aes(x = horizon, y = mean),
      colour = "red", size = 3, shape = 18
    )
  }

  print(p)
  return(invisible(p))
}


# ───────────────────────────────────────────────────────────────────────────────
# 7.  Kalman Temporal Disaggregation  (annual → quarterly)
# ───────────────────────────────────────────────────────────────────────────────

#' Disaggregate an annual time series to quarterly frequency using a
#' Kalman filter / smoother.
#'
#' State-space formulation (Chow-Lin in state-space form):
#'
#'   State (quarterly):  x_t = rho * x_{t-1} + sigma_q * eps_t,   eps_t ~ N(0,1)
#'   Observation (annual): Y_y = sum_{q=1}^{4} x_{4(y-1)+q}  + nu_y,  nu_y ~ N(0, sigma_obs^2)
#'
#' When an indicator series (e.g. quarterly trade data) is available it can be
#' supplied via `indicator`; the AR(1) state is then replaced by a regression
#' on the indicator in the observation model.
#'
#' @param annual_values  Numeric vector of length n_years; the annual series
#'                       (e.g. annual GDP level or annual average CPI).
#' @param n_years        Integer; number of annual observations.
#' @param start_year     Integer; first year (used for output labelling).
#' @param aggregation    "sum" (default) or "average": how four quarters
#'                       aggregate to a year.
#' @param rho_init       Starting value for AR(1) parameter (optimised by ML).
#'                       Default 0.9.
#' @param sigma_q_init   Starting value for quarterly innovation s.d.
#' @param obs_noise_frac Variance of annual observation noise as a fraction of
#'                       the annual series variance.  Default 1e-4 (near hard).
#' @return  A list with:
#'   \item{quarterly}{Named numeric vector, length 4 * n_years, labelled
#'                    "YYYY-QN".}
#'   \item{annual_fitted}{Annual aggregates from the smoothed quarterly series.}
#'   \item{params}{Estimated rho and sigma_q (via grid search ML).}
kalman_temporal_disaggregate <- function(annual_values,
                                         n_years        = length(annual_values),
                                         start_year     = 1990L,
                                         aggregation    = c("average", "sum"),
                                         rho_init       = 0.9,
                                         sigma_q_init   = NULL,
                                         obs_noise_frac = 1e-4) {

  aggregation <- match.arg(aggregation)
  n_q <- n_years * 4L   # total quarterly periods

  # ── Helper: run the Kalman filter for given (rho, sigma_q) ──────────────
  kf_disagg <- function(rho, sigma_q) {
    sigma_obs2 <- obs_noise_frac * var(annual_values, na.rm = TRUE)

    # Aggregation row  (1 × 4 or 4 × 4 depending on position within year)
    agg_w <- if (aggregation == "sum") 1 else 0.25

    # Storage
    x_pred  <- numeric(n_q)
    P_pred  <- numeric(n_q)
    x_filt  <- numeric(n_q)
    P_filt  <- numeric(n_q)
    ll      <- 0

    # Initialise at unconditional mean
    x_prev <- mean(annual_values, na.rm = TRUE) /
      if (aggregation == "sum") 4 else 1
    P_prev <- sigma_q^2 / (1 - rho^2 + 1e-12)

    # Accumulate within-year state for observation update
    x_acc <- 0   # running sum (or average) of quarterly states in current year
    P_acc <- 0   # associated variance (sum of P's, ignoring covariances approx)

    year_q <- 0L   # quarter within year counter

    for (t in seq_len(n_q)) {
      year_q <- year_q + 1L

      # Predict
      x_pr <- rho * x_prev
      P_pr <- rho^2 * P_prev + sigma_q^2
      x_pred[t] <- x_pr
      P_pred[t] <- P_pr

      # Accumulate within-year (for the end-of-year observation)
      x_acc <- x_acc + agg_w * x_pr
      P_acc <- P_acc + agg_w^2 * P_pr

      if (year_q == 4L) {
        # Year complete: update with annual observation
        y_idx <- (t - 3L) %/% 4L + 1L
        y_obs <- annual_values[y_idx]

        if (!is.na(y_obs)) {
          innov  <- y_obs - x_acc
          S      <- P_acc + sigma_obs2
          K      <- P_pr * agg_w / S   # Kalman gain for the last quarter only
          x_filt[t] <- x_pr + K * innov
          P_filt[t] <- (1 - K * agg_w) * P_pr
          ll <- ll - 0.5 * (log(2 * pi * S) + innov^2 / S)
        } else {
          x_filt[t] <- x_pr
          P_filt[t] <- P_pr
        }

        # Back-fill the earlier quarters in this year (approx, no smoother here)
        if (t > 3L) {
          for (s in (t - 3L):(t - 1L)) {
            x_filt[s] <- x_pred[s]
            P_filt[s] <- P_pred[s]
          }
        }

        # Reset accumulator
        x_acc  <- 0
        P_acc  <- 0
        year_q <- 0L

      } else {
        x_filt[t] <- x_pr
        P_filt[t] <- P_pr
      }

      x_prev <- x_filt[t]
      P_prev <- P_filt[t]
    }

    list(x_filt = x_filt, P_filt = P_filt,
         x_pred = x_pred, P_pred = P_pred, ll = ll)
  }

  # ── Grid search over rho ─────────────────────────────────────────────────
  if (is.null(sigma_q_init)) {
    sigma_q_init <- sd(annual_values, na.rm = TRUE) / 4
  }
  rho_grid <- seq(0.5, 0.99, by = 0.05)
  best_ll  <- -Inf
  best_rho <- rho_init

  for (rho_try in rho_grid) {
    res <- tryCatch(kf_disagg(rho_try, sigma_q_init), error = function(e) NULL)
    if (!is.null(res) && is.finite(res$ll) && res$ll > best_ll) {
      best_ll  <- res$ll
      best_rho <- rho_try
    }
  }

  # ── Forward filter + backward smoother (RTS) with best rho ──────────────
  kf_res <- kf_disagg(best_rho, sigma_q_init)

  # RTS smoother
  x_sm <- kf_res$x_filt
  P_sm <- kf_res$P_filt

  for (t in (n_q - 1L):1L) {
    if (kf_res$P_pred[t + 1L] > 1e-12) {
      J_t   <- best_rho * kf_res$P_filt[t] / kf_res$P_pred[t + 1L]
      x_sm[t] <- kf_res$x_filt[t] + J_t * (x_sm[t + 1L] - kf_res$x_pred[t + 1L])
      P_sm[t] <- kf_res$P_filt[t] + J_t^2 * (P_sm[t + 1L] - kf_res$P_pred[t + 1L])
    }
  }

  # ── Label the quarterly series ───────────────────────────────────────────
  qtr_labels <- character(n_q)
  idx <- 1L
  for (y in start_year:(start_year + n_years - 1L)) {
    for (q in 1:4) {
      qtr_labels[idx] <- paste0(y, "-Q", q)
      idx <- idx + 1L
    }
  }
  names(x_sm) <- qtr_labels

  # Annual aggregates from smoothed series
  annual_fitted <- numeric(n_years)
  agg_w <- if (aggregation == "sum") 1 else 0.25
  for (y in seq_len(n_years)) {
    idx_q <- ((y - 1L) * 4L + 1L):(y * 4L)
    annual_fitted[y] <- sum(agg_w * x_sm[idx_q])
  }

  list(
    quarterly    = x_sm,
    annual_fitted = annual_fitted,
    params       = list(rho = best_rho, sigma_q = sigma_q_init)
  )
}
