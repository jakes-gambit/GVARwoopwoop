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
    unit_names <- names(data_list)
    TT_data <- nrow(data_list[[unit_names[1]]])

    # Stack the global vector for the last p periods
    alpha_init <- c()
    for (l in 0:(p_global - 1)) {
      t_idx <- TT_data - l
      x_t <- c()
      for (u in unit_names) {
        x_t <- c(x_t, as.numeric(data_list[[u]][t_idx, ]))
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
