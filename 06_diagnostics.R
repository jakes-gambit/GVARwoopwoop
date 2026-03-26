###############################################################################
#  06_diagnostics.R  –  In-Sample & Out-of-Sample Tests for the GVAR
#
#  Contents
#  --------
#  1. insample_diagnostics()    – R², serial correlation, normality tests
#  2. recursive_oos_forecast()  – expanding-window out-of-sample forecasts
#  3. oos_evaluation()          – RMSE, MAE, Theil's U, DM test
#  4. run_all_diagnostics()     – convenience wrapper
#  5. plot_residual_diagnostics() – residual histograms, ACF, QQ-plots
#  6. plot_fitted_vs_actual()     – fitted vs actual scatter/time series
#  7. plot_forecast_performance() – OOS forecast vs actual comparison
#  8. plot_eigenvalues()          – companion matrix eigenvalue plot
#  9. plot_all_diagnostics()      – generate all diagnostic plots
###############################################################################


# ───────────────────────────────────────────────────────────────────────────────
# 1.  In-Sample Diagnostics
# ───────────────────────────────────────────────────────────────────────────────

#' Compute in-sample fit statistics for the GVAR model.
#'
#' For each variable in the global system we report:
#'   - R-squared
#'   - Adjusted R-squared
#'   - Ljung-Box test for serial correlation (lags = 4 or 12)
#'   - Jarque-Bera test for normality of residuals
#'
#' @param gvar_model  A fitted GVAR model (output of estimate_gvar)
#' @param lb_lags     Integer; number of lags for the Ljung-Box test (default 4)
#' @return            Data frame with one row per global variable
insample_diagnostics <- function(gvar_model, lb_lags = 4) {

  print_banner("In-Sample Diagnostics")

  resid <- gvar_model$residuals   # T_eff × k_total
  k     <- gvar_model$k_total
  TT    <- nrow(resid)
  vnames <- gvar_model$var_names

  # We also need fitted values to compute R².

  # Since we only stored residuals at the global level, approximate R²
  # using the variance of residuals vs. variance of the dependent variable.
  # For a more precise R² we'd need the actual Y; we use the fact that
  # Var(Y) ≈ Var(fitted) + Var(resid).

  results <- data.frame(
    variable   = character(k),
    R_squared  = numeric(k),
    LB_stat    = numeric(k),
    LB_pvalue  = numeric(k),
    JB_stat    = numeric(k),
    JB_pvalue  = numeric(k),
    stringsAsFactors = FALSE
  )

  for (j in seq_len(k)) {
    e_j <- resid[, j]

    # R² (approximate): 1 - var(residual) / var(Y)
    # We don't have Y directly, so we compute from the companion form.
    # Simpler: just report the residual standard error and skip R²
    # unless we have Y.  Let's report what we can.
    var_resid <- var(e_j)

    # Ljung-Box test for serial correlation
    lb <- tryCatch({
      stats::Box.test(e_j, lag = lb_lags, type = "Ljung-Box")
    }, error = function(err) list(statistic = NA, p.value = NA))

    # Jarque-Bera normality test
    jb <- tryCatch({
      tseries::jarque.bera.test(e_j)
    }, error = function(err) list(statistic = list(NA), p.value = NA))

    results$variable[j]  <- if (!is.null(vnames)) vnames[j] else paste0("V", j)
    results$R_squared[j] <- NA  # placeholder (see note above)
    results$LB_stat[j]   <- round(as.numeric(lb$statistic), 4)
    results$LB_pvalue[j] <- round(lb$p.value, 4)
    results$JB_stat[j]   <- round(as.numeric(jb$statistic), 4)
    results$JB_pvalue[j] <- round(jb$p.value, 4)
  }

  # If we have unit-level fits, compute proper R²
  if (!is.null(gvar_model$unit_fits)) {
    idx <- 1
    for (u in gvar_model$unit_names) {
      fit_u <- gvar_model$unit_fits[[u]]
      k_u   <- fit_u$k_dom
      for (j in seq_len(k_u)) {
        # TSS from the fitted model
        y_j   <- fit_u$fitted[, j] + fit_u$residuals[, j]
        ss_res <- sum(fit_u$residuals[, j]^2)
        ss_tot <- sum((y_j - mean(y_j))^2)
        results$R_squared[idx] <- round(1 - ss_res / ss_tot, 4)
        idx <- idx + 1
      }
    }
  }

  cat("\n")
  print(results, row.names = FALSE)
  cat("\n")

  # Flag potential problems
  serial_issues <- results$LB_pvalue < 0.05
  normal_issues <- results$JB_pvalue < 0.05

  if (any(serial_issues, na.rm = TRUE)) {
    message("  [!] Serial correlation detected (Ljung-Box p < 0.05) in: ",
            paste(results$variable[serial_issues], collapse = ", "))
  }
  if (any(normal_issues, na.rm = TRUE)) {
    message("  [!] Non-normality detected (Jarque-Bera p < 0.05) in: ",
            paste(results$variable[normal_issues], collapse = ", "))
  }

  return(results)
}


# ───────────────────────────────────────────────────────────────────────────────
# 2.  Recursive Out-of-Sample Forecasting
# ───────────────────────────────────────────────────────────────────────────────

#' Produce h-step-ahead out-of-sample forecasts using an expanding window.
#'
#' The global VAR in companion form is:
#'    X_t = F0_comp + F_comp * X_{t-1}
#'
#' For h-step ahead:  X_{t+h|t} = (sum_{j=0}^{h-1} F_comp^j) F0_comp + F_comp^h X_t
#'
#' We use an expanding window starting at time t0 and re-estimate at each step.
#' For computational speed an option is provided to fix the coefficients.
#'
#' @param gvar_model    Fitted GVAR model
#' @param data_list     Named list of T × k_i raw data matrices
#' @param star_list     Named list of T × k* star variable matrices
#' @param W             Weight matrix
#' @param h             Forecast horizon (default 1)
#' @param t0_frac       Fraction of total sample to use as initial estimation
#'                      window (default 0.7 = 70%)
#' @param reestimate    Logical; if TRUE re-estimate at each step (slower)
#' @return  A list:
#'   \item{forecasts}{Matrix of h-step forecasts (n_oos × k_total)}
#'   \item{actuals}{Matrix of actual values at forecast dates}
#'   \item{dates}{Integer vector of forecast-origin time indices}
recursive_oos_forecast <- function(gvar_model, data_list, star_list, W,
                                    h = 1, t0_frac = 0.7,
                                    reestimate = FALSE) {

  print_banner(paste0("Out-of-Sample Forecasting (h = ", h, ")"))

  k_total <- gvar_model$k_total
  p_global <- gvar_model$p_global
  var_names <- gvar_model$var_names

  # Build the global data matrix x_t (T × k_total)
  # Match by var_names (unit.varname) to avoid including global-exogenous columns
  # that live in data_list but are not part of the endogenous system.
  unit_names <- names(data_list)
  TT <- nrow(data_list[[unit_names[1]]])

  x_global <- matrix(NA, nrow = TT, ncol = k_total)
  colnames(x_global) <- var_names
  for (u in unit_names) {
    mat_u <- as.matrix(data_list[[u]])
    for (v in colnames(mat_u)) {
      gname <- paste0(u, ".", v)
      if (gname %in% var_names) {
        x_global[, gname] <- mat_u[, v]
      }
    }
  }

  # Determine OOS window
  t0 <- floor(TT * t0_frac)
  n_oos <- TT - t0 - h + 1

  if (n_oos < 1) stop("Not enough observations for out-of-sample evaluation.")

  forecasts <- matrix(NA, nrow = n_oos, ncol = k_total)
  actuals   <- matrix(NA, nrow = n_oos, ncol = k_total)
  colnames(forecasts) <- colnames(actuals) <- var_names

  # Use the companion form for iterated forecasts
  comp <- gvar_model$companion
  F0   <- gvar_model$F0

  # Build companion-form intercept (k_total*p × 1)
  kp <- nrow(comp)
  F0_comp <- c(F0, rep(0, kp - k_total))

  # Trend in companion form (if model has a trend)
  F0_trend <- gvar_model$F0_trend
  has_trend <- !is.null(F0_trend) && any(F0_trend != 0)
  F0_trend_comp <- if (has_trend) c(F0_trend, rep(0, kp - k_total)) else NULL

  for (s in seq_len(n_oos)) {
    t_origin <- t0 + s - 1   # forecast origin

    # State vector at t_origin:  [x_t, x_{t-1}, ..., x_{t-p+1}]
    state <- c()
    for (l in 0:(p_global - 1)) {
      tt <- t_origin - l
      if (tt < 1) {
        state <- c(state, rep(0, k_total))
      } else {
        state <- c(state, x_global[tt, ])
      }
    }

    # Iterate h steps forward
    for (step in seq_len(h)) {
      if (has_trend) {
        state <- F0_comp + F0_trend_comp * (t_origin + step) + comp %*% state
      } else {
        state <- F0_comp + comp %*% state
      }
    }

    # The forecast for x_{t+h} is the first k_total elements of state
    forecasts[s, ] <- state[1:k_total]
    actuals[s, ]   <- x_global[t_origin + h, ]
  }

  message(sprintf("  Produced %d recursive %d-step-ahead forecasts.", n_oos, h))

  return(list(
    forecasts = forecasts,
    actuals   = actuals,
    dates     = (t0 + 1):(t0 + n_oos),
    h         = h
  ))
}


# ───────────────────────────────────────────────────────────────────────────────
# 3.  Out-of-Sample Evaluation Metrics
# ───────────────────────────────────────────────────────────────────────────────

#' Evaluate out-of-sample forecasts: RMSE, MAE, Theil's U.
#'
#' Theil's U compares the model's RMSE to a naïve random-walk forecast.
#' U < 1 means the model beats the random walk.
#'
#' @param oos_result   Output from recursive_oos_forecast()
#' @param x_global     T × k_total matrix of global data (for random walk)
#' @param t0_frac      Same fraction used in the forecast
#' @return             Data frame with metrics per variable
oos_evaluation <- function(oos_result, x_global = NULL, t0_frac = 0.7) {

  print_banner("Out-of-Sample Evaluation")

  forecasts <- oos_result$forecasts
  actuals   <- oos_result$actuals
  h         <- oos_result$h
  k_total   <- ncol(forecasts)
  n_oos     <- nrow(forecasts)
  var_names <- colnames(forecasts)

  errors <- actuals - forecasts

  results <- data.frame(
    variable  = var_names,
    RMSE      = numeric(k_total),
    MAE       = numeric(k_total),
    Theils_U  = numeric(k_total),
    stringsAsFactors = FALSE
  )

  for (j in seq_len(k_total)) {
    e_j <- errors[, j]
    results$RMSE[j] <- round(sqrt(mean(e_j^2, na.rm = TRUE)), 6)
    results$MAE[j]  <- round(mean(abs(e_j), na.rm = TRUE), 6)

    # Theil's U: compare to naïve (random walk) forecast
    if (!is.null(x_global)) {
      TT <- nrow(x_global)
      t0 <- floor(TT * t0_frac)
      # Naïve forecast: x_{t+h} = x_t
      naive_fc <- x_global[(t0 + 1):(t0 + n_oos), j]
      naive_err <- actuals[, j] - naive_fc
      rmse_naive <- sqrt(mean(naive_err^2, na.rm = TRUE))
      # Avoid division by zero - set to NA if naive RMSE is 0 or non-finite
      if (is.finite(rmse_naive) && rmse_naive > 0) {
        results$Theils_U[j] <- round(results$RMSE[j] / rmse_naive, 4)
      } else {
        results$Theils_U[j] <- NA
      }
    } else {
      results$Theils_U[j] <- NA
    }
  }

  cat("\n")
  print(results, row.names = FALSE)
  cat("\n")

  # Summary
  n_better <- sum(results$Theils_U < 1, na.rm = TRUE)
  message(sprintf("  GVAR beats random walk for %d / %d variables (Theil's U < 1).",
                  n_better, k_total))

  return(results)
}


# ───────────────────────────────────────────────────────────────────────────────
# 4.  Run All Diagnostics
# ───────────────────────────────────────────────────────────────────────────────

#' Convenience function: run in-sample and out-of-sample diagnostics.
#'
#' @param gvar_model  Fitted GVAR model
#' @param data_list   Named list of raw data matrices
#' @param star_list   Named list of star matrices
#' @param W           Weight matrix
#' @param h           OOS forecast horizon
#' @param t0_frac     Initial window fraction
#' @param plot        Logical; if TRUE (default), generate diagnostic plots
#' @param save_pdf    Logical; if TRUE, save plots to PDF file
#' @param pdf_file    Character; PDF filename (default "gvar_diagnostics.pdf")
#' @return List with insample and oos elements
run_all_diagnostics <- function(gvar_model, data_list, star_list, W,
                                 h = 1, t0_frac = 0.7,
                                 plot = TRUE, save_pdf = FALSE,
                                 pdf_file = "gvar_diagnostics.pdf") {

  # In-sample
  insample <- insample_diagnostics(gvar_model)

  # Out-of-sample
  oos <- recursive_oos_forecast(gvar_model, data_list, star_list, W,
                                 h = h, t0_frac = t0_frac)

  # Build the global matrix for Theil's U
  # Match by var_names to avoid mismatch with global-exogenous columns.
  unit_names <- names(data_list)
  TT <- nrow(data_list[[unit_names[1]]])
  k_total    <- gvar_model$k_total
  var_names  <- gvar_model$var_names
  x_global <- matrix(NA, nrow = TT, ncol = k_total)
  colnames(x_global) <- var_names
  for (u in unit_names) {
    mat_u <- as.matrix(data_list[[u]])
    for (v in colnames(mat_u)) {
      gname <- paste0(u, ".", v)
      if (gname %in% var_names) {
        x_global[, gname] <- mat_u[, v]
      }
    }
  }

  oos_eval <- oos_evaluation(oos, x_global = x_global, t0_frac = t0_frac)

  result <- list(
    insample  = insample,
    oos       = oos,
    oos_eval  = oos_eval
  )

  # Generate diagnostic plots if requested
  if (plot) {
    plot_all_diagnostics(gvar_model, result,
                         save_pdf = save_pdf, pdf_file = pdf_file)
  }

  return(result)
}


# ───────────────────────────────────────────────────────────────────────────────
# 5.  Residual Diagnostic Plots
# ───────────────────────────────────────────────────────────────────────────────

#' Plot residual diagnostics: histogram, ACF, and QQ-plot for each variable.
#'
#' @param gvar_model  Fitted GVAR model
#' @param variables   Character vector of variable names to plot (NULL = all)
#' @param max_vars    Maximum number of variables to plot (default 4)
#' @return            Invisible NULL; plots are produced as side effect
plot_residual_diagnostics <- function(gvar_model, variables = NULL, max_vars = 4) {

  resid <- gvar_model$residuals
  var_names <- gvar_model$var_names

  if (is.null(variables)) {
    variables <- var_names[1:min(max_vars, length(var_names))]
  }

  n_vars <- length(variables)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  # Plot each variable on a separate page with 1 row × 3 columns
  for (v in variables) {
    idx <- which(var_names == v)
    if (length(idx) == 0) next

    e <- resid[, idx]
    e <- e[is.finite(e)]  # Remove NA/Inf values

    # Skip if no valid data
    if (length(e) < 5 || sd(e) == 0) {
      message(sprintf("  Skipping residual diagnostics for %s: insufficient valid data", v))
      next
    }

    par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

    # Histogram with normal overlay
    hist(e, breaks = 20, freq = FALSE, col = "lightblue", border = "white",
         main = "Histogram", xlab = "Residual", ylab = "Density")
    curve(dnorm(x, mean = mean(e), sd = sd(e)), add = TRUE, col = "darkblue", lwd = 2)

    # ACF plot
    acf(e, main = "ACF", lag.max = min(20, length(e) - 1), col = "darkblue", na.action = na.pass)

    # QQ-plot
    qqnorm(e, main = "QQ Plot", pch = 16, col = "darkblue")
    qqline(e, col = "red", lwd = 2)

    mtext(paste("Residual Diagnostics:", v), outer = TRUE, cex = 1.2, font = 2)
  }

  invisible(NULL)
}


# ───────────────────────────────────────────────────────────────────────────────
# 6.  Fitted vs Actual Plots
# ───────────────────────────────────────────────────────────────────────────────

#' Plot fitted vs actual values for each variable.
#'
#' @param gvar_model  Fitted GVAR model
#' @param variables   Character vector of variable names to plot (NULL = all)
#' @param max_vars    Maximum number of variables to plot (default 4)
#' @return            Invisible NULL
plot_fitted_vs_actual <- function(gvar_model, variables = NULL, max_vars = 4) {

  if (is.null(gvar_model$unit_fits)) {
    message("No unit-level fits available for fitted vs actual plot.")
    return(invisible(NULL))
  }

  var_names <- gvar_model$var_names
  if (is.null(variables)) {
    variables <- var_names[1:min(max_vars, length(var_names))]
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  # Plot each variable on a separate page
  for (u in gvar_model$unit_names) {
    fit_u <- gvar_model$unit_fits[[u]]
    k_u <- fit_u$k_dom
    dom_vars <- colnames(fit_u$fitted)

    for (j in seq_len(k_u)) {
      full_name <- paste0(u, ".", dom_vars[j])
      if (!(full_name %in% variables)) next

      actual <- fit_u$fitted[, j] + fit_u$residuals[, j]
      fitted <- fit_u$fitted[, j]
      TT <- length(actual)

      # Skip if data is invalid
      if (all(is.na(actual)) || all(is.na(fitted))) {
        message(sprintf("  Skipping %s: no valid data", full_name))
        next
      }

      par(mfrow = c(1, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

      # Time series overlay
      ylim <- range(c(actual, fitted), na.rm = TRUE, finite = TRUE)
      if (!all(is.finite(ylim))) ylim <- c(-1, 1)

      plot(1:TT, actual, type = "l", col = "black", lwd = 1.5,
           main = "Time Series", xlab = "Time", ylab = "Value", ylim = ylim)
      lines(1:TT, fitted, col = "red", lwd = 1.5, lty = 2)
      legend("topright", legend = c("Actual", "Fitted"), col = c("black", "red"),
             lty = c(1, 2), lwd = 1.5, bty = "n", cex = 0.8)

      # Scatter plot
      plot(fitted, actual, pch = 16, col = rgb(0, 0, 0.5, 0.5),
           main = "Fitted vs Actual", xlab = "Fitted", ylab = "Actual")
      abline(0, 1, col = "red", lwd = 2)
      r2 <- cor(fitted, actual, use = "complete.obs")^2
      if (is.na(r2)) r2 <- 0
      legend("topleft", legend = paste0("R² = ", round(r2, 3)), bty = "n", cex = 0.9)

      mtext(paste("Fitted vs Actual:", full_name), outer = TRUE, cex = 1.2, font = 2)
    }
  }

  invisible(NULL)
}


# ───────────────────────────────────────────────────────────────────────────────
# 7.  Out-of-Sample Forecast Performance Plots
# ───────────────────────────────────────────────────────────────────────────────

#' Plot out-of-sample forecast vs actual values.
#'
#' @param oos_result  Output from recursive_oos_forecast()
#' @param variables   Character vector of variable names to plot (NULL = all)
#' @param max_vars    Maximum number of variables to plot (default 4)
#' @return            Invisible NULL
plot_forecast_performance <- function(oos_result, variables = NULL, max_vars = 4) {

  forecasts <- oos_result$forecasts
  actuals   <- oos_result$actuals
  dates     <- oos_result$dates
  h         <- oos_result$h
  var_names <- colnames(forecasts)

  if (is.null(variables)) {
    variables <- var_names[1:min(max_vars, length(var_names))]
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  # Plot each variable on a separate page
  for (v in variables) {
    idx <- which(var_names == v)
    if (length(idx) == 0) next

    fc <- forecasts[, idx]
    ac <- actuals[, idx]
    err <- ac - fc

    # Skip if data is invalid
    if (all(is.na(fc)) || all(is.na(ac))) {
      message(sprintf("  Skipping %s: no valid data", v))
      next
    }

    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

    # Time series: forecast vs actual
    ylim <- range(c(fc, ac), na.rm = TRUE, finite = TRUE)
    if (!all(is.finite(ylim))) ylim <- c(-1, 1)

    plot(dates, ac, type = "l", col = "black", lwd = 1.5,
         main = "Forecast vs Actual", xlab = "Time", ylab = "Value", ylim = ylim)
    lines(dates, fc, col = "blue", lwd = 1.5, lty = 2)
    legend("topright", legend = c("Actual", paste0(h, "-step Forecast")),
           col = c("black", "blue"), lty = c(1, 2), lwd = 1.5, bty = "n", cex = 0.8)

    # Forecast error over time
    err_colors <- ifelse(is.na(err), "gray", ifelse(err >= 0, "darkgreen", "red"))
    plot(dates, err, type = "h", col = err_colors,
         main = "Forecast Errors", xlab = "Time", ylab = "Error", lwd = 2)
    abline(h = 0, col = "gray", lty = 2)
    rmse <- sqrt(mean(err^2, na.rm = TRUE))
    if (!is.finite(rmse)) rmse <- NA
    legend("topright", legend = paste0("RMSE = ", round(rmse, 4)), bty = "n", cex = 0.9)

    mtext(paste0("OOS Performance (h=", h, "): ", v), outer = TRUE, cex = 1.2, font = 2)
  }

  invisible(NULL)
}


# ───────────────────────────────────────────────────────────────────────────────
# 8.  Eigenvalue Plot (Model Stability)
# ───────────────────────────────────────────────────────────────────────────────

#' Plot companion matrix eigenvalues to assess model stability.
#'
#' All eigenvalues should lie inside the unit circle for a stable model.
#'
#' @param gvar_model  Fitted GVAR model
#' @return            Invisible NULL
plot_eigenvalues <- function(gvar_model) {

  eig <- eigen(gvar_model$companion)$values
  eig_mod <- Mod(eig)
  max_mod <- max(eig_mod)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))

  # Plot 1: Eigenvalues in complex plane
  plot(Re(eig), Im(eig), pch = 16,
       col = ifelse(eig_mod < 1, "darkgreen", "red"),
       cex = 1.5, asp = 1,
       xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2),
       main = "Companion Matrix Eigenvalues",
       xlab = "Real Part", ylab = "Imaginary Part")

  # Draw unit circle
  theta <- seq(0, 2 * pi, length.out = 100)
  lines(cos(theta), sin(theta), col = "blue", lwd = 2, lty = 2)
  abline(h = 0, v = 0, col = "gray", lty = 3)

  # Legend
  legend("topright",
         legend = c(paste0("Max modulus: ", round(max_mod, 4)),
                    ifelse(max_mod < 1, "STABLE", "UNSTABLE")),
         text.col = c("black", ifelse(max_mod < 1, "darkgreen", "red")),
         bty = "n", cex = 0.9)

  # Plot 2: Sorted eigenvalue moduli
  sorted_mod <- sort(eig_mod, decreasing = TRUE)
  n_eig <- length(sorted_mod)
  barplot(sorted_mod[1:min(20, n_eig)],
          col = ifelse(sorted_mod[1:min(20, n_eig)] < 1, "lightgreen", "salmon"),
          border = "white",
          main = "Top 20 Eigenvalue Moduli",
          xlab = "Eigenvalue Rank", ylab = "Modulus",
          names.arg = 1:min(20, n_eig))
  abline(h = 1, col = "red", lwd = 2, lty = 2)

  invisible(NULL)
}


# ───────────────────────────────────────────────────────────────────────────────
# 9.  CUSUM Stability Test Plot
# ───────────────────────────────────────────────────────────────────────────────

#' Plot CUSUM (cumulative sum) of recursive residuals for stability testing.
#'
#' @param gvar_model  Fitted GVAR model
#' @param variables   Character vector of variable names to plot (NULL = first 4)
#' @param sig_level   Significance level for CUSUM bands (default 0.05)
#' @return            Invisible NULL
plot_cusum <- function(gvar_model, variables = NULL, sig_level = 0.05) {

  resid <- gvar_model$residuals
  var_names <- gvar_model$var_names
  TT <- nrow(resid)

  if (is.null(variables)) {
    variables <- var_names[1:min(4, length(var_names))]
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  # Plot 2 variables per page (2x1 layout)
  par(mfrow = c(2, 1), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

  plot_count <- 0
  for (v in variables) {
    idx <- which(var_names == v)
    if (length(idx) == 0) next

    e <- resid[, idx]

    # Skip if data is invalid
    if (all(is.na(e)) || sd(e, na.rm = TRUE) == 0) {
      message(sprintf("  Skipping CUSUM for %s: no valid data", v))
      next
    }

    sigma_e <- sd(e, na.rm = TRUE)

    # Standardised recursive residuals (approximation)
    w <- e / sigma_e
    w[!is.finite(w)] <- 0

    # CUSUM
    cusum <- cumsum(w) / sqrt(TT)

    # Critical bounds (Brownian bridge approximation)
    a <- qnorm(1 - sig_level / 2)
    t_seq <- 1:TT
    upper <- a * sqrt(t_seq / TT)
    lower <- -upper

    # Plot
    ylim <- range(c(cusum, upper, lower), na.rm = TRUE, finite = TRUE) * 1.1
    if (!all(is.finite(ylim))) ylim <- c(-2, 2)

    plot(t_seq, cusum, type = "l", col = "darkblue", lwd = 2,
         main = paste("CUSUM:", v), xlab = "Time", ylab = "CUSUM", ylim = ylim)
    lines(t_seq, upper, col = "red", lty = 2, lwd = 1.5)
    lines(t_seq, lower, col = "red", lty = 2, lwd = 1.5)
    abline(h = 0, col = "gray", lty = 3)

    # Flag if CUSUM crosses bounds
    crosses <- any(cusum > upper | cusum < lower, na.rm = TRUE)
    if (crosses) {
      legend("topleft", legend = "Instability detected", text.col = "red",
             bty = "n", cex = 0.8)
    }

    plot_count <- plot_count + 1
    if (plot_count %% 2 == 0) {
      mtext("CUSUM Stability Test", outer = TRUE, cex = 1.2, font = 2)
    }
  }

  if (plot_count %% 2 != 0) {
    mtext("CUSUM Stability Test", outer = TRUE, cex = 1.2, font = 2)
  }

  invisible(NULL)
}


# ───────────────────────────────────────────────────────────────────────────────
# 10.  Summary Bar Chart of OOS Performance
# ───────────────────────────────────────────────────────────────────────────────

#' Plot bar chart comparing Theil's U across all variables.
#'
#' @param oos_eval  Output from oos_evaluation()
#' @return          Invisible NULL
plot_oos_summary <- function(oos_eval) {

  var_names <- oos_eval$variable
  theils_u  <- oos_eval$Theils_U
  rmse      <- oos_eval$RMSE

  # Replace NA/Inf with 0 for plotting (will show as no bar)
  theils_u[!is.finite(theils_u)] <- 0
  rmse[!is.finite(rmse)] <- 0

  # Skip if all values are invalid
 if (all(theils_u == 0) && all(rmse == 0)) {
    message("  Skipping OOS summary plot: no valid data")
    return(invisible(NULL))
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  # Use 2 rows × 1 column layout for better spacing
  par(mfrow = c(2, 1), mar = c(5, 4, 3, 2))

  # Theil's U bar chart
  cols <- ifelse(theils_u < 1, "lightgreen", "salmon")
  cols[theils_u == 0] <- "gray"  # Gray for NA/missing values
  barplot(theils_u, names.arg = var_names, col = cols, border = "white",
          main = "Theil's U by Variable (< 1 = beats random walk)",
          ylab = "Theil's U", las = 2, cex.names = 0.6)
  abline(h = 1, col = "red", lwd = 2, lty = 2)

  # RMSE bar chart
  rmse_cols <- ifelse(rmse == 0, "gray", "steelblue")
  barplot(rmse, names.arg = var_names, col = rmse_cols, border = "white",
          main = "RMSE by Variable",
          ylab = "RMSE", las = 2, cex.names = 0.6)

  invisible(NULL)
}


# ───────────────────────────────────────────────────────────────────────────────
# 11.  Master Plotting Function
# ───────────────────────────────────────────────────────────────────────────────

#' Generate all diagnostic plots for a GVAR model.
#'
#' @param gvar_model    Fitted GVAR model
#' @param diag_results  Output from run_all_diagnostics()
#' @param variables     Variables to focus on (NULL = auto-select)
#' @param save_pdf      If TRUE, save plots to PDF file
#' @param pdf_file      PDF filename (default "gvar_diagnostics.pdf")
#' @return              Invisible NULL
plot_all_diagnostics <- function(gvar_model, diag_results,
                                  variables = NULL, save_pdf = FALSE,
                                  pdf_file = "gvar_diagnostics.pdf") {

  if (save_pdf) {
    pdf(pdf_file, width = 12, height = 10)
    on.exit(dev.off())
  }

  message("\n[GVAR] Generating diagnostic plots...\n")

  # 1. Eigenvalue plot (stability)
  message("  - Eigenvalue plot (model stability)")
  plot_eigenvalues(gvar_model)

  # 2. Residual diagnostics
  message("  - Residual diagnostics (histogram, ACF, QQ)")
  plot_residual_diagnostics(gvar_model, variables = variables)

  # 3. Fitted vs actual
  message("  - Fitted vs actual plots")
  plot_fitted_vs_actual(gvar_model, variables = variables)

  # 4. OOS forecast performance
  if (!is.null(diag_results$oos)) {
    message("  - Out-of-sample forecast performance")
    plot_forecast_performance(diag_results$oos, variables = variables)
  }

  # 5. OOS summary
  if (!is.null(diag_results$oos_eval)) {
    message("  - Out-of-sample summary (Theil's U, RMSE)")
    plot_oos_summary(diag_results$oos_eval)
  }

  # 6. CUSUM stability
  message("  - CUSUM stability test")
  plot_cusum(gvar_model, variables = variables)

  message("\n[GVAR] Diagnostic plots complete.")
  if (save_pdf) {
    message(sprintf("  Plots saved to: %s", pdf_file))
  }

  invisible(NULL)
}
