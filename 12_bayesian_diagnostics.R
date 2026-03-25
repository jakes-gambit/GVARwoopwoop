###############################################################################
#  12_bayesian_diagnostics.R  –  Bayesian Diagnostics and Plotting
#
#  Diagnostic tools for the Bayesian GVAR:
#    1. DIC (Deviance Information Criterion)
#    2. Log marginal likelihood (closed-form for Normal-IW)
#    3. Posterior predictive checks
#    4. Plotting: Bayesian IRFs, conditional forecasts, prior vs posterior
#
#  Contents
#  --------
#  1. bayesian_diagnostics()     – DIC, marginal likelihood, PPC
#  2. plot_bayesian_irf()        – GIRF plot with credible intervals
#  3. plot_bayesian_forecast()   – Conditional forecast plot (Bayesian)
#  4. plot_prior_posterior()     – Prior vs posterior density comparison
#  5. bayesian_model_comparison() – Compare Bayesian vs frequentist
###############################################################################


# ───────────────────────────────────────────────────────────────────────────────
# 1.  Bayesian Diagnostics
# ───────────────────────────────────────────────────────────────────────────────

#' Compute DIC, log marginal likelihood, and posterior predictive checks.
#'
#' @param bayesian_gvar  Output from bayesian_estimate_gvar()
#' @param gvar_data      Output from prepare_gvar_dataset()
#' @return List with dic, log_marginal_lik, effective_params, ppc_summary
bayesian_diagnostics <- function(bayesian_gvar, gvar_data) {

  print_banner("Bayesian Diagnostics")

  unit_names <- bayesian_gvar$point_estimate$unit_names
  N <- length(unit_names)
  deterministic <- if (!is.null(bayesian_gvar$deterministic)) bayesian_gvar$deterministic else "intercept"

  # ---- DIC (Deviance Information Criterion) ----
  # D(θ) = -2 * Σ_i log L(Y_i | B_i, Σ_i)
  # DIC = D_bar + p_D = 2 * D_bar - D(θ_bar)

  message("  Computing DIC ...")

  # D at posterior mean (θ_bar)
  D_bar_theta <- 0
  for (u in unit_names) {
    cd <- gvar_data$unit_data[[u]]
    Y  <- cd$Y
    det <- build_deterministic_columns(nrow(Y), deterministic)
    X  <- if (!is.null(det$D_mat)) cbind(det$D_mat, cd$X) else cd$X
    post <- bayesian_gvar$unit_posteriors[[u]]
    B_n <- post$B_n
    k_dom <- post$k_dom
    TT <- nrow(Y)

    Sigma_bar <- post$S_n / max(post$v_n - k_dom - 1, 1)
    resid <- Y - X %*% B_n

    # Log-likelihood for MN: -T/2 * log|Σ| - 1/2 * tr(Σ^{-1} (Y-XB)'(Y-XB))
    Sigma_inv <- tryCatch(solve(Sigma_bar), error = function(e) MASS::ginv(Sigma_bar))
    log_det <- determinant(Sigma_bar, logarithm = TRUE)$modulus[1]

    loglik <- -TT / 2 * k_dom * log(2 * pi) - TT / 2 * log_det -
      0.5 * sum(diag(Sigma_inv %*% crossprod(resid)))

    D_bar_theta <- D_bar_theta + (-2 * loglik)
  }

  # D_bar: mean of D(θ^(s)) over draws
  # Use a subset of stable draws for efficiency
  stable_idx <- which(bayesian_gvar$draws$is_stable)
  n_use <- min(200, length(stable_idx))
  if (n_use < 10) {
    message("  [Warning] Fewer than 10 stable draws for DIC computation.")
    dic <- NA
    p_D <- NA
  } else {
    use_idx <- sample(stable_idx, n_use)
    D_draws <- numeric(n_use)

    for (i in seq_len(n_use)) {
      s <- use_idx[i]
      D_s <- 0

      for (u in unit_names) {
        cd <- gvar_data$unit_data[[u]]
        Y  <- cd$Y
        det_u <- build_deterministic_columns(nrow(Y), deterministic)
        X  <- if (!is.null(det_u$D_mat)) cbind(det_u$D_mat, cd$X) else cd$X
        post <- bayesian_gvar$unit_posteriors[[u]]
        k_dom <- post$k_dom
        TT <- nrow(Y)

        # Draw for this unit
        draw <- draw_niw_posterior(post$B_n, post$V_n, post$S_n, post$v_n, k_dom)
        resid <- Y - X %*% draw$B_draw

        Sigma_inv <- tryCatch(solve(draw$Sigma_draw),
                               error = function(e) MASS::ginv(draw$Sigma_draw))
        log_det <- determinant(draw$Sigma_draw, logarithm = TRUE)$modulus[1]

        loglik_s <- -TT / 2 * k_dom * log(2 * pi) - TT / 2 * log_det -
          0.5 * sum(diag(Sigma_inv %*% crossprod(resid)))

        D_s <- D_s + (-2 * loglik_s)
      }

      D_draws[i] <- D_s
    }

    D_bar <- mean(D_draws)
    p_D <- D_bar - D_bar_theta
    dic <- D_bar + p_D   # = 2 * D_bar - D(θ_bar)

    message(sprintf("  DIC = %.2f  (D_bar = %.2f, p_D = %.2f)", dic, D_bar, p_D))
  }

  # ---- Log Marginal Likelihood (per unit, closed-form for Normal-IW) ----
  message("  Computing log marginal likelihood ...")

  log_ml_total <- 0

  for (u in unit_names) {
    post <- bayesian_gvar$unit_posteriors[[u]]
    prior <- post$prior

    k_dom <- post$k_dom
    cd    <- gvar_data$unit_data[[u]]
    TT    <- nrow(cd$Y)
    m     <- nrow(post$B_n)

    v_0 <- prior$v_0
    v_n <- post$v_n
    S_0 <- prior$S_0
    S_n <- post$S_n
    V_0 <- prior$V_0
    V_n <- post$V_n

    # log p(Y) = const + (k/2)*log|V_n|/|V_0|
    #          + (v_0/2)*log|S_0| - (v_n/2)*log|S_n|
    #          + log Gamma_k(v_n/2) - log Gamma_k(v_0/2)
    #          - (T*k/2)*log(pi)

    log_det_V0 <- determinant(V_0, logarithm = TRUE)$modulus[1]
    log_det_Vn <- determinant(V_n, logarithm = TRUE)$modulus[1]
    log_det_S0 <- determinant(S_0, logarithm = TRUE)$modulus[1]
    log_det_Sn <- determinant(S_n, logarithm = TRUE)$modulus[1]

    # Multivariate log gamma function
    mvlgamma <- function(p, a) {
      p * (p - 1) / 4 * log(pi) + sum(lgamma(a + (1 - (1:p)) / 2))
    }

    log_ml_u <- -TT * k_dom / 2 * log(pi) +
      k_dom / 2 * (log_det_Vn - log_det_V0) +
      v_0 / 2 * log_det_S0 - v_n / 2 * log_det_Sn +
      mvlgamma(k_dom, v_n / 2) - mvlgamma(k_dom, v_0 / 2)

    log_ml_total <- log_ml_total + log_ml_u
  }

  message(sprintf("  Log marginal likelihood = %.2f", log_ml_total))

  # ---- Posterior Predictive Checks ----
  message("  Running posterior predictive checks ...")

  ppc_rows <- list()

  for (u in unit_names) {
    cd <- gvar_data$unit_data[[u]]
    Y_actual <- cd$Y
    det_u <- build_deterministic_columns(nrow(cd$Y), deterministic)
    X <- if (!is.null(det_u$D_mat)) cbind(det_u$D_mat, cd$X) else cd$X
    post <- bayesian_gvar$unit_posteriors[[u]]
    k_dom <- post$k_dom
    TT <- nrow(Y_actual)

    # Simulate from posterior predictive
    n_ppc <- 100
    Y_rep_mean <- matrix(0, nrow = TT, ncol = k_dom)
    Y_rep_var  <- matrix(0, nrow = TT, ncol = k_dom)

    for (r in seq_len(n_ppc)) {
      draw <- draw_niw_posterior(post$B_n, post$V_n, post$S_n, post$v_n, k_dom)
      Y_hat <- X %*% draw$B_draw

      # Add noise from Σ
      L_Sigma <- tryCatch(chol(draw$Sigma_draw), error = function(e) {
        eig <- eigen(draw$Sigma_draw, symmetric = TRUE)
        eig$values <- pmax(eig$values, 1e-10)
        chol(eig$vectors %*% diag(eig$values) %*% t(eig$vectors))
      })

      noise <- matrix(stats::rnorm(TT * k_dom), nrow = TT) %*% L_Sigma
      Y_rep <- Y_hat + noise

      Y_rep_mean <- Y_rep_mean + Y_rep / n_ppc
      Y_rep_var  <- Y_rep_var + Y_rep^2 / n_ppc
    }
    Y_rep_var <- Y_rep_var - Y_rep_mean^2

    var_names_u <- colnames(Y_actual)
    if (is.null(var_names_u)) var_names_u <- paste0("V", 1:k_dom)

    for (v in seq_len(k_dom)) {
      ppc_rows[[length(ppc_rows) + 1]] <- data.frame(
        unit         = u,
        variable     = var_names_u[v],
        actual_mean  = round(mean(Y_actual[, v]), 4),
        ppc_mean     = round(mean(Y_rep_mean[, v]), 4),
        actual_sd    = round(sd(Y_actual[, v]), 4),
        ppc_sd       = round(sqrt(mean(Y_rep_var[, v])), 4),
        stringsAsFactors = FALSE
      )
    }
  }

  ppc_summary <- do.call(rbind, ppc_rows)
  rownames(ppc_summary) <- NULL

  message("\n  Posterior Predictive Check Summary:")
  print(ppc_summary, row.names = FALSE)

  return(list(
    dic              = dic,
    p_D              = p_D,
    log_marginal_lik = log_ml_total,
    ppc_summary      = ppc_summary
  ))
}


# ───────────────────────────────────────────────────────────────────────────────
# 2.  Plot Bayesian IRFs
# ───────────────────────────────────────────────────────────────────────────────

#' Plot Bayesian GIRFs with posterior credible intervals.
#'
#' @param bayesian_irf  Output from bayesian_girf()
#' @param response_vars Variables to plot (NULL = all)
#' @param shock_label   Label for the shock
#' @return              A ggplot2 object (also printed)
plot_bayesian_irf <- function(bayesian_irf, response_vars = NULL,
                               shock_label = "Shock") {

  var_names <- bayesian_irf$var_names
  horizon   <- nrow(bayesian_irf$point) - 1
  ci_label  <- paste0(bayesian_irf$ci_level * 100, "% Posterior Credible Interval")

  if (is.null(response_vars)) response_vars <- var_names

  # Build long-format data frame
  df_list <- list()
  for (v in response_vars) {
    v_idx <- match(v, var_names)
    if (is.na(v_idx)) next

    df_list[[length(df_list) + 1]] <- data.frame(
      horizon  = 0:horizon,
      variable = v,
      point    = bayesian_irf$point[, v_idx],
      median   = bayesian_irf$median[, v_idx],
      lower    = bayesian_irf$lower[, v_idx],
      upper    = bayesian_irf$upper[, v_idx],
      stringsAsFactors = FALSE
    )
  }
  df <- do.call(rbind, df_list)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = horizon)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                          fill = "steelblue", alpha = 0.25) +
    ggplot2::geom_line(ggplot2::aes(y = median), color = "steelblue", linewidth = 1) +
    ggplot2::geom_line(ggplot2::aes(y = point), color = "darkred",
                        linewidth = 0.7, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
    ggplot2::facet_wrap(~ variable, scales = "free_y") +
    ggplot2::labs(
      title    = paste0("Bayesian GIRF: ", shock_label),
      subtitle = paste0(ci_label, " | ", bayesian_irf$n_draws, " posterior draws"),
      x = "Horizon", y = "Response"
    ) +
    ggplot2::theme_minimal()

  print(p)
  invisible(p)
}


# ───────────────────────────────────────────────────────────────────────────────
# 3.  Plot Bayesian Conditional Forecast
# ───────────────────────────────────────────────────────────────────────────────

#' Plot Bayesian conditional forecasts with credible intervals.
#'
#' @param bayesian_cf  Output from bayesian_conditional_forecast()
#' @param plot_vars    Variables to plot (NULL = all)
#' @return             A ggplot2 object (also printed)
plot_bayesian_forecast <- function(bayesian_cf, plot_vars = NULL) {

  var_names <- bayesian_cf$var_names
  max_h     <- nrow(bayesian_cf$forecast_mean)
  ci_label  <- paste0(bayesian_cf$ci_level * 100, "% Posterior Credible Interval")

  if (is.null(plot_vars)) plot_vars <- var_names

  # Build long-format data frame
  df_list <- list()
  for (v in plot_vars) {
    v_idx <- match(v, var_names)
    if (is.na(v_idx)) next

    df_v <- data.frame(
      horizon  = 1:max_h,
      variable = v,
      mean     = bayesian_cf$forecast_mean[, v_idx],
      lower    = bayesian_cf$forecast_lower[, v_idx],
      upper    = bayesian_cf$forecast_upper[, v_idx],
      stringsAsFactors = FALSE
    )

    # Mark conditions
    conds <- bayesian_cf$conditions
    conds_v <- conds[conds$variable == v, ]
    if (nrow(conds_v) > 0) {
      df_v$condition_value <- NA
      for (r in seq_len(nrow(conds_v))) {
        h_idx <- which(df_v$horizon == conds_v$horizon[r])
        if (length(h_idx) > 0) {
          df_v$condition_value[h_idx] <- conds_v$value[r]
        }
      }
    } else {
      df_v$condition_value <- NA
    }

    df_list[[length(df_list) + 1]] <- df_v
  }
  df <- do.call(rbind, df_list)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = horizon)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                          fill = "steelblue", alpha = 0.25) +
    ggplot2::geom_line(ggplot2::aes(y = mean), color = "steelblue", linewidth = 1) +
    ggplot2::geom_point(ggplot2::aes(y = condition_value), color = "red",
                         size = 3, na.rm = TRUE) +
    ggplot2::facet_wrap(~ variable, scales = "free_y") +
    ggplot2::labs(
      title    = paste0("Bayesian Conditional Forecast (", bayesian_cf$method, ")"),
      subtitle = paste0(ci_label, " | ", bayesian_cf$n_draws_used, " posterior draws"),
      x = "Horizon", y = "Forecast"
    ) +
    ggplot2::theme_minimal()

  print(p)
  invisible(p)
}


# ───────────────────────────────────────────────────────────────────────────────
# 4.  Prior vs Posterior Density Plot
# ───────────────────────────────────────────────────────────────────────────────

#' Plot prior vs posterior density for selected coefficients of a unit.
#'
#' @param bayesian_gvar  Output from bayesian_estimate_gvar()
#' @param unit_name      Name of the unit to inspect
#' @param coef_indices   Integer vector of coefficient row indices to plot
#'                       (rows of B_n). If NULL, plots first 4 diagonal elements.
#' @param n_draws        Number of draws for density estimation
#' @return               A ggplot2 object (also printed)
plot_prior_posterior <- function(bayesian_gvar, unit_name, coef_indices = NULL,
                                 n_draws = 2000) {

  post <- bayesian_gvar$unit_posteriors[[unit_name]]
  if (is.null(post)) stop("Unit '", unit_name, "' not found.")

  k_dom <- post$k_dom
  m     <- nrow(post$B_n)

  # Default: plot the own-first-lag coefficients (diagonal of A_1)
  if (is.null(coef_indices)) {
    lag1_start <- 1 + post$k_star + post$k_global + 1
    coef_indices <- lag1_start + (seq_len(min(4, k_dom)) - 1)
    coef_indices <- coef_indices[coef_indices <= m]
  }

  # Generate prior and posterior draws
  df_list <- list()

  for (ci in coef_indices) {
    for (eq in seq_len(min(2, k_dom))) {
      # Prior: N(B_0[ci, eq], V_0[ci, ci] * sigma^2)
      prior_mean <- post$prior$B_0[ci, eq]
      prior_var  <- diag(post$prior$V_0)[ci] * post$S_n[eq, eq] / max(post$v_n - k_dom - 1, 1)
      prior_sd   <- sqrt(max(prior_var, 1e-12))
      prior_draws <- stats::rnorm(n_draws, mean = prior_mean, sd = prior_sd)

      # Posterior: N(B_n[ci, eq], V_n[ci, ci] * Sigma_post[eq, eq])
      post_mean <- post$B_n[ci, eq]
      Sigma_post <- post$S_n / max(post$v_n - k_dom - 1, 1)
      post_var  <- post$V_n[ci, ci] * Sigma_post[eq, eq]
      post_sd   <- sqrt(max(post_var, 1e-12))
      post_draws <- stats::rnorm(n_draws, mean = post_mean, sd = post_sd)

      label <- paste0("B[", ci, ",", eq, "]")

      df_list[[length(df_list) + 1]] <- data.frame(
        coef  = label,
        value = c(prior_draws, post_draws),
        type  = rep(c("Prior", "Posterior"), each = n_draws),
        stringsAsFactors = FALSE
      )
    }
  }

  df <- do.call(rbind, df_list)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = value, fill = type)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::facet_wrap(~ coef, scales = "free") +
    ggplot2::scale_fill_manual(values = c("Prior" = "grey70", "Posterior" = "steelblue")) +
    ggplot2::labs(
      title    = paste0("Prior vs Posterior: Unit '", unit_name, "'"),
      subtitle = paste0("Minnesota prior (lambda_1 = ", bayesian_gvar$lambda$lambda_1, ")"),
      x = "Coefficient value", y = "Density", fill = ""
    ) +
    ggplot2::theme_minimal()

  print(p)
  invisible(p)
}


# ───────────────────────────────────────────────────────────────────────────────
# 5.  Bayesian vs Frequentist Model Comparison
# ───────────────────────────────────────────────────────────────────────────────

#' Side-by-side comparison of Bayesian and frequentist GVAR results.
#'
#' Overlays Bayesian posterior credible intervals with frequentist bootstrap
#' confidence bands for the GIRF.
#'
#' @param bayesian_irf   Output from bayesian_girf()
#' @param frequentist_irf Output from bootstrap_girf()
#' @param response_vars   Variables to plot (NULL = all)
#' @return                A ggplot2 object (also printed)
bayesian_model_comparison <- function(bayesian_irf, frequentist_irf,
                                       response_vars = NULL) {

  var_names <- bayesian_irf$var_names
  horizon   <- nrow(bayesian_irf$point) - 1

  if (is.null(response_vars)) response_vars <- var_names

  df_list <- list()
  for (v in response_vars) {
    v_idx <- match(v, var_names)
    if (is.na(v_idx)) next

    # Bayesian
    df_list[[length(df_list) + 1]] <- data.frame(
      horizon  = 0:horizon,
      variable = v,
      method   = "Bayesian",
      median   = bayesian_irf$median[, v_idx],
      lower    = bayesian_irf$lower[, v_idx],
      upper    = bayesian_irf$upper[, v_idx],
      stringsAsFactors = FALSE
    )

    # Frequentist
    freq_med <- if (!is.null(frequentist_irf$median)) {
      frequentist_irf$median[, v_idx]
    } else {
      frequentist_irf$point[, v_idx]
    }

    n_freq <- min(horizon + 1, nrow(frequentist_irf$lower))
    df_list[[length(df_list) + 1]] <- data.frame(
      horizon  = 0:(n_freq - 1),
      variable = v,
      method   = "Frequentist",
      median   = freq_med[1:n_freq],
      lower    = frequentist_irf$lower[1:n_freq, v_idx],
      upper    = frequentist_irf$upper[1:n_freq, v_idx],
      stringsAsFactors = FALSE
    )
  }

  df <- do.call(rbind, df_list)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = horizon, color = method, fill = method)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.15) +
    ggplot2::geom_line(ggplot2::aes(y = median), linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
    ggplot2::facet_wrap(~ variable, scales = "free_y") +
    ggplot2::scale_color_manual(values = c("Bayesian" = "steelblue",
                                            "Frequentist" = "darkred")) +
    ggplot2::scale_fill_manual(values = c("Bayesian" = "steelblue",
                                           "Frequentist" = "darkred")) +
    ggplot2::labs(
      title = "Bayesian vs Frequentist GIRF Comparison",
      x = "Horizon", y = "Response",
      color = "Method", fill = "Method"
    ) +
    ggplot2::theme_minimal()

  print(p)
  invisible(p)
}


# ───────────────────────────────────────────────────────────────────────────────
# 6.  Bayesian In-Sample Diagnostics (Wrapper)
# ───────────────────────────────────────────────────────────────────────────────

#' In-sample diagnostics for the Bayesian GVAR posterior mean.
#'
#' Calls the frequentist insample_diagnostics() on the Bayesian point
#' estimate (posterior mean), which is a gvar_model object.
#'
#' @param bayesian_gvar  Output from bayesian_estimate_gvar()
#' @param lb_lags        Ljung-Box lag count (default 4)
#' @return               Output from insample_diagnostics()
bayesian_insample_diagnostics <- function(bayesian_gvar, lb_lags = 4) {
  print_banner("Bayesian In-Sample Diagnostics (Posterior Mean)")
  insample_diagnostics(bayesian_gvar$point_estimate, lb_lags = lb_lags)
}


# ───────────────────────────────────────────────────────────────────────────────
# 7.  Bayesian Out-of-Sample Forecast (Wrapper)
# ───────────────────────────────────────────────────────────────────────────────

#' Recursive out-of-sample forecast using the Bayesian posterior mean model.
#'
#' @param bayesian_gvar  Output from bayesian_estimate_gvar()
#' @param data_list      Named list of T x k_i matrices
#' @param star_list      Named list of T x k* matrices
#' @param W              Weight matrix
#' @param h              Forecast horizon
#' @param t0_frac        Fraction of data used for initial estimation
#' @return               Output from recursive_oos_forecast()
bayesian_oos_forecast <- function(bayesian_gvar, data_list, star_list, W,
                                   h = 1, t0_frac = 0.7) {
  print_banner("Bayesian Out-of-Sample Forecast (Posterior Mean)")
  recursive_oos_forecast(bayesian_gvar$point_estimate,
                          data_list = data_list,
                          star_list = star_list,
                          W = W, h = h, t0_frac = t0_frac)
}


# ───────────────────────────────────────────────────────────────────────────────
# 8.  Bayesian Out-of-Sample Evaluation (Wrapper)
# ───────────────────────────────────────────────────────────────────────────────

#' Evaluate OOS forecast performance for the Bayesian model.
#'
#' @param oos_result  Output from bayesian_oos_forecast()
#' @param t0_frac     Fraction used for initial estimation
#' @return            Output from oos_evaluation()
bayesian_oos_evaluation <- function(oos_result, t0_frac = 0.7) {
  print_banner("Bayesian Out-of-Sample Evaluation")
  oos_evaluation(oos_result, t0_frac = t0_frac)
}


# ───────────────────────────────────────────────────────────────────────────────
# 9.  Bayesian Residual Diagnostics Plot (Wrapper)
# ───────────────────────────────────────────────────────────────────────────────

#' Plot residual diagnostics (histogram, ACF, QQ) for the Bayesian point estimate.
#'
#' @param bayesian_gvar  Output from bayesian_estimate_gvar()
#' @param variables      Variables to plot (NULL = auto-select)
#' @param max_vars       Max number of variables per plot
#' @return               ggplot2 object
plot_bayesian_residual_diagnostics <- function(bayesian_gvar, variables = NULL,
                                                max_vars = 4) {
  plot_residual_diagnostics(bayesian_gvar$point_estimate,
                             variables = variables, max_vars = max_vars)
}


# ───────────────────────────────────────────────────────────────────────────────
# 10. Bayesian Fitted vs Actual Plot (Wrapper)
# ───────────────────────────────────────────────────────────────────────────────

#' Plot fitted vs actual values from the Bayesian posterior mean.
#'
#' @param bayesian_gvar  Output from bayesian_estimate_gvar()
#' @param variables      Variables to plot (NULL = auto-select)
#' @param max_vars       Max number of variables per plot
#' @return               ggplot2 object
plot_bayesian_fitted_vs_actual <- function(bayesian_gvar, variables = NULL,
                                            max_vars = 4) {
  plot_fitted_vs_actual(bayesian_gvar$point_estimate,
                         variables = variables, max_vars = max_vars)
}


# ───────────────────────────────────────────────────────────────────────────────
# 11. Bayesian Eigenvalue Plot (Wrapper)
# ───────────────────────────────────────────────────────────────────────────────

#' Plot eigenvalues (stability check) from the Bayesian posterior mean.
#'
#' @param bayesian_gvar  Output from bayesian_estimate_gvar()
#' @return               ggplot2 object
plot_bayesian_eigenvalues <- function(bayesian_gvar) {
  plot_eigenvalues(bayesian_gvar$point_estimate)
}


# ───────────────────────────────────────────────────────────────────────────────
# 12. Bayesian CUSUM Plot (Wrapper)
# ───────────────────────────────────────────────────────────────────────────────

#' Plot CUSUM statistics for the Bayesian posterior mean residuals.
#'
#' @param bayesian_gvar  Output from bayesian_estimate_gvar()
#' @param variables      Variables to plot (NULL = auto-select)
#' @param sig_level      Significance level for CUSUM bands
#' @return               ggplot2 object
plot_bayesian_cusum <- function(bayesian_gvar, variables = NULL,
                                 sig_level = 0.05) {
  plot_cusum(bayesian_gvar$point_estimate,
              variables = variables, sig_level = sig_level)
}


# ───────────────────────────────────────────────────────────────────────────────
# 13. Posterior Eigenvalue Distribution Plot (New)
# ───────────────────────────────────────────────────────────────────────────────

#' Plot the posterior distribution of the maximum eigenvalue modulus.
#'
#' Uses all stable posterior draws of the companion matrix to compute
#' the max eigenvalue modulus at each draw. Displays a histogram + density
#' with a vertical line at the unit circle threshold.
#'
#' @param bayesian_gvar  Output from bayesian_estimate_gvar()
#' @return               ggplot2 object
plot_bayesian_eigenvalue_distribution <- function(bayesian_gvar) {

  stable_idx <- which(bayesian_gvar$draws$is_stable)
  n_use <- length(stable_idx)

  if (n_use < 5) {
    message("[Bayesian] Fewer than 5 stable draws – skipping eigenvalue distribution plot.")
    return(invisible(NULL))
  }

  message(sprintf("  Computing max eigenvalue modulus across %d stable draws ...", n_use))

  max_eig <- numeric(n_use)
  for (i in seq_len(n_use)) {
    s <- stable_idx[i]
    eig_mod <- Mod(eigen(bayesian_gvar$draws$companion[, , s], only.values = TRUE)$values)
    max_eig[i] <- max(eig_mod)
  }

  # Point estimate max eigenvalue
  pe_eig <- max(Mod(eigen(bayesian_gvar$point_estimate$companion,
                            only.values = TRUE)$values))

  df <- data.frame(max_eigenvalue = max_eig)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = max_eigenvalue)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                             bins = 40, fill = "steelblue", alpha = 0.6) +
    ggplot2::geom_density(color = "darkblue", linewidth = 1) +
    ggplot2::geom_vline(xintercept = 1, color = "red", linetype = "dashed",
                         linewidth = 1) +
    ggplot2::geom_vline(xintercept = pe_eig, color = "darkgreen", linetype = "solid",
                         linewidth = 0.8) +
    ggplot2::annotate("text", x = 1.002, y = Inf, label = "Unit Circle",
                       color = "red", hjust = 0, vjust = 2, size = 3.5) +
    ggplot2::annotate("text", x = pe_eig, y = Inf,
                       label = sprintf("Point est. = %.4f", pe_eig),
                       color = "darkgreen", hjust = 0, vjust = 4, size = 3.5) +
    ggplot2::labs(
      title    = "Posterior Distribution of Maximum Eigenvalue Modulus",
      subtitle = sprintf("%d stable draws | %.1f%% < 1.0",
                          n_use, 100 * mean(max_eig < 1)),
      x = "Max Eigenvalue Modulus", y = "Density"
    ) +
    ggplot2::theme_minimal()

  print(p)
  invisible(p)
}


# ───────────────────────────────────────────────────────────────────────────────
# 14. Master Bayesian Diagnostics Runner
# ───────────────────────────────────────────────────────────────────────────────

#' Run comprehensive diagnostics for the Bayesian GVAR.
#'
#' Combines Bayesian-specific diagnostics (DIC, marginal likelihood, PPC)
#' with the full frequentist diagnostic suite applied to the posterior mean.
#'
#' @param bayesian_gvar  Output from bayesian_estimate_gvar()
#' @param gvar_data      Output from prepare_gvar_dataset()
#' @param data_list      Named list of T x k_i matrices
#' @param star_list      Named list of T x k* matrices
#' @param W              Weight matrix
#' @param h              OOS forecast horizon
#' @param t0_frac        OOS initial estimation fraction
#' @param plot           Logical; generate diagnostic plots
#' @param save_pdf       Logical; save plots to PDF
#' @param pdf_file       PDF file name
#' @return               List with all diagnostic results
run_bayesian_all_diagnostics <- function(bayesian_gvar, gvar_data,
                                          data_list, star_list, W,
                                          h = 1, t0_frac = 0.7,
                                          plot = TRUE,
                                          save_pdf = FALSE,
                                          pdf_file = "bayesian_diagnostics.pdf") {

  print_banner("Comprehensive Bayesian GVAR Diagnostics")

  results <- list()

  if (save_pdf) {
    grDevices::pdf(pdf_file, width = 12, height = 8)
    on.exit(grDevices::dev.off(), add = TRUE)
    message(sprintf("  Saving all plots to '%s' ...", pdf_file))
  }

  # ---- 1. Bayesian-specific: DIC, marginal likelihood, PPC ----
  message("\n  [1/5] Bayesian-specific diagnostics (DIC, ML, PPC) ...")
  results$bayesian_specific <- tryCatch(
    bayesian_diagnostics(bayesian_gvar, gvar_data),
    error = function(e) {
      message("  [Warning] Bayesian-specific diagnostics failed: ", e$message)
      NULL
    }
  )

  # ---- 2. In-sample diagnostics (R^2, Ljung-Box, Jarque-Bera) ----
  message("\n  [2/5] In-sample diagnostics (posterior mean) ...")
  results$insample <- tryCatch(
    bayesian_insample_diagnostics(bayesian_gvar),
    error = function(e) {
      message("  [Warning] In-sample diagnostics failed: ", e$message)
      NULL
    }
  )

  # ---- 3. Out-of-sample forecasting ----
  message("\n  [3/5] Out-of-sample recursive forecast ...")
  results$oos_forecast <- tryCatch(
    bayesian_oos_forecast(bayesian_gvar, data_list, star_list, W,
                           h = h, t0_frac = t0_frac),
    error = function(e) {
      message("  [Warning] OOS forecast failed: ", e$message)
      NULL
    }
  )

  # ---- 4. OOS evaluation (RMSE, MAE, Theil's U) ----
  if (!is.null(results$oos_forecast)) {
    message("\n  [4/5] OOS evaluation metrics ...")
    results$oos_eval <- tryCatch(
      bayesian_oos_evaluation(results$oos_forecast, t0_frac = t0_frac),
      error = function(e) {
        message("  [Warning] OOS evaluation failed: ", e$message)
        NULL
      }
    )
  }

  # ---- 5. Plots ----
  if (plot) {
    message("\n  [5/5] Generating diagnostic plots ...")

    # Residual diagnostics
    tryCatch(
      plot_bayesian_residual_diagnostics(bayesian_gvar),
      error = function(e) message("  [Warning] Residual plot failed: ", e$message)
    )

    # Fitted vs actual
    tryCatch(
      plot_bayesian_fitted_vs_actual(bayesian_gvar),
      error = function(e) message("  [Warning] Fitted vs actual plot failed: ", e$message)
    )

    # Eigenvalues (point estimate)
    tryCatch(
      plot_bayesian_eigenvalues(bayesian_gvar),
      error = function(e) message("  [Warning] Eigenvalue plot failed: ", e$message)
    )

    # Posterior eigenvalue distribution
    tryCatch(
      plot_bayesian_eigenvalue_distribution(bayesian_gvar),
      error = function(e) message("  [Warning] Eigenvalue distribution plot failed: ", e$message)
    )

    # CUSUM
    tryCatch(
      plot_bayesian_cusum(bayesian_gvar),
      error = function(e) message("  [Warning] CUSUM plot failed: ", e$message)
    )
  }

  message("\n[Bayesian Diagnostics] Complete.")
  return(results)
}
