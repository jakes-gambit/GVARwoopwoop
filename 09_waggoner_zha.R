###############################################################################
#  09_waggoner_zha.R  –  Conditional Forecasting à la Waggoner & Zha (1999)
#
#  This module implements the conditional forecasting methodology of
#  Waggoner & Zha (1999, "Conditional Forecasts in Dynamic Multivariate
#  Models", Review of Economics and Statistics).
#
#  CONCEPTUAL OVERVIEW
#  -------------------
#  The Waggoner-Zha (WZ) approach differs fundamentally from the Kalman
#  filter method (08_kalman_filter.R):
#
#  • The Kalman filter treats conditions as "observations" of the state
#    and updates beliefs via Bayes' rule in state space.
#
#  • The WZ approach works in the SHOCK SPACE:
#      1. Express the VAR in its Moving-Average (MA) form so that future
#         values are explicit functions of future structural shocks.
#      2. Conditions (imposed future paths) become LINEAR CONSTRAINTS on
#         the sequence of future shocks.
#      3. Find the MINIMUM-NORM (least-energy) shock sequence satisfying
#         the constraints → the "most likely" conditional forecast.
#      4. Use a GIBBS SAMPLER to draw from the conditional distribution
#         of shocks, yielding proper posterior uncertainty bands.
#
#  WHY USE WZ INSTEAD OF / IN ADDITION TO THE KALMAN FILTER?
#  ----------------------------------------------------------
#  • WZ gives forecasts that are explicitly "closest" to the unconditional
#    forecast (in the Mahalanobis-distance sense), which is a natural
#    economic interpretation.
#  • The Gibbs sampler produces exact finite-sample posterior bands
#    (under Gaussianity), not the approximate Gaussian bands of the KF.
#  • WZ naturally handles structural identification (Cholesky, sign
#    restrictions, etc.) since it operates on structural shocks.
#  • Soft constraints can be implemented by adding slack to the constraint
#    system rather than inflating measurement noise.
#
#  MATHEMATICAL FORMULATION
#  ------------------------
#  The reduced-form global VAR:
#
#      x_t = F_0 + F_1 x_{t-1} + ... + F_p x_{t-p} + e_t
#
#  where e_t ~ N(0, Sigma).  Let P be the lower-Cholesky factor of Sigma
#  (Sigma = P P'), so e_t = P eps_t  with  eps_t ~ N(0, I_k).
#
#  The h-step-ahead forecast can be written in MA form:
#
#      x_{T+h} = mu_h  +  sum_{j=0}^{h-1}  Psi_j  P  eps_{T+h-j}
#
#  where mu_h is the unconditional (no-shock) forecast at horizon h and
#  Psi_j are the MA coefficient matrices (Psi_0 = I_k).
#
#  Stack the forecasts for h = 1,...,H:
#
#      x_tilde = mu_tilde  +  M  eps_tilde
#
#  where x_tilde = [x_{T+1}; ...; x_{T+H}]  is (kH x 1),
#        eps_tilde = [eps_{T+1}; ...; eps_{T+H}]  is (kH x 1),
#        M  is a (kH x kH) lower-block-triangular matrix built from Psi_j P.
#
#  Conditions select certain elements of x_tilde:
#
#      R  x_tilde  =  r        (m constraints)
#
#  where R is an (m x kH) selection matrix and r holds the imposed values.
#  Substituting the MA form:
#
#      R M eps_tilde  =  r  -  R mu_tilde   equiv   Phi eps_tilde  =  delta
#
#  The minimum-norm solution:
#
#      eps_star  =  Phi' (Phi Phi')^{-1}  delta
#
#  and the conditional forecast is:
#
#      x_tilde_cond  =  mu_tilde  +  M  eps_star
#
#  GIBBS SAMPLER FOR UNCERTAINTY BANDS
#  ------------------------------------
#  Any eps_tilde satisfying the constraints can be written as:
#
#      eps_tilde  =  eps_star  +  (I - Phi' (Phi Phi')^{-1} Phi)  nu
#
#  where nu ~ N(0, I_{kH}).  The sampler draws nu, projects out the
#  constrained directions, and adds the particular solution eps_star.
#
#  Contents
#  --------
#  1.  compute_ma_coefficients()      - Psi_0, Psi_1, ..., Psi_{H-1}
#  2.  build_M_matrix()               - the (kH x kH) MA shock-mapping matrix
#  3.  unconditional_forecast_path()  - mu_tilde (the no-shock baseline)
#  4.  build_constraint_system()      - R, r -> Phi, delta
#  5.  wz_point_forecast()            - analytical minimum-norm solution
#  6.  wz_gibbs_sampler()             - Gibbs draws for uncertainty bands
#  7.  conditional_forecast_wz()      - master function
#  8.  plot_conditional_forecast_wz() - visualisation
#  9.  compare_kf_wz()               - side-by-side comparison with KF
###############################################################################


# =============================================================================
# 1.  MA Coefficient Matrices
# =============================================================================

#' Compute the Moving-Average coefficient matrices Psi_0, Psi_1, ..., Psi_{H-1}
#' from the companion-form VAR.
#'
#' The MA representation of x_{T+h} is:
#'    x_{T+h}  =  deterministic  +  sum_{j=0}^{h-1}  Psi_j  e_{T+h-j}
#'
#' where Psi_0 = I_k and Psi_j = J * C^j * J'  with C being the companion
#' matrix and J = [I_k  0  ...  0] the selector.
#'
#' @param gvar_model  Fitted GVAR model
#' @param max_h       Maximum horizon
#' @return            List of (max_h) matrices, each k x k.
#'                    Element [[j+1]] = Psi_j  (j = 0,...,max_h-1)
compute_ma_coefficients <- function(gvar_model, max_h) {

  k_total <- gvar_model$k_total
  comp    <- gvar_model$companion
  kp      <- nrow(comp)

  # Selector matrix J:  (k_total x kp)
  # Picks out the first k_total rows of the companion state
  J <- cbind(diag(k_total), matrix(0, k_total, kp - k_total))

  Psi <- vector("list", max_h)

  # Psi_0 = I_k
  Psi[[1]] <- diag(k_total)

  # Psi_j = J * C^j * J'    for j = 1, ..., H-1
  C_power <- comp   # C^1
  for (j in seq_len(max_h - 1)) {
    Psi[[j + 1]] <- J %*% C_power %*% t(J)
    C_power <- C_power %*% comp   # C^{j+1}
  }

  return(Psi)
}


# =============================================================================
# 2.  Build the M Matrix  (Shock-to-Forecast Mapping)
# =============================================================================

#' Construct the (kH x kH) lower-block-triangular matrix M that maps the
#' stacked structural shocks eps_tilde to the stacked forecast deviations.
#'
#' M has block structure:
#'
#'   | Psi_0 P        0             0          ...  0          |
#'   | Psi_1 P        Psi_0 P       0          ...  0          |
#'   | Psi_2 P        Psi_1 P       Psi_0 P    ...  0          |
#'   |   :              :             :              :          |
#'   | Psi_{H-1} P    Psi_{H-2} P   ...             Psi_0 P    |
#'
#' @param Psi   List of MA coefficient matrices (from compute_ma_coefficients)
#' @param P     k x k lower-Cholesky factor of Sigma  (Sigma = P P')
#' @param max_h Forecast horizon
#' @return      (k*max_h) x (k*max_h) matrix
build_M_matrix <- function(Psi, P, max_h) {

  k <- nrow(P)
  kH <- k * max_h
  M <- matrix(0, nrow = kH, ncol = kH)

  for (row_h in seq_len(max_h)) {
    for (col_h in seq_len(row_h)) {
      # Block (row_h, col_h) corresponds to Psi_{row_h - col_h} * P
      j <- row_h - col_h   # MA index (0-based)
      rows <- ((row_h - 1) * k + 1):(row_h * k)
      cols <- ((col_h - 1) * k + 1):(col_h * k)
      M[rows, cols] <- Psi[[j + 1]] %*% P
    }
  }

  return(M)
}


# =============================================================================
# 3.  Unconditional Forecast Path  (mu_tilde)
# =============================================================================

#' Compute the unconditional (baseline, no-shock) forecast path mu_tilde.
#'
#' This is obtained by iterating the companion form forward with zero shocks:
#'    alpha_{T+h} = T * alpha_{T+h-1} + c      (no shock term)
#'
#' @param gvar_model  Fitted GVAR model
#' @param data_list   Named list of raw data matrices (for initial state)
#' @param max_h       Forecast horizon
#' @return            (k*max_h) x 1 vector  mu_tilde = [mu_1; mu_2; ...; mu_H]
unconditional_forecast_path <- function(gvar_model, data_list, max_h) {

  k_total  <- gvar_model$k_total
  p_global <- gvar_model$p_global
  kp       <- k_total * p_global
  comp     <- gvar_model$companion
  F0       <- gvar_model$F0

  # Companion intercept
  cc <- c(F0, rep(0, kp - k_total))

  # Build initial companion state from the last p observations.
  # Use gvar_model$var_names to select only the model variables in the right
  # order (non-dominant units may have fewer columns than data_list[[unit]]).
  TT_data <- nrow(data_list[[names(data_list)[1]]])
  vnames  <- gvar_model$var_names

  state <- c()
  for (l in 0:(p_global - 1)) {
    t_idx <- TT_data - l
    x_t   <- numeric(k_total)
    for (v in seq_len(k_total)) {
      vname  <- vnames[v]
      dot    <- regexpr("\\.", vname)[1]
      u      <- substr(vname, 1, dot - 1)
      varstr <- substr(vname, dot + 1, nchar(vname))
      if (u %in% names(data_list) && varstr %in% colnames(data_list[[u]])) {
        x_t[v] <- data_list[[u]][t_idx, varstr]
      }
    }
    state <- c(state, x_t)
  }

  # Iterate forward, collecting the first k_total elements at each step
  mu_tilde <- numeric(k_total * max_h)

  for (h in seq_len(max_h)) {
    state <- as.numeric(comp %*% state + cc)
    idx <- ((h - 1) * k_total + 1):(h * k_total)
    mu_tilde[idx] <- state[1:k_total]
  }

  return(mu_tilde)
}


# =============================================================================
# 4.  Build the Constraint System  (Phi, delta)
# =============================================================================

#' Translate user-supplied conditions into the constraint matrices Phi and delta.
#'
#' The conditions select specific elements of the stacked forecast x_tilde:
#'
#'    R  x_tilde  =  r
#'
#' Substituting x_tilde = mu_tilde + M eps_tilde :
#'
#'    R M eps_tilde  =  r - R mu_tilde    <=>    Phi eps_tilde = delta
#'
#' @param conditions   Data frame with columns: variable, horizon, value.
#'                     Optionally: type ("hard" or "soft"), tolerance (for soft).
#' @param gvar_model   Fitted GVAR model
#' @param M            (kH x kH) MA mapping matrix
#' @param mu_tilde     (kH x 1)  unconditional forecast path
#' @param max_h        Forecast horizon
#' @return  List with Phi, delta, R_mat, r_vec, and soft-constraint metadata
build_constraint_system <- function(conditions, gvar_model, M, mu_tilde, max_h) {

  k_total <- gvar_model$k_total
  kH      <- k_total * max_h
  vnames  <- gvar_model$var_names

  conditions <- as.data.frame(conditions, stringsAsFactors = FALSE)
  m <- nrow(conditions)

  # Default: all hard constraints
  if (!"type" %in% names(conditions)) conditions$type <- "hard"
  if (!"tolerance" %in% names(conditions)) conditions$tolerance <- 0

  # Build R (selection matrix) and r (values vector)
  R_mat <- matrix(0, nrow = m, ncol = kH)
  r_vec <- numeric(m)

  for (i in seq_len(m)) {
    var_name <- conditions$variable[i]
    h_cond   <- conditions$horizon[i]
    val      <- conditions$value[i]

    var_idx <- which(vnames == var_name)
    if (length(var_idx) == 0)
      stop("[WZ] Variable '", var_name, "' not found in the GVAR model.")
    if (h_cond < 1 || h_cond > max_h)
      stop("[WZ] Horizon ", h_cond, " is outside the range [1, ", max_h, "].")

    # Position in the stacked vector x_tilde
    pos <- (h_cond - 1) * k_total + var_idx
    R_mat[i, pos] <- 1
    r_vec[i] <- val
  }

  # Phi = R M
  Phi <- R_mat %*% M

  # delta = r - R mu_tilde
  delta <- r_vec - as.numeric(R_mat %*% mu_tilde)

  # Identify soft constraints
  soft_idx <- which(conditions$type == "soft")
  n_soft <- length(soft_idx)
  n_hard <- m - n_soft
  soft_tol <- conditions$tolerance[soft_idx]

  return(list(
    Phi      = Phi,
    delta    = delta,
    R_mat    = R_mat,
    r_vec    = r_vec,
    n_hard   = n_hard,
    n_soft   = n_soft,
    soft_idx = soft_idx,
    soft_tol = soft_tol
  ))
}


# =============================================================================
# 5.  Analytical Point Forecast  (Minimum-Norm Solution)
# =============================================================================

#' Compute the WZ minimum-norm conditional forecast.
#'
#' The particular solution is:
#'
#'    eps_star  =  Phi' (Phi Phi')^{-1}  delta
#'
#' and the conditional forecast path:
#'
#'    x_tilde_cond  =  mu_tilde  +  M eps_star
#'
#' This gives the forecast that satisfies ALL hard constraints while
#' requiring the smallest possible (Euclidean-norm) structural shocks --
#' i.e. the forecast "closest" to the unconditional baseline.
#'
#' @param cs          Constraint system (from build_constraint_system)
#' @param M           MA mapping matrix
#' @param mu_tilde    Unconditional forecast path
#' @param k_total     Number of variables
#' @param max_h       Forecast horizon
#' @return List with eps_star, x_cond (stacked), x_cond_mat (matrix form)
wz_point_forecast <- function(cs, M, mu_tilde, k_total, max_h) {

  Phi   <- cs$Phi
  delta <- cs$delta

  # Phi' (Phi Phi')^{-1} delta
  PhiPhiT <- Phi %*% t(Phi)

  # Use generalised inverse for robustness
  PhiPhiT_inv <- tryCatch(
    solve(PhiPhiT),
    error = function(e) {
      message("[WZ] Constraint matrix near-singular -- using generalised inverse.")
      MASS::ginv(PhiPhiT)
    }
  )

  eps_star <- as.numeric(t(Phi) %*% PhiPhiT_inv %*% delta)

  # Conditional forecast
  x_cond <- mu_tilde + as.numeric(M %*% eps_star)

  # Reshape to matrix
  x_cond_mat <- matrix(x_cond, nrow = max_h, ncol = k_total, byrow = TRUE)

  return(list(
    eps_star   = eps_star,
    x_cond     = x_cond,
    x_cond_mat = x_cond_mat
  ))
}


# =============================================================================
# 6.  Gibbs Sampler for Uncertainty Bands
# =============================================================================

#' Gibbs sampler for the WZ conditional forecast distribution.
#'
#' Any shock sequence satisfying the constraints can be written as:
#'
#'    eps_tilde  =  eps_star  +  N_Phi  nu
#'
#' where N_Phi is the null-space projector  N_Phi = I - Phi' (Phi Phi')^{-1} Phi
#' and nu ~ N(0, I_{kH}).
#'
#' The sampler:
#'   1. Draw  nu ~ N(0, I_{kH})
#'   2. Project:  eps^{(s)} = eps_star + N_Phi nu
#'   3. Compute:  x^{(s)} = mu_tilde + M eps^{(s)}
#'   4. Repeat and collect quantiles.
#'
#' For SOFT constraints, the sampler additionally includes a
#' Metropolis-Hastings accept/reject step: draws that violate the soft
#' constraint beyond its tolerance are accepted with probability
#' proportional to exp(-0.5 * violation^2 / tolerance^2).
#'
#' @param cs          Constraint system
#' @param M           MA mapping matrix
#' @param mu_tilde    Unconditional forecast path
#' @param k_total     Number of variables
#' @param max_h       Forecast horizon
#' @param n_draws     Number of Gibbs draws (default 2000)
#' @param burn_in     Number of burn-in draws to discard (default 500)
#' @param ci_level    Confidence level (default 0.90)
#' @param seed        Random seed
#' @return List with draws array, mean, median, lower, upper, point
wz_gibbs_sampler <- function(cs, M, mu_tilde, k_total, max_h,
                              n_draws = 2000, burn_in = 500,
                              ci_level = 0.90, seed = 42) {

  set.seed(seed)

  kH    <- k_total * max_h
  Phi   <- cs$Phi
  delta <- cs$delta

  # ---- Compute the null-space projector N_Phi ----
  # N_Phi = I_{kH} - Phi' (Phi Phi')^{-1} Phi
  PhiPhiT <- Phi %*% t(Phi)
  PhiPhiT_inv <- tryCatch(solve(PhiPhiT), error = function(e) MASS::ginv(PhiPhiT))

  N_Phi <- diag(kH) - t(Phi) %*% PhiPhiT_inv %*% Phi

  # ---- Particular solution eps_star ----
  eps_star <- as.numeric(t(Phi) %*% PhiPhiT_inv %*% delta)

  # ---- Pre-compute for speed ----
  # x_tilde = mu_tilde + M (eps_star + N_Phi nu)  =  (mu_tilde + M eps_star) + M N_Phi nu
  x_base <- mu_tilde + as.numeric(M %*% eps_star)   # conditional point forecast
  MN <- M %*% N_Phi                                  # pre-multiply once

  # ---- Storage ----
  total_draws <- n_draws + burn_in
  draws_array <- array(NA, dim = c(n_draws, max_h, k_total))

  # Soft-constraint bookkeeping
  has_soft <- cs$n_soft > 0
  if (has_soft) {
    soft_rows <- cs$soft_idx
    soft_tol  <- cs$soft_tol
    R_soft <- cs$R_mat[soft_rows, , drop = FALSE]
    r_soft <- cs$r_vec[soft_rows]
  }

  accepted <- 0
  draw_idx <- 0

  message(sprintf("  [WZ Gibbs] Running %d draws (%d burn-in) ...",
                  n_draws, burn_in))

  # Over-iterate to account for rejections in soft-constraint case
  max_iter <- total_draws * ifelse(has_soft, 5, 1)

  for (s in seq_len(max_iter)) {

    # Step 1: draw standard normal
    nu <- rnorm(kH)

    # Step 2: project into the feasible subspace
    x_draw_vec <- x_base + as.numeric(MN %*% nu)

    # Step 3: soft-constraint accept/reject (if applicable)
    accept <- TRUE
    if (has_soft) {
      violations <- as.numeric(R_soft %*% x_draw_vec) - r_soft
      log_prob <- -0.5 * sum((violations / soft_tol)^2)
      if (log(runif(1)) > log_prob) accept <- FALSE
    }

    if (!accept) next

    accepted <- accepted + 1
    if (accepted <= burn_in) next   # discard burn-in

    draw_idx <- draw_idx + 1
    if (draw_idx > n_draws) break

    # Reshape to max_h x k_total and store
    draws_array[draw_idx, , ] <- matrix(x_draw_vec, nrow = max_h,
                                         ncol = k_total, byrow = TRUE)
  }

  if (draw_idx < n_draws) {
    warning("[WZ Gibbs] Only ", draw_idx, " / ", n_draws,
            " draws accepted. Consider relaxing soft constraints.")
    draws_array <- draws_array[seq_len(draw_idx), , , drop = FALSE]
  }

  # ---- Summarise ----
  alpha <- 1 - ci_level
  mean_fc   <- apply(draws_array, c(2, 3), mean, na.rm = TRUE)
  median_fc <- apply(draws_array, c(2, 3), median, na.rm = TRUE)
  lower_fc  <- apply(draws_array, c(2, 3), quantile, probs = alpha / 2, na.rm = TRUE)
  upper_fc  <- apply(draws_array, c(2, 3), quantile, probs = 1 - alpha / 2, na.rm = TRUE)

  # Point forecast (minimum-norm)
  point_fc <- matrix(x_base, nrow = max_h, ncol = k_total, byrow = TRUE)

  message(sprintf("  [WZ Gibbs] Done. %d draws retained.", draw_idx))

  return(list(
    draws    = draws_array,
    mean     = mean_fc,
    median   = median_fc,
    lower    = lower_fc,
    upper    = upper_fc,
    point    = point_fc,
    ci_level = ci_level,
    n_draws  = draw_idx
  ))
}


# =============================================================================
# 7.  Master Function: Conditional Forecast (Waggoner-Zha)
# =============================================================================

#' Produce conditional forecasts using the Waggoner & Zha (1999) methodology.
#'
#' @param gvar_model   Fitted GVAR model (from estimate_gvar)
#' @param conditions   Data frame with columns:
#'   \describe{
#'     \item{variable}{Character; global name e.g. "US.gdp"}
#'     \item{horizon}{Integer; 1 = one step ahead}
#'     \item{value}{Numeric; the imposed value}
#'     \item{type}{Character; "hard" (default) or "soft"}
#'     \item{tolerance}{Numeric; s.d. of allowed deviation for soft constraints}
#'   }
#' @param max_h        Maximum forecast horizon (default: max of conditions$horizon)
#' @param data_list    Named list of raw data matrices (for initial state)
#' @param n_draws      Gibbs sampler draws (default 2000)
#' @param burn_in      Burn-in draws (default 500)
#' @param ci_level     Confidence level for bands (default 0.90)
#' @param seed         Random seed
#' @return A list (class "wz_conditional_forecast") with forecasts,
#'         conditions, Gibbs output, MA structures, etc.
conditional_forecast_wz <- function(gvar_model, conditions, max_h = NULL,
                                     data_list = NULL,
                                     n_draws = 2000, burn_in = 500,
                                     ci_level = 0.90, seed = 42) {

  print_banner("Conditional Forecasting (Waggoner & Zha 1999)")

  conditions <- as.data.frame(conditions, stringsAsFactors = FALSE)

  if (is.null(max_h)) max_h <- max(conditions$horizon)
  if (is.null(data_list)) stop("[WZ] data_list is required for the initial state.")

  k_total  <- gvar_model$k_total
  vnames   <- gvar_model$var_names

  # -- Step 1: MA coefficients --
  message("  Step 1/5: Computing MA coefficients ...")
  Psi <- compute_ma_coefficients(gvar_model, max_h)

  # -- Step 2: Cholesky factor of Sigma --
  message("  Step 2/5: Cholesky decomposition of Sigma ...")
  Sigma <- gvar_model$Sigma

  # Ensure Sigma is positive-definite (add small ridge if needed)
  eig_min <- min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values)
  if (eig_min <= 0) {
    ridge <- abs(eig_min) + 1e-6
    message(sprintf("    [!] Sigma not PD (min eigenvalue = %.2e). Adding ridge %.2e.",
                    eig_min, ridge))
    Sigma <- Sigma + ridge * diag(k_total)
  }
  P <- t(chol(Sigma))   # lower-triangular Cholesky factor

  # -- Step 3: Build M matrix --
  message("  Step 3/5: Building shock-mapping matrix M ...")
  M <- build_M_matrix(Psi, P, max_h)

  # -- Step 4: Unconditional forecast path --
  message("  Step 4/5: Computing unconditional forecast baseline ...")
  mu_tilde <- unconditional_forecast_path(gvar_model, data_list, max_h)

  # -- Step 5: Constraint system --
  message("  Step 5/5: Building constraint system and solving ...")
  cs <- build_constraint_system(conditions, gvar_model, M, mu_tilde, max_h)

  message(sprintf("    %d hard constraints, %d soft constraints.",
                  cs$n_hard, cs$n_soft))

  # -- Point forecast (analytical minimum-norm) --
  pf <- wz_point_forecast(cs, M, mu_tilde, k_total, max_h)

  # -- Gibbs sampler --
  gibbs <- wz_gibbs_sampler(
    cs, M, mu_tilde, k_total, max_h,
    n_draws = n_draws, burn_in = burn_in,
    ci_level = ci_level, seed = seed
  )

  # -- Assemble output --
  uncond_mat <- matrix(mu_tilde, nrow = max_h, ncol = k_total, byrow = TRUE)
  colnames(uncond_mat)       <- vnames
  colnames(pf$x_cond_mat)   <- vnames
  colnames(gibbs$mean)       <- vnames
  colnames(gibbs$median)     <- vnames
  colnames(gibbs$lower)      <- vnames
  colnames(gibbs$upper)      <- vnames

  rownames(uncond_mat)       <- paste0("h=", 1:max_h)
  rownames(pf$x_cond_mat)   <- paste0("h=", 1:max_h)
  rownames(gibbs$mean)       <- paste0("h=", 1:max_h)

  result <- list(
    forecast_point = pf$x_cond_mat,
    forecast_mean  = gibbs$mean,
    forecast_lower = gibbs$lower,
    forecast_upper = gibbs$upper,
    unconditional  = uncond_mat,
    conditions     = conditions,
    var_names      = vnames,
    ci_level       = ci_level,
    max_h          = max_h,
    gibbs          = gibbs,
    eps_star       = pf$eps_star,
    M              = M,
    Psi            = Psi,
    P_chol         = P
  )

  class(result) <- "wz_conditional_forecast"

  # Print summary
  cat("\n  WZ Conditional Forecast Summary\n")
  cat("  --------------------------------\n")
  cat(sprintf("  Horizon:          %d steps\n", max_h))
  cat(sprintf("  Variables:        %d\n", k_total))
  cat(sprintf("  Hard constraints: %d\n", cs$n_hard))
  cat(sprintf("  Soft constraints: %d\n", cs$n_soft))
  cat(sprintf("  Gibbs draws:      %d (burn-in: %d)\n", gibbs$n_draws, burn_in))
  cat(sprintf("  Confidence level: %.0f%%\n", ci_level * 100))
  cat(sprintf("  Shock norm:       %.4f\n", sqrt(sum(pf$eps_star^2))))
  cat("\n")

  cat("  Conditional forecast (minimum-norm point estimate):\n")
  print(round(pf$x_cond_mat, 4))
  cat("\n")

  return(result)
}


# =============================================================================
# 8.  Plotting
# =============================================================================

#' Visualise Waggoner-Zha conditional forecasts.
#'
#' Shows: unconditional baseline (dashed), conditional point forecast (solid),
#' Gibbs uncertainty bands (shaded), and condition points (red diamonds).
#'
#' @param wz_result    Output from conditional_forecast_wz()
#' @param plot_vars    Character or integer vector of variables to plot
#' @param show_uncond  Logical; if TRUE, overlay the unconditional forecast
#' @return  ggplot object
plot_conditional_forecast_wz <- function(wz_result, plot_vars = NULL,
                                          show_uncond = TRUE) {

  fp     <- wz_result$forecast_point
  fl     <- wz_result$forecast_lower
  fu     <- wz_result$forecast_upper
  uncond <- wz_result$unconditional
  conds  <- wz_result$conditions
  vnames <- wz_result$var_names
  max_h  <- wz_result$max_h

  if (is.null(plot_vars)) {
    plot_vars <- seq_len(min(ncol(fp), 9))
  }
  if (is.character(plot_vars)) {
    plot_vars <- match(plot_vars, vnames)
  }

  # Build long-format data
  plot_df <- data.frame()
  for (j in plot_vars) {
    vname <- vnames[j]

    df_j <- data.frame(
      horizon     = 1:max_h,
      point       = fp[, j],
      lower       = fl[, j],
      upper       = fu[, j],
      uncond      = uncond[, j],
      variable    = vname,
      stringsAsFactors = FALSE
    )
    plot_df <- rbind(plot_df, df_j)
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = horizon)) +

    # Uncertainty band (Gibbs)
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         fill = "#2E86AB", alpha = 0.20) +

    # Conditional forecast (point, solid)
    ggplot2::geom_line(ggplot2::aes(y = point), colour = "#2E86AB",
                       linewidth = 1.0) +

    # Zero line
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50")

  # Unconditional baseline (dashed)
  if (show_uncond) {
    p <- p + ggplot2::geom_line(ggplot2::aes(y = uncond),
                                colour = "grey40", linetype = "dashed",
                                linewidth = 0.6)
  }

  # Condition points (red diamonds)
  cond_points <- data.frame()
  for (j in plot_vars) {
    vname <- vnames[j]
    cond_j <- conds[conds$variable == vname, ]
    if (nrow(cond_j) > 0) {
      cond_points <- rbind(cond_points, data.frame(
        horizon  = cond_j$horizon,
        value    = cond_j$value,
        variable = vname,
        stringsAsFactors = FALSE
      ))
    }
  }

  if (nrow(cond_points) > 0) {
    p <- p + ggplot2::geom_point(
      data = cond_points,
      ggplot2::aes(x = horizon, y = value),
      colour = "#E83151", size = 3.5, shape = 18
    )
  }

  p <- p +
    ggplot2::facet_wrap(~ variable, scales = "free_y") +
    ggplot2::labs(
      title    = "Conditional Forecast (Waggoner & Zha 1999)",
      subtitle = paste0(
        wz_result$ci_level * 100, "% Gibbs bands  |  ",
        "Dashed = unconditional  |  ",
        "Diamond = imposed conditions"
      ),
      x = "Forecast Horizon", y = "Value"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold")
    )

  print(p)
  return(invisible(p))
}


# =============================================================================
# 9.  Compare Kalman Filter vs. Waggoner-Zha Forecasts
# =============================================================================

#' Side-by-side comparison of conditional forecasts from the Kalman filter
#' (08_kalman_filter.R) and the Waggoner-Zha approach.
#'
#' @param kf_result    Output from conditional_forecast()  (Kalman filter)
#' @param wz_result    Output from conditional_forecast_wz()
#' @param plot_vars    Variables to compare
#' @return  ggplot object showing both methods
compare_kf_wz <- function(kf_result, wz_result, plot_vars = NULL) {

  vnames <- wz_result$var_names
  max_h  <- min(kf_result$max_h, wz_result$max_h)

  if (is.null(plot_vars)) {
    plot_vars <- seq_len(min(length(vnames), 6))
  }
  if (is.character(plot_vars)) {
    plot_vars <- match(plot_vars, vnames)
  }

  plot_df <- data.frame()

  for (j in plot_vars) {
    vname <- vnames[j]

    # KF results
    df_kf <- data.frame(
      horizon  = 1:max_h,
      mean     = kf_result$forecast_mean[1:max_h, j],
      lower    = kf_result$forecast_lower[1:max_h, j],
      upper    = kf_result$forecast_upper[1:max_h, j],
      method   = "Kalman Filter",
      variable = vname,
      stringsAsFactors = FALSE
    )

    # WZ results
    df_wz <- data.frame(
      horizon  = 1:max_h,
      mean     = wz_result$forecast_point[1:max_h, j],
      lower    = wz_result$forecast_lower[1:max_h, j],
      upper    = wz_result$forecast_upper[1:max_h, j],
      method   = "Waggoner-Zha",
      variable = vname,
      stringsAsFactors = FALSE
    )

    plot_df <- rbind(plot_df, df_kf, df_wz)
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = horizon, fill = method,
                                              colour = method)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         alpha = 0.15, colour = NA) +
    ggplot2::geom_line(ggplot2::aes(y = mean), linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
    ggplot2::facet_wrap(~ variable, scales = "free_y") +
    ggplot2::scale_colour_manual(values = c("Kalman Filter" = "darkorange",
                                             "Waggoner-Zha"  = "#2E86AB")) +
    ggplot2::scale_fill_manual(values = c("Kalman Filter" = "darkorange",
                                           "Waggoner-Zha"  = "#2E86AB")) +
    ggplot2::labs(
      title    = "Conditional Forecast Comparison: Kalman Filter vs. Waggoner-Zha",
      subtitle = paste0(wz_result$ci_level * 100, "% confidence bands"),
      x = "Forecast Horizon", y = "Value",
      colour = "Method", fill = "Method"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  print(p)
  return(invisible(p))
}
