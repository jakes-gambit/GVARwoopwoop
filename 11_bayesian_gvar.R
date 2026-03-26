###############################################################################
#  11_bayesian_gvar.R  –  Bayesian GVAR with Minnesota / Litterman Prior
#
#  Implements Bayesian estimation of the GVAR using the conjugate
#  Normal-Inverse-Wishart (NIW) prior with Minnesota-style shrinkage.
#
#  The Minnesota prior shrinks each unit's VARX* coefficients toward:
#    - A random-walk (own first lag = 1) or white-noise prior
#    - All other coefficients toward zero
#    - Tighter shrinkage for cross-variable, lagged, and foreign variables
#
#  The NIW conjugate posterior has a closed-form solution:
#    Σ_i | Y ~ IW(S_n, v_n)
#    vec(B_i) | Σ_i, Y ~ N(vec(B_n), Σ_i ⊗ V_n)
#
#  This allows direct independent draws (no MCMC burn-in needed).
#  At each draw, unit-level coefficients are drawn, partitioned into
#  A/B/D matrices, and stacked via stack_to_gvar() to produce a global
#  companion matrix.  This generates a posterior distribution over the
#  global GVAR, enabling:
#    - Bayesian GIRFs with credible intervals
#    - Bayesian conditional forecasting
#    - DIC and marginal likelihood computation
#
#  Contents
#  --------
#  1. minnesota_prior()             – set up Minnesota prior parameters
#  2. bayesian_estimate_unit()      – closed-form posterior for one unit
#  3. draw_niw_posterior()          – single draw from NIW posterior
#  4. bayesian_estimate_gvar()      – master: posterior draws + global stacking
#  5. bayesian_girf()               – GIRFs with posterior credible intervals
#  6. bayesian_conditional_forecast() – conditional forecasts from posterior
###############################################################################


# ───────────────────────────────────────────────────────────────────────────────
# 1.  Minnesota Prior Specification
# ───────────────────────────────────────────────────────────────────────────────

#' Set up the Minnesota / Litterman prior for a single unit's VARX* model.
#'
#' The prior on vec(B_i) is N(vec(B_0), V_0) where V_0 is diagonal
#' with entries controlling shrinkage.  The prior on Σ_i is IW(S_0, v_0).
#'
#' Hyperparameters:
#'   λ_1 : overall tightness (typ. 0.1-0.2)
#'   λ_2 : cross-variable relative tightness (typ. 0.5-1.0)
#'   λ_3 : lag decay rate (typ. 1-2)
#'   λ_4 : star variable relative tightness (typ. 0.5-1.0)
#'
#' @param k_dom     Number of domestic variables
#' @param k_star    Number of star variables
#' @param p_lag     Domestic lag order
#' @param q_lag     Foreign lag order
#' @param lambda_1  Overall tightness
#' @param lambda_2  Cross-variable tightness
#' @param lambda_3  Lag decay rate
#' @param lambda_4  Star variable tightness
#' @param sigma_sq  Named numeric vector of univariate AR residual variances
#'                  (length k_dom + k_star). If NULL, all set to 1.
#' @param rw_prior  Logical; if TRUE, own first lag prior mean = 1 (random walk).
#'                  If FALSE, all prior means = 0 (white noise).
#' @param k_global  Number of global exogenous variables (default 0)
#' @param d_lag     Number of global exogenous lags (default 0)
#' @return List with B_0 (m x k_dom), V_0 (m x m diagonal), S_0 (k_dom x k_dom), v_0
minnesota_prior <- function(k_dom, k_star, p_lag, q_lag,
                             lambda_1 = 0.1,
                             lambda_2 = 0.5,
                             lambda_3 = 1.0,
                             lambda_4 = 0.5,
                             sigma_sq = NULL,
                             rw_prior = TRUE,
                             k_global = 0,
                             d_lag = 0,
                             n_det = 1) {

  # Total number of regressors per equation (matches prepare_country_data layout):
  # n_det (deterministic) + k_star (contemp star) + k_global (contemp global)
  # + p_lag * k_dom (lagged domestic) + q_lag * k_star (lagged star)
  # + d_lag * k_global (lagged global)
  m <- n_det + k_star + k_global + p_lag * k_dom + q_lag * k_star + d_lag * k_global

  # Default sigma_sq: univariate AR residual variances
  if (is.null(sigma_sq)) {
    sigma_sq <- rep(1, k_dom + k_star + k_global)
  }
  sigma_dom  <- sigma_sq[1:k_dom]
  sigma_star <- sigma_sq[(k_dom + 1):(k_dom + k_star)]
  sigma_glob <- if (k_global > 0) sigma_sq[(k_dom + k_star + 1):(k_dom + k_star + k_global)] else numeric(0)

  # ---- Prior mean B_0 (m x k_dom) ----
  B_0 <- matrix(0, nrow = m, ncol = k_dom)

  if (rw_prior && p_lag >= 1) {
    # The lag1 domestic block starts at row (n_det + k_star + k_global + 1)
    lag1_start <- n_det + k_star + k_global + 1
    for (j in seq_len(k_dom)) {
      if (lag1_start + (j - 1) <= m) B_0[lag1_start + (j - 1), j] <- 1
    }
  }

  # ---- Prior covariance V_0 (m x m diagonal) ----
  v0_diag <- numeric(m)

  # Deterministic terms: large variance (non-informative)
  if (n_det > 0) {
    for (dd in seq_len(n_det)) {
      v0_diag[dd] <- (100 * lambda_1)^2
    }
  }

  idx <- n_det + 1

  # Contemporaneous star variables
  for (s in seq_len(k_star)) {
    # Treat as lag 0 with star tightness
    for (j in seq_len(k_dom)) {
      v0_diag[idx] <- (lambda_1 * lambda_4)^2 * sigma_dom[j] / sigma_star[s]
    }
    idx <- idx + 1
  }

  # Contemporaneous global exogenous variables
  if (k_global > 0) {
    for (g in seq_len(k_global)) {
      v0_diag[idx] <- (lambda_1 * lambda_4)^2
      idx <- idx + 1
    }
  }

  # Lagged domestic variables: l = 1, ..., p_lag
  for (l in seq_len(p_lag)) {
    for (j in seq_len(k_dom)) {
      # For equation i, coefficient on lag l of variable j:
      # If j == i (own lag): (λ_1 / l^λ_3)^2
      # If j != i (cross):   (λ_1 * λ_2 / l^λ_3)^2 * (σ_i / σ_j)
      # Since V_0 is m x m (shared across equations), we use an average σ_i/σ_j.
      # More precisely, for the diagonal Minnesota prior, each equation shares
      # the same V_0 (applied to each column of B).  So the variance for row
      # corresponding to lag l of variable j is:
      #   own:   (λ_1 / l^λ_3)^2
      #   cross: (λ_1 * λ_2 / l^λ_3)^2
      # We use the "own" entry for the j-th row (will be own for equation j,
      # cross for others).  A common simplification is to use the average.
      # Here we use the cross-variable formula with sigma scaling as an overall entry.
      v0_diag[idx] <- (lambda_1 * lambda_2 / l^lambda_3)^2
      idx <- idx + 1
    }
  }

  # Lagged star variables: l = 1, ..., q_lag
  for (l in seq_len(q_lag)) {
    for (s in seq_len(k_star)) {
      v0_diag[idx] <- (lambda_1 * lambda_4 / (l + 1)^lambda_3)^2
      idx <- idx + 1
    }
  }

  # Lagged global exogenous: l = 1, ..., d_lag
  if (k_global > 0 && d_lag > 0) {
    for (l in seq_len(d_lag)) {
      for (g in seq_len(k_global)) {
        v0_diag[idx] <- (lambda_1 * lambda_4 / (l + 1)^lambda_3)^2
        idx <- idx + 1
      }
    }
  }

  # Fix own-lag diagonal entries (tighter than cross-variable)
  if (p_lag >= 1) {
    lag1_start <- n_det + k_star + k_global + 1
    for (j in seq_len(k_dom)) {
      row_j <- lag1_start + (j - 1)
      v0_diag[row_j] <- (lambda_1 / 1^lambda_3)^2   # own first lag: no lambda_2 penalty
    }
    # For higher lags, own entries
    for (l in 2:p_lag) {
      for (j in seq_len(k_dom)) {
        row_j <- n_det + k_star + k_global + (l - 1) * k_dom + j
        # This is already set to the cross formula; overwrite own-variable entry
        # For the j-th equation, variable j is own. Since V_0 is shared
        # across equations, we keep the cross-variable formula as a compromise.
        # (The equation-specific adjustment comes through the Kronecker product Σ ⊗ V_n.)
      }
    }
  }

  V_0 <- diag(v0_diag)

  # ---- Prior for Σ_i: IW(S_0, v_0) ----
  S_0 <- diag(sigma_dom)   # k_dom x k_dom
  v_0 <- k_dom + 2         # minimally informative

  return(list(
    B_0      = B_0,
    V_0      = V_0,
    S_0      = S_0,
    v_0      = v_0,
    m        = m,
    k_dom    = k_dom,
    lambda_1 = lambda_1,
    lambda_2 = lambda_2,
    lambda_3 = lambda_3,
    lambda_4 = lambda_4
  ))
}


# ───────────────────────────────────────────────────────────────────────────────
# 2.  Bayesian Estimation for a Single Unit
# ───────────────────────────────────────────────────────────────────────────────

#' Compute the closed-form Normal-Inverse-Wishart posterior for one unit.
#'
#' Model: Y = X B + E, E ~ MN(0, I_T, Σ)
#' Prior: vec(B) ~ N(vec(B_0), V_0 ⊗ I_k)  [equation-by-equation],
#'        Σ ~ IW(S_0, v_0)
#'
#' Posterior:
#'   V_n = (V_0^{-1} + X'X)^{-1}
#'   B_n = V_n (V_0^{-1} B_0 + X'Y)
#'   S_n = S_0 + Y'Y + B_0' V_0^{-1} B_0 - B_n' V_n^{-1} B_n
#'   v_n = v_0 + T
#'
#' @param domestic    T x k_dom matrix
#' @param star        T x k_star matrix
#' @param p_lag       Domestic lag order
#' @param q_lag       Foreign lag order
#' @param prior       Output of minnesota_prior(). If NULL, computed internally.
#' @param lambda_1    Overall tightness (used if prior = NULL)
#' @param lambda_2    Cross-variable tightness
#' @param lambda_3    Lag decay
#' @param lambda_4    Star tightness
#' @param rw_prior    Random-walk prior flag
#' @param global_exog T x k_global matrix or NULL
#' @param d_lag       Global exogenous lags
#' @return List with posterior parameters and partitioned coefficients
bayesian_estimate_unit <- function(domestic, star, p_lag, q_lag,
                                    prior = NULL,
                                    lambda_1 = 0.1, lambda_2 = 0.5,
                                    lambda_3 = 1.0, lambda_4 = 0.5,
                                    rw_prior = TRUE,
                                    global_exog = NULL, d_lag = NULL,
                                    deterministic = "intercept") {

  # Prepare regression data using existing machinery
  cd  <- prepare_country_data(domestic, star, p_lag, q_lag,
                               global_exog = global_exog, d_lag = d_lag)
  Y <- cd$Y
  det <- build_deterministic_columns(nrow(Y), deterministic)
  X <- if (!is.null(det$D_mat)) cbind(det$D_mat, cd$X) else cd$X

  k_dom    <- cd$k_dom
  k_star   <- cd$k_star
  k_global <- cd$k_global
  d_lag_i  <- cd$d_lag
  TT       <- nrow(Y)
  m        <- ncol(X)

  # ---- Compute univariate AR residual variances for scaling ----
  sigma_sq <- numeric(k_dom + k_star + k_global)
  all_vars <- cbind(domestic, star)
  if (k_global > 0) all_vars <- cbind(all_vars, global_exog)

  for (v in seq_len(ncol(all_vars))) {
    sv <- all_vars[, v]
    sv <- sv[!is.na(sv)]
    if (length(sv) > 5) {
      # Fit AR(1)
      ar_fit <- tryCatch(stats::ar(sv, aic = FALSE, order.max = 1, method = "ols"),
                          error = function(e) NULL)
      if (!is.null(ar_fit) && !is.null(ar_fit$var.pred)) {
        sigma_sq[v] <- ar_fit$var.pred
      } else {
        sigma_sq[v] <- var(sv, na.rm = TRUE)
      }
    } else {
      sigma_sq[v] <- 1
    }
    sigma_sq[v] <- max(sigma_sq[v], 1e-8)
  }

  # ---- Set up prior ----
  if (is.null(prior)) {
    prior <- minnesota_prior(
      k_dom    = k_dom,
      k_star   = k_star,
      p_lag    = p_lag,
      q_lag    = q_lag,
      lambda_1 = lambda_1,
      lambda_2 = lambda_2,
      lambda_3 = lambda_3,
      lambda_4 = lambda_4,
      sigma_sq = sigma_sq,
      rw_prior = rw_prior,
      k_global = k_global,
      d_lag    = d_lag_i,
      n_det    = det$n_det
    )
  }

  B_0 <- prior$B_0
  V_0 <- prior$V_0
  S_0 <- prior$S_0
  v_0 <- prior$v_0

  # Ensure dimensions match
  if (nrow(B_0) != m || ncol(B_0) != k_dom) {
    # Dimension mismatch: recompute prior with correct dimensions
    prior <- minnesota_prior(
      k_dom = k_dom, k_star = k_star, p_lag = p_lag, q_lag = q_lag,
      lambda_1 = lambda_1, lambda_2 = lambda_2, lambda_3 = lambda_3,
      lambda_4 = lambda_4, sigma_sq = sigma_sq, rw_prior = rw_prior,
      k_global = k_global, d_lag = d_lag_i, n_det = det$n_det
    )
    B_0 <- prior$B_0
    V_0 <- prior$V_0
    S_0 <- prior$S_0
    v_0 <- prior$v_0
  }

  # ---- Compute posterior ----
  # The Minnesota prior ensures V_0 is positive definite (diagonal with positive entries),
  # so V_0_inv and V_0_inv + X'X are always invertible — no fallback needed.
  V_0_inv <- solve(V_0)

  XtX <- crossprod(X)
  XtY <- crossprod(X, Y)

  V_n     <- solve(V_0_inv + XtX)
  B_n     <- V_n %*% (V_0_inv %*% B_0 + XtY)
  V_n_inv <- solve(V_n)

  # Posterior scale for Inverse-Wishart
  S_n <- S_0 + crossprod(Y) + t(B_0) %*% V_0_inv %*% B_0 - t(B_n) %*% V_n_inv %*% B_n

  # Ensure S_n is symmetric and positive-definite
  S_n <- (S_n + t(S_n)) / 2
  eig_S <- eigen(S_n, symmetric = TRUE)
  if (any(eig_S$values < 0)) {
    eig_S$values <- pmax(eig_S$values, 1e-10)
    S_n <- eig_S$vectors %*% diag(eig_S$values) %*% t(eig_S$vectors)
  }

  v_n <- v_0 + TT

  # ---- Posterior mean of Σ ----
  sigma_post <- S_n / max(v_n - k_dom - 1, 1)

  # ---- Partition B_n into deterministic + A, B, D ----
  dp <- partition_deterministic(B_n, det$n_det, det$has_intercept,
                                 det$has_trend, k_dom)
  intercept <- dp$intercept
  trend     <- dp$trend
  idx       <- dp$start_idx

  B_coef <- list()
  B_coef[[1]] <- B_n[idx:(idx + k_star - 1), , drop = FALSE]
  idx <- idx + k_star

  D_coef <- list()
  if (k_global > 0) {
    D_coef[[1]] <- B_n[idx:(idx + k_global - 1), , drop = FALSE]
    idx <- idx + k_global
  }

  A_coef <- list()
  for (l in seq_len(p_lag)) {
    A_coef[[l]] <- B_n[idx:(idx + k_dom - 1), , drop = FALSE]
    idx <- idx + k_dom
  }

  for (l in seq_len(q_lag)) {
    B_coef[[l + 1]] <- B_n[idx:(idx + k_star - 1), , drop = FALSE]
    idx <- idx + k_star
  }

  if (k_global > 0 && d_lag_i > 0) {
    for (l in seq_len(d_lag_i)) {
      D_coef[[l + 1]] <- B_n[idx:(idx + k_global - 1), , drop = FALSE]
      idx <- idx + k_global
    }
  }

  # Residuals at posterior mean
  fitted_vals <- X %*% B_n
  residuals   <- Y - fitted_vals

  return(list(
    # Posterior parameters
    B_n        = B_n,
    V_n        = V_n,
    S_n        = S_n,
    v_n        = v_n,

    # Point estimates (posterior mean)
    intercept  = intercept,
    trend      = trend,
    A          = A_coef,
    B          = B_coef,
    D          = D_coef,
    residuals  = residuals,
    sigma      = sigma_post,
    fitted     = fitted_vals,

    # Metadata
    k_dom      = k_dom,
    k_star     = k_star,
    k_global   = k_global,
    p_lag      = p_lag,
    q_lag      = q_lag,
    d_lag      = d_lag_i,

    # Prior used
    prior      = prior,
    sigma_sq   = sigma_sq
  ))
}


# ───────────────────────────────────────────────────────────────────────────────
# 3.  Draw from NIW Posterior
# ───────────────────────────────────────────────────────────────────────────────

#' Draw a single sample from the Normal-Inverse-Wishart posterior.
#'
#' Σ_draw ~ IW(S_n, v_n)
#' B_draw | Σ_draw ~ MN(B_n, V_n, Σ_draw)
#'
#' Uses an efficient factored draw: B = B_n + chol(V_n) %*% Z %*% chol(Σ)'
#' where Z is m x k_dom standard normal.
#'
#' @param B_n    m x k_dom posterior mean
#' @param V_n    m x m posterior row covariance
#' @param S_n    k_dom x k_dom posterior scale
#' @param v_n    Posterior degrees of freedom
#' @param k_dom  Number of domestic variables
#' @return List with Sigma_draw and B_draw
draw_niw_posterior <- function(B_n, V_n, S_n, v_n, k_dom) {

  m <- nrow(B_n)

  # Draw Σ from Inverse-Wishart(S_n, v_n)
  # IW(S_n, v_n) = solve(Wishart(solve(S_n), v_n))
  S_n_inv <- solve(S_n)
  S_n_inv <- (S_n_inv + t(S_n_inv)) / 2   # symmetrize

  W_draw <- stats::rWishart(1, df = v_n, Sigma = S_n_inv)[, , 1]

  Sigma_draw <- solve(W_draw)
  Sigma_draw <- (Sigma_draw + t(Sigma_draw)) / 2

  # Draw B from MN(B_n, V_n, Σ_draw)
  # Efficient: B = B_n + L_V %*% Z %*% L_Sigma'
  # where L_V = chol(V_n), L_Sigma = chol(Σ_draw), Z ~ MN(0, I, I)

  L_V <- tryCatch(chol(V_n), error = function(e) {
    # V_n may not be PD; regularize
    eig <- eigen(V_n, symmetric = TRUE)
    eig$values <- pmax(eig$values, 1e-10)
    V_reg <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
    chol(V_reg)
  })

  L_Sigma <- tryCatch(chol(Sigma_draw), error = function(e) {
    eig <- eigen(Sigma_draw, symmetric = TRUE)
    eig$values <- pmax(eig$values, 1e-10)
    S_reg <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
    chol(S_reg)
  })

  Z <- matrix(stats::rnorm(m * k_dom), nrow = m, ncol = k_dom)
  B_draw <- B_n + t(L_V) %*% Z %*% L_Sigma

  return(list(
    Sigma_draw = Sigma_draw,
    B_draw     = B_draw
  ))
}


# ───────────────────────────────────────────────────────────────────────────────
# 4.  Master Bayesian GVAR Estimation
# ───────────────────────────────────────────────────────────────────────────────

#' Estimate the Bayesian GVAR: unit-level posterior + global stacking draws.
#'
#' @param gvar_data   Output from prepare_gvar_dataset()
#' @param n_draws     Number of posterior draws to retain
#' @param lambda_1    Overall tightness (default 0.1)
#' @param lambda_2    Cross-variable tightness (default 0.5)
#' @param lambda_3    Lag decay (default 1.0)
#' @param lambda_4    Star variable tightness (default 0.5)
#' @param rw_prior    Logical; random-walk prior (default TRUE)
#' @param seed        Random seed
#' @param regularise  Logical; if TRUE use rejection sampling to retain only
#'                    stable draws (max eigenvalue modulus < 1).  The sampler
#'                    will attempt up to \code{max_attempts} draws to collect
#'                    \code{n_draws} stable ones.  Unstable draws are discarded.
#'                    Useful when default lambdas yield too many explosive draws.
#'                    Default FALSE (backward compatible).
#' @param max_attempts Maximum total draw attempts when regularise = TRUE.
#'                    Defaults to 10 * n_draws.  A warning is issued if fewer
#'                    than n_draws stable draws are obtained within this budget.
#' @return A bayesian_gvar object
bayesian_estimate_gvar <- function(gvar_data,
                                    n_draws = 1000,
                                    lambda_1 = 0.1,
                                    lambda_2 = 0.5,
                                    lambda_3 = 1.0,
                                    lambda_4 = 0.5,
                                    rw_prior = TRUE,
                                    seed = 42,
                                    deterministic = "intercept",
                                    regularise = FALSE,
                                    max_attempts = NULL) {

  print_banner("Bayesian GVAR Estimation (Minnesota Prior)")

  set.seed(seed)

  unit_names <- gvar_data$unit_names
  N <- length(unit_names)

  # ---- Stage 1: Compute posterior parameters for each unit ----
  message("  Stage 1: Computing unit-level posteriors ...")

  unit_posteriors <- setNames(vector("list", N), unit_names)

  for (u in unit_names) {
    cd <- gvar_data$unit_data[[u]]

    # Determine global exogenous
    g_exog  <- NULL
    g_d_lag <- NULL
    if (!is.null(gvar_data$global_config) &&
        u != gvar_data$global_config$dominant_unit) {
      g_exog  <- gvar_data$global_config$global_series
      g_d_lag <- gvar_data$global_config$d_lag
    }

    # We need the original domestic and star matrices
    # Reconstruct from the data preparation
    # cd contains Y (trimmed), X (regressors), k_dom, k_star, etc.
    # For bayesian_estimate_unit we need the full domestic and star
    # But actually, bayesian_estimate_unit calls prepare_country_data internally
    # which needs the full T series.
    # WORKAROUND: pass Y and X directly instead of calling prepare_country_data.

    # Direct posterior computation using already-prepared Y and X:
    Y <- cd$Y
    det <- build_deterministic_columns(nrow(Y), deterministic)
    X <- if (!is.null(det$D_mat)) cbind(det$D_mat, cd$X) else cd$X

    k_dom    <- cd$k_dom
    k_star   <- cd$k_star
    k_global <- if (!is.null(cd$k_global)) cd$k_global else 0
    d_lag_i  <- if (!is.null(cd$d_lag)) cd$d_lag else 0
    n_det    <- det$n_det
    TT       <- nrow(Y)
    m        <- ncol(X)

    # Compute sigma_sq from residuals of a simple OLS
    X_ols <- if (!is.null(det$D_mat)) cbind(det$D_mat, cd$X) else cd$X
    ols_fit <- ols_estimate(cd$Y, X_ols, intercept = FALSE)
    sigma_sq <- numeric(k_dom + k_star + k_global)
    for (v in seq_len(k_dom)) {
      sigma_sq[v] <- max(var(ols_fit$residuals[, v]), 1e-8)
    }
    for (v in (k_dom + 1):length(sigma_sq)) {
      sigma_sq[v] <- 1   # default for star/global
    }

    # Minnesota prior
    prior <- minnesota_prior(
      k_dom = k_dom, k_star = k_star, p_lag = cd$p_lag, q_lag = cd$q_lag,
      lambda_1 = lambda_1, lambda_2 = lambda_2, lambda_3 = lambda_3,
      lambda_4 = lambda_4, sigma_sq = sigma_sq, rw_prior = rw_prior,
      k_global = k_global, d_lag = d_lag_i, n_det = n_det
    )

    B_0 <- prior$B_0
    V_0 <- prior$V_0
    S_0 <- prior$S_0
    v_0 <- prior$v_0

    # Ensure prior dimensions match data
    if (nrow(B_0) != m) {
      # Mismatch: create simple diffuse prior
      message(sprintf("  [%s] Dimension mismatch (prior %d vs data %d) – using diffuse prior.",
                      u, nrow(B_0), m))
      B_0 <- matrix(0, nrow = m, ncol = k_dom)
      if (rw_prior && cd$p_lag >= 1) {
        lag1_start <- n_det + k_star + k_global + 1
        for (j in seq_len(k_dom)) {
          if (lag1_start + j - 1 <= m) B_0[lag1_start + j - 1, j] <- 1
        }
      }
      V_0 <- diag(lambda_1^2, m)
      S_0 <- diag(sigma_sq[1:k_dom])
      v_0 <- k_dom + 2
    }

    # Posterior computation — Minnesota prior ensures V_0 is PD, so no fallback needed
    V_0_inv <- solve(V_0)
    XtX <- crossprod(X)
    XtY <- crossprod(X, Y)

    V_n     <- solve(V_0_inv + XtX)
    B_n     <- V_n %*% (V_0_inv %*% B_0 + XtY)
    V_n_inv <- solve(V_n)
    S_n <- S_0 + crossprod(Y) + t(B_0) %*% V_0_inv %*% B_0 - t(B_n) %*% V_n_inv %*% B_n
    S_n <- (S_n + t(S_n)) / 2
    eig_S <- eigen(S_n, symmetric = TRUE)
    if (any(eig_S$values < 0)) {
      eig_S$values <- pmax(eig_S$values, 1e-10)
      S_n <- eig_S$vectors %*% diag(eig_S$values) %*% t(eig_S$vectors)
    }

    v_n <- v_0 + TT

    message(sprintf("  [%s] Posterior: v_n = %d, m = %d, k_dom = %d",
                    u, v_n, m, k_dom))

    unit_posteriors[[u]] <- list(
      B_n = B_n, V_n = V_n, S_n = S_n, v_n = v_n,
      k_dom = k_dom, k_star = k_star, k_global = k_global,
      p_lag = cd$p_lag, q_lag = cd$q_lag, d_lag = d_lag_i,
      prior = prior, sigma_sq = sigma_sq
    )
  }

  # ---- Stage 2: Point estimate (posterior mean) ----
  message("\n  Stage 2: Computing point estimate (posterior mean) ...")

  unit_fits_mean <- setNames(vector("list", N), unit_names)

  for (u in unit_names) {
    post <- unit_posteriors[[u]]
    B_n  <- post$B_n

    # Partition B_n using helpers
    k_dom <- post$k_dom; k_star <- post$k_star; k_global <- post$k_global
    p_lag <- post$p_lag; q_lag <- post$q_lag; d_lag_i <- post$d_lag

    det <- build_deterministic_columns(1, deterministic)  # just for metadata
    dp <- partition_deterministic(B_n, det$n_det, det$has_intercept,
                                   det$has_trend, k_dom)
    intercept <- dp$intercept
    trend     <- dp$trend
    idx       <- dp$start_idx

    B_coef <- list()
    B_coef[[1]] <- B_n[idx:(idx + k_star - 1), , drop = FALSE]
    idx <- idx + k_star

    D_coef <- list()
    if (k_global > 0) {
      D_coef[[1]] <- B_n[idx:(idx + k_global - 1), , drop = FALSE]
      idx <- idx + k_global
    }

    A_coef <- list()
    for (l in seq_len(p_lag)) {
      A_coef[[l]] <- B_n[idx:(idx + k_dom - 1), , drop = FALSE]
      idx <- idx + k_dom
    }
    for (l in seq_len(q_lag)) {
      B_coef[[l + 1]] <- B_n[idx:(idx + k_star - 1), , drop = FALSE]
      idx <- idx + k_star
    }
    if (k_global > 0 && d_lag_i > 0) {
      for (l in seq_len(d_lag_i)) {
        D_coef[[l + 1]] <- B_n[idx:(idx + k_global - 1), , drop = FALSE]
        idx <- idx + k_global
      }
    }

    sigma_post <- post$S_n / max(post$v_n - k_dom - 1, 1)

    # Compute residuals from posterior mean
    cd <- gvar_data$unit_data[[u]]
    det_full <- build_deterministic_columns(nrow(cd$Y), deterministic)
    X_full <- if (!is.null(det_full$D_mat)) cbind(det_full$D_mat, cd$X) else cd$X
    residuals <- cd$Y - X_full %*% B_n

    unit_fits_mean[[u]] <- list(
      intercept = intercept, trend = trend,
      A = A_coef, B = B_coef, D = D_coef,
      residuals = residuals, sigma = sigma_post, fitted = X_full %*% B_n,
      k_dom = k_dom, k_star = k_star, k_global = k_global,
      p_lag = p_lag, q_lag = q_lag, d_lag = d_lag_i
    )
  }

  # Stack point estimate
  data_list_for_stack <- list()
  for (u in unit_names) {
    data_list_for_stack[[u]] <- gvar_data$unit_data[[u]]$Y
  }

  point_estimate <- stack_to_gvar(unit_fits_mean, data_list_for_stack, gvar_data$W,
                                   global_config = gvar_data$global_config)
  point_estimate <- c(point_estimate, list(
    unit_fits = unit_fits_mean, unit_names = unit_names,
    freq = gvar_data$freq, W = gvar_data$W
  ))
  class(point_estimate) <- "gvar_model"

  k_total  <- point_estimate$k_total
  p_global <- point_estimate$p_global
  kp       <- k_total * p_global

  # ---- Stage 3: Posterior draws ----
  if (regularise) {
    message(sprintf("\n  Stage 3: Drawing posterior samples (regularised – retaining only stable draws) ..."))
  } else {
    message(sprintf("\n  Stage 3: Drawing %d posterior samples ...", n_draws))
  }

  if (is.null(max_attempts)) max_attempts <- 10L * n_draws

  companion_draws <- array(NA, dim = c(kp, kp, n_draws))
  Sigma_draws     <- array(NA, dim = c(k_total, k_total, n_draws))
  F0_draws        <- matrix(NA, nrow = k_total, ncol = n_draws)
  stable_flags    <- logical(n_draws)

  # Helper: draw one GVAR sample from unit posteriors
  .one_draw <- function() {
    unit_fits_draw <- setNames(vector("list", N), unit_names)

    for (u in unit_names) {
      post  <- unit_posteriors[[u]]
      draw  <- draw_niw_posterior(post$B_n, post$V_n, post$S_n, post$v_n, post$k_dom)
      B_d   <- draw$B_draw
      k_dom <- post$k_dom; k_star <- post$k_star; k_global <- post$k_global
      p_lag <- post$p_lag; q_lag  <- post$q_lag;  d_lag_i  <- post$d_lag

      det_d <- build_deterministic_columns(1, deterministic)
      dp_d  <- partition_deterministic(B_d, det_d$n_det, det_d$has_intercept,
                                        det_d$has_trend, k_dom)
      intercept_d <- dp_d$intercept
      trend_d     <- dp_d$trend
      idx         <- dp_d$start_idx

      B_coef_d <- list()
      B_coef_d[[1]] <- B_d[idx:(idx + k_star - 1), , drop = FALSE]
      idx <- idx + k_star

      D_coef_d <- list()
      if (k_global > 0) {
        D_coef_d[[1]] <- B_d[idx:(idx + k_global - 1), , drop = FALSE]
        idx <- idx + k_global
      }

      A_coef_d <- list()
      for (l in seq_len(p_lag)) {
        A_coef_d[[l]] <- B_d[idx:(idx + k_dom - 1), , drop = FALSE]
        idx <- idx + k_dom
      }
      for (l in seq_len(q_lag)) {
        B_coef_d[[l + 1]] <- B_d[idx:(idx + k_star - 1), , drop = FALSE]
        idx <- idx + k_star
      }
      if (k_global > 0 && d_lag_i > 0) {
        for (l in seq_len(d_lag_i)) {
          D_coef_d[[l + 1]] <- B_d[idx:(idx + k_global - 1), , drop = FALSE]
          idx <- idx + k_global
        }
      }

      unit_fits_draw[[u]] <- list(
        intercept = intercept_d, trend = trend_d,
        A = A_coef_d, B = B_coef_d, D = D_coef_d,
        sigma = draw$Sigma_draw,
        residuals = unit_fits_mean[[u]]$residuals,
        k_dom = k_dom, k_star = k_star, k_global = k_global,
        p_lag = p_lag, q_lag = q_lag, d_lag = d_lag_i
      )
    }

    tryCatch(
      stack_to_gvar(unit_fits_draw, data_list_for_stack, gvar_data$W,
                     global_config = gvar_data$global_config),
      error = function(e) NULL
    )
  }

  pb_interval <- max(1, floor(n_draws / 10))

  if (!regularise) {
    # Standard loop: collect exactly n_draws (stable or not)
    for (s in seq_len(n_draws)) {
      if (s %% pb_interval == 0)
        message(sprintf("    Draw %d / %d ...", s, n_draws))

      gd <- .one_draw()
      if (!is.null(gd)) {
        companion_draws[, , s] <- gd$companion
        Sigma_draws[, , s]     <- gd$Sigma
        F0_draws[, s]          <- gd$F0
        stable_flags[s]        <- gd$is_stable
      }
    }
  } else {
    # Rejection sampling: only retain stable draws
    s        <- 0L      # number of retained draws
    attempts <- 0L

    while (s < n_draws && attempts < max_attempts) {
      attempts <- attempts + 1L
      gd <- .one_draw()
      if (is.null(gd) || !isTRUE(gd$is_stable)) next

      s <- s + 1L
      companion_draws[, , s] <- gd$companion
      Sigma_draws[, , s]     <- gd$Sigma
      F0_draws[, s]          <- gd$F0
      stable_flags[s]        <- TRUE

      if (s %% pb_interval == 0)
        message(sprintf("    Stable draw %d / %d (attempts so far: %d) ...",
                        s, n_draws, attempts))
    }

    if (s < n_draws) {
      warning(sprintf(
        "[Bayesian] regularise=TRUE: only %d / %d stable draws obtained after %d attempts.\n  Consider tightening lambda_1 or increasing max_attempts.",
        s, n_draws, attempts))
      # Trim arrays to the draws actually obtained
      if (s == 0) {
        companion_draws <- array(NA, dim = c(kp, kp, 0))
        Sigma_draws     <- array(NA, dim = c(k_total, k_total, 0))
        F0_draws        <- matrix(NA, nrow = k_total, ncol = 0)
        stable_flags    <- logical(0)
      } else {
        companion_draws <- companion_draws[, , seq_len(s), drop = FALSE]
        Sigma_draws     <- Sigma_draws[, , seq_len(s), drop = FALSE]
        F0_draws        <- F0_draws[, seq_len(s), drop = FALSE]
        stable_flags    <- stable_flags[seq_len(s)]
      }
    }
    message(sprintf("    Rejection sampling: %d stable draws from %d attempts (%.1f%% acceptance).",
                    s, attempts, 100 * s / max(attempts, 1)))
  }

  n_stable <- sum(stable_flags, na.rm = TRUE)
  n_valid  <- sum(!is.na(stable_flags))

  message(sprintf("\n  [Bayesian] %d / %d draws are stable (%d valid).",
                  n_stable, n_draws, n_valid))

  result <- list(
    point_estimate  = point_estimate,
    draws           = list(
      companion     = companion_draws,
      Sigma         = Sigma_draws,
      F0            = F0_draws,
      is_stable     = stable_flags
    ),
    unit_posteriors = unit_posteriors,
    n_draws         = n_draws,
    n_stable        = n_stable,
    var_names       = point_estimate$var_names,
    k_total         = k_total,
    p_global        = p_global,
    lambda          = list(lambda_1 = lambda_1, lambda_2 = lambda_2,
                           lambda_3 = lambda_3, lambda_4 = lambda_4),
    rw_prior        = rw_prior,
    deterministic   = deterministic
  )

  class(result) <- "bayesian_gvar"
  return(result)
}


# ───────────────────────────────────────────────────────────────────────────────
# 5.  Bayesian GIRFs with Posterior Credible Intervals
# ───────────────────────────────────────────────────────────────────────────────

#' Compute Generalised IRFs with posterior credible intervals.
#'
#' At each stable posterior draw, computes the GIRF using that draw's
#' companion matrix and Σ.  The posterior distribution of the IRF is
#' summarised by the median and percentile-based credible bands.
#'
#' @param bayesian_gvar  Output from bayesian_estimate_gvar()
#' @param shock_var      Index or name of the shocked variable
#' @param horizon        Number of periods
#' @param ci_level       Credible interval level (default 0.90)
#' @return List with point, median, lower, upper matrices
bayesian_girf <- function(bayesian_gvar, shock_var, horizon = 20,
                           ci_level = 0.90) {

  print_banner("Bayesian GIRFs")

  pe        <- bayesian_gvar$point_estimate
  k_total   <- bayesian_gvar$k_total
  p_global  <- bayesian_gvar$p_global
  var_names <- bayesian_gvar$var_names
  kp        <- k_total * p_global

  # Resolve shock variable
  if (is.character(shock_var)) {
    shock_idx <- match(shock_var, var_names)
    if (is.na(shock_idx)) stop("Shock variable '", shock_var, "' not found.")
  } else {
    shock_idx <- shock_var
  }

  # Point estimate GIRF
  point_irf <- girf_analytic(pe, shock_var = shock_idx, horizon = horizon)

  # Collect draws
  stable_idx <- which(bayesian_gvar$draws$is_stable)
  n_use <- length(stable_idx)

  if (n_use < 10) {
    warning("[Bayesian GIRF] Fewer than 10 stable draws available.")
  }

  irf_array <- array(NA, dim = c(horizon + 1, k_total, n_use))

  message(sprintf("  Computing GIRFs across %d stable draws ...", n_use))

  for (i in seq_len(n_use)) {
    s <- stable_idx[i]

    comp_s  <- bayesian_gvar$draws$companion[, , s]
    Sigma_s <- bayesian_gvar$draws$Sigma[, , s]

    # Build temporary model for girf_analytic
    tmp_model <- list(
      companion = comp_s,
      Sigma     = Sigma_s,
      k_total   = k_total,
      p_global  = p_global,
      var_names = var_names
    )
    class(tmp_model) <- "gvar_model"

    irf_s <- tryCatch(
      girf_analytic(tmp_model, shock_var = shock_idx, horizon = horizon),
      error = function(e) NULL
    )

    if (!is.null(irf_s)) {
      irf_array[, , i] <- irf_s
    }
  }

  # Compute summary statistics
  alpha <- 1 - ci_level
  lower <- apply(irf_array, c(1, 2), quantile, probs = alpha / 2, na.rm = TRUE)
  upper <- apply(irf_array, c(1, 2), quantile, probs = 1 - alpha / 2, na.rm = TRUE)
  median_irf <- apply(irf_array, c(1, 2), median, na.rm = TRUE)

  colnames(point_irf) <- var_names
  colnames(median_irf) <- var_names
  colnames(lower) <- var_names
  colnames(upper) <- var_names

  message("[Bayesian GIRF] Complete.")

  return(list(
    point     = point_irf,
    median    = median_irf,
    lower     = lower,
    upper     = upper,
    ci_level  = ci_level,
    n_draws   = n_use,
    shock_var = if (is.character(shock_var)) shock_var else var_names[shock_idx],
    var_names = var_names
  ))
}


# ───────────────────────────────────────────────────────────────────────────────
# 6.  Bayesian Conditional Forecasting
# ───────────────────────────────────────────────────────────────────────────────

#' Conditional forecasting from the Bayesian posterior.
#'
#' At each posterior draw, constructs a gvar_model and applies the chosen
#' conditional forecasting method (Kalman or WZ).  The posterior distribution
#' of the forecast is summarised by credible intervals.
#'
#' @param bayesian_gvar  Output from bayesian_estimate_gvar()
#' @param conditions     Data frame with variable, horizon, value, type, tolerance
#' @param max_h          Forecast horizon (default: max horizon in conditions + 4)
#' @param data_list      Named list of T x k_i matrices (for initial state)
#' @param method         "kalman" or "wz"
#' @param ci_level       Credible interval level
#' @param n_cf_draws     Number of posterior draws to use
#' @return List with forecast_mean, forecast_lower, forecast_upper, etc.
bayesian_conditional_forecast <- function(bayesian_gvar,
                                           conditions,
                                           max_h = NULL,
                                           data_list = NULL,
                                           method = c("kalman", "wz"),
                                           ci_level = 0.90,
                                           n_cf_draws = 500) {

  print_banner("Bayesian Conditional Forecasting")

  method <- match.arg(method)

  if (is.null(max_h)) {
    max_h <- max(conditions$horizon) + 4
  }

  k_total   <- bayesian_gvar$k_total
  var_names <- bayesian_gvar$var_names

  # Select stable draws
  stable_idx <- which(bayesian_gvar$draws$is_stable)
  n_use <- min(n_cf_draws, length(stable_idx))
  use_idx <- sample(stable_idx, n_use)

  message(sprintf("  Using %d posterior draws with method '%s' ...", n_use, method))

  forecast_array <- array(NA, dim = c(max_h, k_total, n_use))

  for (i in seq_len(n_use)) {
    s <- use_idx[i]

    # Build gvar_model for this draw
    comp_s  <- bayesian_gvar$draws$companion[, , s]
    Sigma_s <- bayesian_gvar$draws$Sigma[, , s]
    F0_s    <- bayesian_gvar$draws$F0[, s]
    p_global <- bayesian_gvar$p_global

    # Reconstruct F matrices from companion
    F_list <- list()
    for (l in seq_len(p_global)) {
      cols <- ((l - 1) * k_total + 1):(l * k_total)
      F_list[[l]] <- comp_s[1:k_total, cols]
    }

    tmp_model <- list(
      F0        = F0_s,
      F         = F_list,
      Sigma     = Sigma_s,
      companion = comp_s,
      k_total   = k_total,
      p_global  = p_global,
      var_names = var_names
    )
    class(tmp_model) <- "gvar_model"

    # Run conditional forecast
    cf_s <- tryCatch({
      if (method == "kalman") {
        conditional_forecast(
          gvar_model = tmp_model,
          conditions = conditions,
          max_h      = max_h,
          data_list  = data_list,
          ci_level   = ci_level
        )
      } else {
        conditional_forecast_wz(
          gvar_model = tmp_model,
          conditions = conditions,
          max_h      = max_h,
          data_list  = data_list,
          n_draws    = 200,
          ci_level   = ci_level
        )
      }
    }, error = function(e) NULL)

    if (!is.null(cf_s)) {
      fmean <- if (method == "kalman") cf_s$forecast_mean else cf_s$forecast_point
      if (!is.null(fmean)) {
        n_rows <- min(nrow(fmean), max_h)
        forecast_array[1:n_rows, , i] <- fmean[1:n_rows, ]
      }
    }

    if (i %% max(1, floor(n_use / 5)) == 0) {
      message(sprintf("    %d / %d draws processed ...", i, n_use))
    }
  }

  # Compute summary statistics
  alpha <- 1 - ci_level
  forecast_mean  <- apply(forecast_array, c(1, 2), mean, na.rm = TRUE)
  forecast_lower <- apply(forecast_array, c(1, 2), quantile, probs = alpha / 2, na.rm = TRUE)
  forecast_upper <- apply(forecast_array, c(1, 2), quantile, probs = 1 - alpha / 2, na.rm = TRUE)
  forecast_sd    <- apply(forecast_array, c(1, 2), sd, na.rm = TRUE)

  colnames(forecast_mean)  <- var_names
  colnames(forecast_lower) <- var_names
  colnames(forecast_upper) <- var_names
  colnames(forecast_sd)    <- var_names

  message("[Bayesian CF] Complete.")

  return(list(
    forecast_mean  = forecast_mean,
    forecast_lower = forecast_lower,
    forecast_upper = forecast_upper,
    forecast_sd    = forecast_sd,
    conditions     = conditions,
    method         = method,
    ci_level       = ci_level,
    var_names      = var_names,
    n_draws_used   = n_use
  ))
}
