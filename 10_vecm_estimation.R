###############################################################################
#  10_vecm_estimation.R  –  VECM Estimation for the GVAR Framework
#
#  When unit-root tests indicate I(1) variables and Johansen tests detect
#  cointegration (r_i > 0), the appropriate specification is a VECX*(p, q)
#  model in error-correction form.
#
#  The VECX* model for unit i:
#
#    Δx_{it} = c_i + α_i β_i' z_{i,t-1}
#              + Σ_{l=1}^{p-1} Γ_{il}^x Δx_{i,t-l}
#              + Σ_{l=0}^{q-1} Γ_{il}^* Δx*_{i,t-l}  +  u_{it}
#
#  where z_{it} = (x_{it}', x*_{it}')' , α_i is k_i × r_i (loadings),
#  β_i is (k_i + k*_i) × r_i (cointegrating vectors).
#
#  After estimation the VECM is converted to a levels-VAR representation
#  so that the existing stack_to_gvar() machinery (and all downstream
#  IRF / Kalman / WZ code) works unchanged.
#
#  Contents
#  --------
#  1. estimate_unit_vecx()    – Johansen-based VECX* estimation for one unit
#  2. vecm_to_levels_var()    – convert VECM parameters to levels-VAR form
#  3. estimate_gvecm()        – master function: VECX* + global stacking
###############################################################################


# ───────────────────────────────────────────────────────────────────────────────
# 1.  Estimate a Single VECX* Model
# ───────────────────────────────────────────────────────────────────────────────

#' Estimate the VECX*(p, q) model for a single unit using the
#' Pesaran-Shin-Smith (2000) conditional approach.
#'
#' Steps:
#'   1. Run Johansen on the full z = (x, x*) system to identify beta.
#'   2. Compute error-correction terms: ect_{t} = beta' z_{t}.
#'   3. Estimate the conditional VECM by OLS with Δx_{it} on the LHS
#'      and (ect_{t-1}, ΔX lags, ΔX* lags) on the RHS.
#'   4. Convert to levels-VAR form for stack_to_gvar() compatibility.
#'
#' @param domestic   T × k_i  matrix of domestic variables
#' @param star       T × k*   matrix of star variables
#' @param rank       Integer r_i; cointegrating rank from johansen_test_single()
#' @param p_lag      Integer; domestic lag order (in levels; VECM uses p-1 diffs)
#' @param q_lag      Integer; foreign lag order (in levels; VECM uses q diffs incl contemp)
#' @param ecdet      Deterministic specification: "const", "trend", or "none"
#' @param global_exog  T × g matrix of global exogenous variables (NULL if none)
#' @param d_lag      Integer; lags for global exogenous (NULL defaults to q_lag)
#' @return  A list with VECM-specific and levels-VAR-compatible components
estimate_unit_vecx <- function(domestic, star, rank, p_lag, q_lag,
                               ecdet = "const",
                               global_exog = NULL, d_lag = NULL,
                               deterministic = "intercept") {
  
  domestic <- as.matrix(domestic)
  star     <- as.matrix(star)
  TT       <- nrow(domestic)
  k_dom    <- ncol(domestic)
  k_star   <- ncol(star)
  m        <- k_dom + k_star
  
  # Global exogenous
  k_global <- 0
  if (!is.null(global_exog)) {
    global_exog <- as.matrix(global_exog)
    k_global <- ncol(global_exog)
    if (is.null(d_lag)) d_lag <- q_lag
  } else {
    d_lag <- 0
  }
  
  # Cap rank: rank must be strictly less than m (full rank = no cointegration)
  if (rank >= m) {
    message(sprintf("  [VECX*] rank (%d) >= system dimension m (%d) – capping at m-1 = %d.",
                    rank, m, m - 1))
    rank <- m - 1
  }
  
  # ---- Step 1: Johansen to get beta (cointegrating vectors) ----
  z_mat <- cbind(domestic, star)
  K_jo  <- max(p_lag, q_lag) + 1   # VAR order in levels for Johansen
  K_jo  <- max(2, min(K_jo, floor(TT / (m + 2))))
  
  ca_obj <- urca::ca.jo(z_mat, type = "trace", ecdet = ecdet, K = K_jo,
                        spec = "transitory")
  
  # Extract beta: includes deterministic rows when ecdet != "none"
  # ca.jo with ecdet="const" adds a restricted-constant row to V;
  # ca.jo with ecdet="trend" adds a restricted-trend row.
  beta_full <- ca_obj@V[, 1:rank, drop = FALSE]
  
  # Separate variable part (m rows) from deterministic part
  beta <- beta_full[1:m, , drop = FALSE]   # m × rank  (variable coefficients)
  beta_det <- NULL
  if (ecdet == "const") {
    beta_det <- beta_full[m + 1, , drop = FALSE]  # 1 × rank
    z_aug <- cbind(z_mat, 1)
  } else if (ecdet == "trend") {
    beta_det <- beta_full[m + 1, , drop = FALSE]
    z_aug <- cbind(z_mat, seq_len(TT))
  } else {
    z_aug <- z_mat
  }
  
  # ---- Step 2: Compute error-correction terms ----
  # ect_t = z_aug_t' beta_full  =>  T × rank  (includes deterministic component)
  ect <- z_aug %*% beta_full
  
  # ---- Step 3: Build VECM regressors and estimate by OLS ----
  # First-differences
  d_dom  <- diff(domestic)        # (T-1) × k_dom
  d_star <- diff(star)            # (T-1) × k_star
  d_glob <- NULL
  if (k_global > 0) {
    d_glob <- diff(global_exog)   # (T-1) × k_global
  }
  
  # Effective number of lagged differences for domestic and star
  p_diff <- max(p_lag - 1, 0)     # VECM has p-1 lagged Δx
  q_diff <- max(q_lag - 1, 0)     # VECM has q-1 lagged Δx* (plus contemporaneous)
  d_diff <- max(d_lag - 1, 0)
  
  max_lag_vecm <- max(1, p_diff, q_diff + 1, d_diff + 1)  # +1 for contemp star diff
  
  # Time indices: we need observations from max_lag_vecm+1 to T-1 (in diff space)
  # In diff space, row t corresponds to original time t+1
  T_diff <- nrow(d_dom)   # = TT - 1
  t_start <- max_lag_vecm + 1
  t_end   <- T_diff
  n_eff   <- t_end - t_start + 1
  
  if (n_eff < m + rank + 5) {
    warning(sprintf("Unit VECX*: Only %d effective obs (need more for reliable estimation).", n_eff))
  }
  
  # LHS: Δx_{it}
  Y_vecm <- d_dom[t_start:t_end, , drop = FALSE]
  
  # RHS components
  X_parts <- list()
  
  # (a) Error-correction terms: ect_{t-1} in original time = ect at row (t_start-1):(t_end-1) in levels
  #     In diff indexing, diff row t corresponds to levels row t+1.
  #     ect_{t-1} for diff row t = ect at levels row t.
  #     So for diff rows t_start:t_end, ect_{t-1} = ect[(t_start):(t_end), ]
  ect_lagged <- ect[t_start:(t_end), , drop = FALSE]
  colnames(ect_lagged) <- paste0("ect", 1:rank)
  X_parts[["ect"]] <- ect_lagged
  
  # (b) Contemporaneous Δx*_{it}
  X_parts[["dstar_contemp"]] <- d_star[t_start:t_end, , drop = FALSE]
  
  # (c) Lagged Δx_{i,t-l} for l = 1, ..., p_diff
  if (p_diff > 0) {
    for (l in seq_len(p_diff)) {
      block <- d_dom[(t_start - l):(t_end - l), , drop = FALSE]
      colnames(block) <- paste0(colnames(domestic), "_dlag", l)
      X_parts[[paste0("ddom_lag", l)]] <- block
    }
  }
  
  # (d) Lagged Δx*_{i,t-l} for l = 1, ..., q_diff
  if (q_diff > 0) {
    for (l in seq_len(q_diff)) {
      block <- d_star[(t_start - l):(t_end - l), , drop = FALSE]
      colnames(block) <- paste0(colnames(star), "_dlag", l)
      X_parts[[paste0("dstar_lag", l)]] <- block
    }
  }
  
  # (e) Global exogenous: contemporaneous Δd_{it} and lags
  if (k_global > 0) {
    X_parts[["dglob_contemp"]] <- d_glob[t_start:t_end, , drop = FALSE]
    if (d_diff > 0) {
      for (l in seq_len(d_diff)) {
        block <- d_glob[(t_start - l):(t_end - l), , drop = FALSE]
        colnames(block) <- paste0(colnames(global_exog), "_dlag", l)
        X_parts[[paste0("dglob_lag", l)]] <- block
      }
    }
  }
  
  X_vecm <- do.call(cbind, X_parts)
  
  # ---- Deterministic terms in the VECM short-run equation ----
  # When ecdet = "const", the restricted constant is already embedded in the ECT
  # (via z_aug = [z_mat, 1] %*% beta_full).  Adding an unrestricted intercept on
  # top creates collinearity when rank is high.  Same logic for ecdet = "trend".
  #
  # Johansen's five cases:
  #   ecdet="none"  + deterministic="none"      → Case 1: no deterministic terms
  #   ecdet="const" + deterministic="none"       → Case 2: restricted intercept only
  #   ecdet="const" + deterministic="intercept"  → Case 3: restricted + unrestricted
  #   ecdet="trend" + deterministic="intercept"  → Case 4: restricted trend + unrest. intercept
  #   ecdet="trend" + deterministic="both"        → Case 5: restricted trend + unrest. both
  #
  # To avoid collinearity, strip the deterministic component that duplicates ecdet:
  vecm_deterministic <- deterministic
  if (ecdet == "const" && deterministic %in% c("intercept", "both")) {
    # The restricted constant in ECT subsumes the intercept role;
    # keep only the trend part if requested
    vecm_deterministic <- if (deterministic == "both") "trend" else "none"
  } else if (ecdet == "trend" && deterministic %in% c("trend", "both")) {
    # The restricted trend in ECT subsumes the trend role;
    # keep only the intercept part if requested
    vecm_deterministic <- if (deterministic == "both") "intercept" else "none"
  }
  
  det_info <- build_deterministic_columns(nrow(Y_vecm), vecm_deterministic)
  X_full <- if (!is.null(det_info$D_mat)) cbind(det_info$D_mat, X_vecm) else X_vecm
  fit <- ols_estimate(Y_vecm, X_full, intercept = FALSE)
  
  # ---- Step 4: Extract VECM parameters ----
  beta_hat <- fit$beta
  
  # Verify dimensions: expected rows = n_det + rank + k_star + p_diff*k_dom + q_diff*k_star + k_global*(1 + d_diff)
  n_expected <- det_info$n_det + rank + k_star + p_diff * k_dom + q_diff * k_star +
    k_global * (1 + d_diff)
  if (nrow(beta_hat) != n_expected) {
    stop(sprintf("VECX* coefficient dimension mismatch: beta_hat has %d rows, expected %d (n_det=%d, rank=%d, k_star=%d, p_diff=%d, k_dom=%d, q_diff=%d, k_global=%d, d_diff=%d)",
                 nrow(beta_hat), n_expected, det_info$n_det, rank, k_star, p_diff, k_dom, q_diff, k_global, d_diff))
  }
  
  dp <- partition_deterministic(beta_hat, det_info$n_det,
                                det_info$has_intercept,
                                det_info$has_trend, k_dom)
  intercept_vecm <- dp$intercept
  trend_vecm     <- dp$trend
  idx            <- dp$start_idx
  
  # Alpha (loading matrix): coefficients on ect_{t-1}
  alpha_block <- beta_hat[idx:(idx + rank - 1), , drop = FALSE]   # rank × k_dom
  alpha <- t(alpha_block)   # k_dom × rank
  idx <- idx + rank
  
  # Levels-VAR intercept: VECM intercept + restricted-constant contribution
  intercept_levels_final <- intercept_vecm
  if (!is.null(beta_det)) {
    intercept_levels_final <- intercept_levels_final + as.vector(alpha %*% t(beta_det))
  }
  
  # Contemporaneous Δx* coefficients
  Gamma_star_0 <- beta_hat[idx:(idx + k_star - 1), , drop = FALSE]   # k_star × k_dom
  idx <- idx + k_star
  
  # Lagged Δx domestic coefficients
  Gamma_dom <- list()
  if (p_diff > 0) {
    for (l in seq_len(p_diff)) {
      Gamma_dom[[l]] <- beta_hat[idx:(idx + k_dom - 1), , drop = FALSE]   # k_dom × k_dom
      idx <- idx + k_dom
    }
  }
  
  # Lagged Δx* coefficients
  Gamma_star <- list()
  if (q_diff > 0) {
    for (l in seq_len(q_diff)) {
      Gamma_star[[l]] <- beta_hat[idx:(idx + k_star - 1), , drop = FALSE]   # k_star × k_dom
      idx <- idx + k_star
    }
  }
  
  # Global exogenous coefficients
  D_gamma_0 <- NULL
  D_gamma   <- list()
  if (k_global > 0) {
    D_gamma_0 <- beta_hat[idx:(idx + k_global - 1), , drop = FALSE]   # k_global × k_dom
    idx <- idx + k_global
    if (d_diff > 0) {
      for (l in seq_len(d_diff)) {
        D_gamma[[l]] <- beta_hat[idx:(idx + k_global - 1), , drop = FALSE]   # k_global × k_dom
        idx <- idx + k_global
      }
    }
  }
  
  # ---- Step 5: Convert to levels-VAR form ----
  levels_var <- vecm_to_levels_var(
    alpha       = alpha,
    beta        = beta,
    Gamma_dom   = Gamma_dom,
    Gamma_star  = Gamma_star,
    Gamma_star_0 = Gamma_star_0,
    k_dom       = k_dom,
    k_star      = k_star,
    p_lag       = p_lag,
    q_lag       = q_lag
  )
  
  # Compute residuals in levels space (approximately)
  # The VECM residuals are in first-difference space; we use them directly
  # as an approximation (they represent the innovation to x_it)
  residuals_vecm <- fit$residuals
  sigma_vecm     <- fit$sigma
  
  # Global exogenous D matrices in levels form
  D_levels <- list()
  if (k_global > 0) {
    # D_0 (contemporaneous) in levels ≈ D_gamma_0
    D_levels[[1]] <- D_gamma_0
    # For lagged D: analogous conversion
    if (d_lag > 0 && d_diff > 0) {
      for (l in seq_len(min(d_lag, d_diff))) {
        D_levels[[l + 1]] <- D_gamma[[l]]
      }
    }
  }
  
  return(list(
    # VECM-specific outputs
    alpha         = alpha,
    beta          = beta,
    Gamma_dom     = Gamma_dom,
    Gamma_star    = Gamma_star,
    Gamma_star_0  = Gamma_star_0,
    Pi            = alpha %*% t(beta),   # k_dom × m error-correction matrix
    rank          = rank,
    model_type    = "vecx",
    
    # Levels-VAR compatible outputs (for stack_to_gvar)
    intercept     = intercept_levels_final,
    trend         = trend_vecm,
    A             = levels_var$A,
    B             = levels_var$B,
    D             = D_levels,
    residuals     = residuals_vecm,
    sigma         = sigma_vecm,
    fitted        = fit$fitted,
    k_dom         = k_dom,
    k_star        = k_star,
    k_global      = k_global,
    p_lag         = p_lag,
    q_lag         = q_lag,
    d_lag         = d_lag,
    
    # Johansen object for reference
    ca_jo         = ca_obj
  ))
}


# ───────────────────────────────────────────────────────────────────────────────
# 2.  Convert VECM Parameters to Levels-VAR Form
# ───────────────────────────────────────────────────────────────────────────────

#' Convert VECM (error-correction) parameters back to a levels-VAR
#' representation so that the standard GVAR stacking machinery applies.
#'
#' The VECM:
#'   Δx_t = c + α β' z_{t-1} + Σ Γ^x_l Δx_{t-l} + Σ Γ^*_l Δx*_{t-l} + u_t
#'
#' Is equivalent to the levels VAR:
#'   x_t = c' + A_1 x_{t-1} + ... + A_p x_{t-p}
#'            + B_0 x*_t + B_1 x*_{t-1} + ... + B_q x*_{t-q} + u_t
#'
#' The conversion uses the standard VECM ↔ VAR relationship:
#'   Pi = α β' = A_1 + ... + A_p - I  (for domestic part)
#'   Gamma^x_l = -(A_{l+1} + ... + A_p)
#'
#' @param alpha        k_dom × rank loading matrix
#' @param beta         (k_dom + k_star) × rank cointegrating vectors
#' @param Gamma_dom    List of Γ^x_l matrices (k_dom × k_dom), l = 1,...,p-1
#' @param Gamma_star   List of Γ^*_l matrices (k_star × k_dom), l = 1,...,q-1
#' @param Gamma_star_0 k_star × k_dom contemporaneous Δx* coefficient
#' @param k_dom        Number of domestic variables
#' @param k_star       Number of star variables
#' @param p_lag        Domestic lag order in levels
#' @param q_lag        Foreign lag order in levels
#' @return List with A (list of A_1..A_p), B (list of B_0..B_q), intercept
vecm_to_levels_var <- function(alpha, beta, Gamma_dom, Gamma_star,
                               Gamma_star_0, k_dom, k_star,
                               p_lag, q_lag) {
  
  # Pi = alpha %*% t(beta)  is k_dom × (k_dom + k_star)
  Pi <- alpha %*% t(beta)
  
  # Partition Pi into domestic and star parts
  Pi_dom  <- Pi[, 1:k_dom, drop = FALSE]          # k_dom × k_dom
  Pi_star <- Pi[, (k_dom + 1):(k_dom + k_star), drop = FALSE]  # k_dom × k_star
  
  p_diff <- max(p_lag - 1, 0)
  q_diff <- max(q_lag - 1, 0)
  
  # ---- Domestic: recover A_1, ..., A_p from Pi_dom and Gamma_dom ----
  # Gamma^x_l = -(A_{l+1} + ... + A_p)   for l = 1, ..., p-1
  # Pi_dom    = (A_1 + ... + A_p) - I
  #
  # Working backwards:
  #   A_p = -Gamma^x_{p-1}  (if p > 1, else A_1 = Pi_dom + I)
  #   A_{p-1} = Gamma^x_{p-2} - Gamma^x_{p-1}  (if p > 2)
  #   ...
  #   A_1 = I + Pi_dom + Gamma^x_1  (if p > 1)
  #   A_1 = I + Pi_dom              (if p = 1)
  
  A <- list()
  
  if (p_lag == 0) {
    # No domestic lags – shouldn't normally happen but handle gracefully
    # A is empty
  } else if (p_lag == 1) {
    A[[1]] <- diag(k_dom) + Pi_dom
  } else {
    # p_lag >= 2
    # First, compute cumulative sums from the back
    # A_p = -Gamma^x_{p-1}
    if (p_diff > 0 && length(Gamma_dom) >= p_diff) {
      A[[p_lag]] <- -Gamma_dom[[p_diff]]
    } else {
      A[[p_lag]] <- matrix(0, k_dom, k_dom)
    }
    
    # A_l = Gamma^x_{l-1} - Gamma^x_l  for l = p, p-1, ..., 2
    if (p_lag >= 3) {
      for (l in (p_lag - 1):2) {
        Gamma_l   <- if (l <= length(Gamma_dom)) Gamma_dom[[l]] else matrix(0, k_dom, k_dom)
        Gamma_lm1 <- if ((l - 1) <= length(Gamma_dom)) Gamma_dom[[l - 1]] else matrix(0, k_dom, k_dom)
        A[[l]] <- Gamma_lm1 - Gamma_l
      }
    }
    
    # A_1 = I + Pi_dom + Gamma^x_1
    Gamma_1 <- if (length(Gamma_dom) >= 1) Gamma_dom[[1]] else matrix(0, k_dom, k_dom)
    A[[1]] <- diag(k_dom) + Pi_dom + Gamma_1
  }
  
  # ---- Star: recover B_0, ..., B_q from Pi_star and Gamma_star ----
  # The VECM star part:
  #   Δx_t includes: Pi_star * x*_{t-1} + Gamma_star_0 * Δx*_t
  #                  + Σ Gamma_star_l * Δx*_{t-l}
  #
  # In levels:
  #   B_0 = Gamma_star_0  (contemporaneous star in levels ≈ contemp diff coeff)
  #   The levels-VAR relationship for lagged stars follows analogously to domestic.
  #
  # More precisely, the levels VAR has:
  #   x_t = ... + B_0 x*_t + B_1 x*_{t-1} + ... + B_q x*_{t-q} + ...
  # The VECM has:
  #   Δx_t = ... + Pi_star x*_{t-1} + Gamma_star_0 Δx*_t + Σ Gamma_star_l Δx*_{t-l} + ...
  #
  # Substituting Δx*_t = x*_t - x*_{t-1} etc and matching:
  #   B_0' = Gamma_star_0'  (k_dom × k_star → we store k_star × k_dom)
  #   B_1' = Pi_star + Gamma_star_0' - Gamma_star_0' + Gamma_star_1'  etc.
  #
  # Simplified approach: match the VECM to levels form term-by-term
  
  B <- list()
  
  # B_0 (contemporaneous star): from Gamma_star_0
  B[[1]] <- Gamma_star_0   # k_star × k_dom
  
  if (q_lag == 0) {
    # No lagged star - just B_0
  } else if (q_lag == 1) {
    # B_1 comes from Pi_star (error-correction) minus the Δ contribution of B_0
    # In the VECM: Pi_star * x*_{t-1} - Gamma_star_0 * x*_{t-1} (from Δx*_t = x*_t - x*_{t-1})
    # So B_1' = Pi_star + Gamma_star_0' - Gamma_star_0' = Pi_star ... but actually:
    # Let's use the direct mapping:
    # Levels: B_0 x*_t + B_1 x*_{t-1}
    # VECM:  Gamma_star_0 (x*_t - x*_{t-1}) + Pi_star x*_{t-1}
    #      = Gamma_star_0 x*_t + (Pi_star - Gamma_star_0) x*_{t-1}
    # So B_0 = Gamma_star_0 and B_1 = Pi_star - Gamma_star_0 (in model convention)
    # Pi_star is k_dom × k_star (model), Gamma_star_0 is k_star × k_dom (OLS).
    # Subtract in model convention (k_dom × k_star), then transpose for storage.
    B[[2]] <- t(Pi_star - t(Gamma_star_0))   # k_star × k_dom
  } else {
    # q_lag >= 2: general case
    # VECM: Gamma_star_0 Δx*_t + Gamma_star_1 Δx*_{t-1} + ... + Pi_star x*_{t-1}
    # Expanding Δx*_{t-l} = x*_{t-l} - x*_{t-l-1}:
    # Coefficient on x*_t:         Gamma_star_0                    → B_0 = Gamma_star_0
    # Coefficient on x*_{t-1}:     Pi_star - Gamma_star_0 + Gamma_star_1  → B_1
    # Coefficient on x*_{t-l}:     Gamma_star_l - Gamma_star_{l-1}  for l = 2,...,q-1  → B_l
    # Coefficient on x*_{t-q}:     -Gamma_star_{q-1}                → B_q
    
    # B_1 (transposed to k_star × k_dom format)
    Gs1 <- if (length(Gamma_star) >= 1) Gamma_star[[1]] else matrix(0, k_star, k_dom)
    B[[2]] <- t(Pi_star) - Gamma_star_0 + Gs1   # k_star × k_dom
    
    # B_l for l = 2, ..., q-1
    if (q_lag >= 3) {
      for (l in 2:(q_lag - 1)) {
        Gs_l   <- if (l <= length(Gamma_star)) Gamma_star[[l]] else matrix(0, k_star, k_dom)
        Gs_lm1 <- if ((l - 1) <= length(Gamma_star)) Gamma_star[[l - 1]] else matrix(0, k_star, k_dom)
        B[[l + 1]] <- Gs_l - Gs_lm1   # k_star × k_dom
      }
    }
    
    # B_q
    Gs_last <- if (length(Gamma_star) >= q_diff) Gamma_star[[q_diff]] else matrix(0, k_star, k_dom)
    B[[q_lag + 1]] <- -Gs_last   # k_star × k_dom
  }
  
  # Intercept: the VECM intercept maps approximately to the levels intercept
  # In a VECM with restricted constant (ecdet = "const"), the constant is
  # inside the cointegrating relation, so the VECM intercept is approximately
  # the levels intercept. For simplicity, pass through.
  intercept_levels <- rep(0, k_dom)   # Will be overwritten from VECM intercept
  
  return(list(
    A         = A,
    B         = B,
    intercept = intercept_levels
  ))
}


# ───────────────────────────────────────────────────────────────────────────────
# 3.  Master GVECM Estimation Function
# ───────────────────────────────────────────────────────────────────────────────

#' Estimate the full Global VECM: individual VECX* + global stacking.
#'
#' For each unit with cointegrating rank r_i > 0, estimates a VECX* model
#' and converts it to levels-VAR form.  Units with r_i = 0 are estimated
#' as standard VARX* models (in first differences).  The levels-VAR
#' representations are then stacked via stack_to_gvar() to produce a
#' standard gvar_model object.
#'
#' @param gvar_data   Output from prepare_gvar_dataset()
#' @param rank_vec    Named integer vector of cointegrating ranks per unit
#'                    (from panel_cointegration())
#' @param ecdet       Deterministic specification for Johansen: "const",
#'                    "trend", or "none"
#' @return            A gvar_model object with additional $vecm_info
estimate_gvecm <- function(gvar_data, rank_vec, ecdet = "const",
                           deterministic = "intercept") {
  
  print_banner("GVECM Model Estimation (Error-Correction Form)")
  
  unit_names <- gvar_data$unit_names
  N <- length(unit_names)
  
  # We need the original (pre-lagged) data for Johansen estimation
  # Reconstruct from gvar_data
  # star_list is stored in gvar_data
  star_list <- gvar_data$star_list
  
  unit_fits  <- setNames(vector("list", N), unit_names)
  alpha_list <- setNames(vector("list", N), unit_names)
  beta_list  <- setNames(vector("list", N), unit_names)
  
  for (u in unit_names) {
    cd   <- gvar_data$unit_data[[u]]
    r_i  <- rank_vec[u]
    
    # Determine global exogenous for this unit
    g_exog  <- NULL
    g_d_lag <- NULL
    if (!is.null(gvar_data$global_config) &&
        u != gvar_data$global_config$dominant_unit) {
      g_exog  <- gvar_data$global_config$global_series
      g_d_lag <- gvar_data$global_config$d_lag
    }
    
    # Cap rank at system dimension - 1 (full rank = no cointegration)
    m_i <- cd$k_dom + ncol(star_list[[u]])
    if (r_i >= m_i) {
      message(sprintf("  [%s] rank (%d) >= system dim (%d) – treating as VARX* (no cointegration).",
                      u, r_i, m_i))
      r_i <- 0
    }
    
    if (r_i > 0) {
      message(sprintf("  Estimating VECX*(%d,%d) with rank %d for unit '%s' ...",
                      cd$p_lag, cd$q_lag, r_i, u))
      
      # Need the original domestic and star matrices (not lagged/trimmed)
      # Reconstruct from Y (which has been trimmed) – we need the full series.
      # The star_list has full T observations.
      # For domestic, we use the Y column names to find the original data.
      # WORKAROUND: use the data stored in gvar_data.
      # We reconstruct domestic from cd$Y plus the max_lag leading observations.
      # Actually, the star_list is stored in gvar_data, and the original
      # data_list can be reconstructed.  For now, use cd$Y dimensions
      # to guide reconstruction.  We need the full T data.
      
      # The cleanest approach: the caller should pass the original data_list.
      # Since gvar_data stores unit_data (trimmed) and star_list (full T),
      # we can infer the original domestic data from the star_list T.
      T_full <- nrow(star_list[[u]])
      T_trim <- nrow(cd$Y)
      max_lag <- T_full - T_trim
      
      # Reconstruct domestic: we need all T_full rows.
      # cd$Y only has T_trim rows.  We can reconstruct from cd$X
      # but that's fragile.  Better: store full data in gvar_data.
      # FALLBACK: use the Y and lagged X to reconstruct.
      # Actually, for Johansen we need the LEVELS data, not the regression format.
      # Let's reconstruct from gvar_data.
      
      # We stored the full domestic data implicitly:
      # cd$Y has the dependent variable (post-lag observations)
      # The lagged domestic in cd$X gives us access to prior observations.
      # For Johansen, we need the full T_full × k_dom matrix.
      
      # Extract from X: first p_lag blocks of k_dom columns are lags
      # Together with Y, we can reconstruct the full series.
      
      # Simple reconstruction:
      # Y[t] = domestic[max_lag + t] for t = 1, ..., T_trim
      # X includes lagged domestic: domestic[max_lag + t - l] for l = 1, ..., p_lag
      # So we can recover domestic[1:T_full].
      
      # But the safest approach: re-read from star_list T dimension
      # and note that estimate_unit_vecx internally needs only the
      # domestic and star matrices.  We should pass the FULL data.
      
      # Since the full data was aligned and stored in star_list,
      # we know T_full.  For the domestic variables, we need to
      # recover them.  Let's build them from Y and the lag structure.
      
      # Actually, let me just store the reconstructed domestic:
      dom_full <- matrix(NA, nrow = T_full, ncol = cd$k_dom)
      colnames(dom_full) <- colnames(cd$Y)
      # Fill from Y (last T_trim rows)
      dom_full[(max_lag + 1):T_full, ] <- cd$Y
      # Fill from lagged X (if available)
      if (cd$p_lag > 0) {
        lag1_cols <- grep("_lag1$", colnames(cd$X))
        if (length(lag1_cols) == cd$k_dom) {
          dom_full[max_lag:(T_full - 1), ] <- cd$X[, lag1_cols]
          # Work backwards for more lags
          for (l in 2:min(cd$p_lag, max_lag)) {
            lag_cols <- grep(paste0("_lag", l, "$"), colnames(cd$X))
            if (length(lag_cols) == cd$k_dom) {
              dom_full[(max_lag + 1 - l):(T_full - l), ] <- cd$X[, lag_cols]
            }
          }
        }
      }
      # Fill any remaining NAs at the start from the contemporaneous star pattern
      # For safety, if NAs remain, fill with the earliest available value
      for (col in seq_len(ncol(dom_full))) {
        first_valid <- which(!is.na(dom_full[, col]))[1]
        if (!is.na(first_valid) && first_valid > 1) {
          dom_full[1:(first_valid - 1), col] <- dom_full[first_valid, col]
        }
      }
      
      star_full <- star_list[[u]]
      
      fit <- estimate_unit_vecx(
        domestic    = dom_full,
        star        = star_full,
        rank        = r_i,
        p_lag       = cd$p_lag,
        q_lag       = cd$q_lag,
        ecdet       = ecdet,
        global_exog   = g_exog,
        d_lag         = g_d_lag,
        deterministic = deterministic
      )
      
      alpha_list[[u]] <- fit$alpha
      beta_list[[u]]  <- fit$beta
      
    } else {
      message(sprintf("  No cointegration for unit '%s' – estimating VARX*(%d,%d) ...",
                      u, cd$p_lag, cd$q_lag))
      
      # Estimate standard VARX* (same as estimate_gvar)
      k_dom    <- cd$k_dom
      k_star   <- cd$k_star
      k_global <- if (!is.null(cd$k_global)) cd$k_global else 0
      
      det <- build_deterministic_columns(nrow(cd$Y), deterministic)
      X_full <- if (!is.null(det$D_mat)) cbind(det$D_mat, cd$X) else cd$X
      fit_ols <- ols_estimate(cd$Y, X_full, intercept = FALSE)
      
      beta_raw <- fit_ols$beta
      dp <- partition_deterministic(beta_raw, det$n_det,
                                    det$has_intercept, det$has_trend, k_dom)
      intercept <- dp$intercept
      trend     <- dp$trend
      idx       <- dp$start_idx
      
      B <- list()
      B[[1]] <- beta_raw[idx:(idx + k_star - 1), , drop = FALSE]
      idx <- idx + k_star
      
      D <- list()
      if (k_global > 0) {
        D[[1]] <- beta_raw[idx:(idx + k_global - 1), , drop = FALSE]
        idx <- idx + k_global
      }
      
      A <- list()
      for (l in seq_len(cd$p_lag)) {
        A[[l]] <- beta_raw[idx:(idx + k_dom - 1), , drop = FALSE]
        idx <- idx + k_dom
      }
      
      for (l in seq_len(cd$q_lag)) {
        B[[l + 1]] <- beta_raw[idx:(idx + k_star - 1), , drop = FALSE]
        idx <- idx + k_star
      }
      
      if (k_global > 0) {
        d_lag_i <- if (!is.null(cd$d_lag)) cd$d_lag else 0
        for (l in seq_len(d_lag_i)) {
          D[[l + 1]] <- beta_raw[idx:(idx + k_global - 1), , drop = FALSE]
          idx <- idx + k_global
        }
      }
      
      fit <- list(
        intercept  = intercept,
        trend      = trend,
        A          = A,
        B          = B,
        D          = D,
        residuals  = fit_ols$residuals,
        sigma      = fit_ols$sigma,
        fitted     = fit_ols$fitted,
        k_dom      = k_dom,
        k_star     = k_star,
        k_global   = k_global,
        p_lag      = cd$p_lag,
        q_lag      = cd$q_lag,
        d_lag      = if (!is.null(cd$d_lag)) cd$d_lag else 0,
        model_type = "varx",
        rank       = 0
      )
      
      alpha_list[[u]] <- NULL
      beta_list[[u]]  <- NULL
    }
    
    unit_fits[[u]] <- fit
  }
  
  # Stage 2: stack into global VAR using levels-form coefficients
  message("\n  Stacking VECX*/VARX* models into the global system ...")
  
  # Build data_list for stacking (column names needed for link matrices)
  data_list_for_stack <- list()
  for (u in unit_names) {
    cd <- gvar_data$unit_data[[u]]
    data_list_for_stack[[u]] <- cd$Y
  }
  
  global <- stack_to_gvar(unit_fits, data_list_for_stack, gvar_data$W,
                          global_config = gvar_data$global_config)
  
  # Merge into output
  result <- c(global, list(
    unit_fits  = unit_fits,
    unit_names = unit_names,
    freq       = gvar_data$freq,
    W          = gvar_data$W,
    vecm_info  = list(
      rank_vec   = rank_vec,
      alpha      = alpha_list,
      beta       = beta_list,
      model_type = "gvecm"
    )
  ))
  
  class(result) <- "gvar_model"
  message("[GVECM] Estimation complete.")
  return(result)
}
