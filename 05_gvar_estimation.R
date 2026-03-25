###############################################################################
#  05_gvar_estimation.R  –  Core GVAR Model Estimation
#
#  The GVAR estimation proceeds in two stages:
#
#  Stage 1 (Country-level):
#    Estimate individual VARX*(p_i, q_i) models for each unit i by OLS.
#
#  Stage 2 (Global stacking):
#    Use the link/weight matrix W to stack the individual models into a
#    single global VAR of the form:
#
#        G_0 * x_t  =  a_0  +  G_1 * x_{t-1}  +  ...  +  G_p * x_{t-p}  +  u_t
#
#    Then solve for the reduced form:
#
#        x_t  =  F_0  +  F_1 x_{t-1}  +  ...  +  F_p x_{t-p}  +  e_t
#
#  Contents
#  --------
#  1. estimate_unit_varx()   – OLS estimation of one VARX* model
#  2. stack_to_gvar()        – combine individual models into the global VAR
#  3. estimate_gvar()        – master function: individual + stacking
###############################################################################


# ───────────────────────────────────────────────────────────────────────────────
# 1.  Estimate a Single VARX* Model
# ───────────────────────────────────────────────────────────────────────────────

#' Estimate the VARX*(p, q) model for a single unit by OLS.
#'
#' Model:  x_{it} = c_i + A_{i1} x_{i,t-1} + ... + A_{ip} x_{i,t-p}
#'                      + B_{i0} x*_{it}    + ... + B_{iq} x*_{i,t-q} + u_{it}
#'
#' After OLS estimation we partition the coefficient matrix into blocks
#' corresponding to the intercept, A matrices, and B matrices.
#'
#' @param domestic   T × k_i  matrix
#' @param star       T × k*   matrix
#' @param p_lag      Integer; domestic lag order
#' @param q_lag      Integer; foreign lag order
#' @return  A list with:
#'   \item{beta}{Full coefficient matrix from OLS}
#'   \item{intercept}{k_i-vector of intercepts}
#'   \item{A}{List of A_1,...,A_p coefficient matrices (k_i × k_i)}
#'   \item{B}{List of B_0,...,B_q coefficient matrices (k_i × k_star)}
#'   \item{residuals}{Residual matrix}
#'   \item{sigma}{Residual covariance}
#'   \item{fitted}{Fitted values}
#'   \item{k_dom, k_star, p_lag, q_lag}{Metadata}
estimate_unit_varx <- function(domestic, star, p_lag, q_lag) {

  cd  <- prepare_country_data(domestic, star, p_lag, q_lag)
  fit <- ols_estimate(cd$Y, cd$X, intercept = TRUE)

  k_dom  <- cd$k_dom
  k_star <- cd$k_star

  # ---- Partition beta (rows) into blocks ----
  # Row layout of beta (with intercept prepended):
  #   Row 1           : intercept
  #   Next k_star     : B_0  (contemporaneous star)
  #   Next k_dom      : A_1
  #   Next k_dom      : A_2  ...  (p_lag blocks)
  #   Next k_star     : B_1
  #   Next k_star     : B_2  ...  (q_lag blocks)

  beta <- fit$beta   # (1 + k_star + p*k_dom + q*k_star) × k_dom
  idx <- 1

  # Intercept
  intercept <- beta[idx, ]
  idx <- idx + 1

  # B_0 (contemporaneous star)
  B <- list()
  B[[1]] <- beta[idx:(idx + k_star - 1), , drop = FALSE]
  idx <- idx + k_star

  # A_1, ..., A_p
  A <- list()
  for (l in seq_len(p_lag)) {
    A[[l]] <- beta[idx:(idx + k_dom - 1), , drop = FALSE]
    idx <- idx + k_dom
  }

  # B_1, ..., B_q
  for (l in seq_len(q_lag)) {
    B[[l + 1]] <- beta[idx:(idx + k_star - 1), , drop = FALSE]
    idx <- idx + k_star
  }

  return(list(
    beta      = beta,
    intercept = intercept,
    A         = A,
    B         = B,
    residuals = fit$residuals,
    sigma     = fit$sigma,
    fitted    = fit$fitted,
    k_dom     = k_dom,
    k_star    = k_star,
    p_lag     = p_lag,
    q_lag     = q_lag
  ))
}


# ───────────────────────────────────────────────────────────────────────────────
# 2.  Stack Individual Models into the Global VAR
# ───────────────────────────────────────────────────────────────────────────────

#' Combine all country VARX* models into a global VAR.
#'
#' The key identity is:
#'   For unit i define z_{it} = (x_{it}, x*_{it})'
#'   The VARX* can be re-written as:
#'     A_{i0} z_{it} = a_{i0} + A_{i1} z_{i,t-1} + ... + u_{it}
#'
#'   where A_{i0} = [I_ki , -B_{i0}]   (left-hand side matrix).
#'
#'   Then z_{it} = W_i * x_t  where W_i is a "link matrix" derived from the
#'   trade weights.  Stacking all units:
#'
#'     G_0 x_t = a_0 + G_1 x_{t-1} + ... + G_p x_{t-p} + u_t
#'
#'   The reduced form is:  x_t = G_0^{-1} a_0 + G_0^{-1} G_1 x_{t-1} + ...
#'
#' @param unit_fits   Named list of fitted VARX* models (from estimate_unit_varx)
#' @param data_list   Named list of T × k_i matrices
#' @param W           N × N weight matrix
#' @return  A list:
#'   \item{F}{List of reduced-form coefficient matrices F_1,...,F_p}
#'   \item{F0}{Reduced-form intercept vector}
#'   \item{G0, G}{The structural-form matrices}
#'   \item{Sigma}{Global residual covariance}
#'   \item{companion}{Companion-form matrix}
#'   \item{k_total}{Total number of endogenous variables across all units}
#'   \item{p_global}{Global lag order (max across units)}
#'   \item{var_names}{Character vector of all variable names in stacking order}
stack_to_gvar <- function(unit_fits, data_list, W, global_config = NULL) {

  unit_names <- names(unit_fits)
  N <- length(unit_names)

  # --- Determine dimensions ---
  k_vec <- sapply(unit_fits, function(f) f$k_dom)     # domestic vars per unit
  k_total <- sum(k_vec)                                # total endogenous vars
  p_global <- max(sapply(unit_fits, function(f) f$p_lag))
  q_global <- max(sapply(unit_fits, function(f) f$q_lag))
  s_global <- max(p_global, q_global)                  # effective global lag

  # Identify the common variables (used for star construction)
  # Exclude global variable names from common_vars (they are handled separately)
  common_vars <- Reduce(intersect, lapply(data_list, colnames))
  if (!is.null(global_config)) {
    common_vars <- setdiff(common_vars, global_config$var_names)
  }

  # Also account for global exogenous lags in the effective global lag
  d_global <- 0
  if (!is.null(global_config)) {
    d_global <- global_config$d_lag
    s_global <- max(s_global, d_global)
  }

  # --- Build link matrices W_i for each unit ---
  # W_i maps the global vector x_t  (k_total × 1)  to z_{it} = (x_{it}, x*_{it}, d_{it})'
  #
  # For unit i:
  #   rows 1..k_i           select x_{it}  from x_t  (identity)
  #   rows k_i+1..k_i+k*    select x*_{it} = W_{ij} * x_{jt} (weighted)
  #   rows k_i+k*+1..end    select d_{it} (global vars from dominant unit)

  # Build an ordered list of ALL variable names in the global vector
  var_names <- c()
  for (u in unit_names) {
    var_names <- c(var_names, paste0(u, ".", colnames(data_list[[u]])))
  }

  # Create a mapping from (unit, var_name) to column index in x_t
  global_idx <- setNames(seq_along(var_names), var_names)

  # --- Construct G_0, G_1, ..., G_s matrices (k_total × k_total) ---
  G0 <- matrix(0, nrow = k_total, ncol = k_total)
  G  <- replicate(s_global, matrix(0, nrow = k_total, ncol = k_total),
                  simplify = FALSE)
  a0_vec     <- numeric(k_total)
  a_trend_vec <- numeric(k_total)

  # Also collect residuals (stack across units)
  T_eff_vec <- sapply(unit_fits, function(f) nrow(f$residuals))
  T_eff <- min(T_eff_vec)

  resid_global <- matrix(0, nrow = T_eff, ncol = k_total)

  row_offset <- 0   # tracks current row block in the stacked system

  for (u in unit_names) {
    fit <- unit_fits[[u]]
    k_i <- fit$k_dom
    k_s <- fit$k_star
    k_g <- if (!is.null(fit$k_global)) fit$k_global else 0
    rows_i <- (row_offset + 1):(row_offset + k_i)

    # ---- Build the link matrix W_i  (k_i + k_s + k_g) × k_total ----
    Wi <- matrix(0, nrow = k_i + k_s + k_g, ncol = k_total)

    # Top block: identity for domestic variables of unit u
    dom_vars <- colnames(data_list[[u]])
    for (j in seq_len(k_i)) {
      gname <- paste0(u, ".", dom_vars[j])
      Wi[j, global_idx[gname]] <- 1
    }

    # Middle block: weighted selection for star variables
    for (v in seq_along(common_vars)) {
      for (j in unit_names) {
        gname <- paste0(j, ".", common_vars[v])
        if (gname %in% names(global_idx)) {
          Wi[k_i + v, global_idx[gname]] <- W[u, j]
        }
      }
    }

    # Bottom block: global exogenous variables (select from dominant unit)
    if (k_g > 0 && !is.null(global_config)) {
      for (g in seq_len(k_g)) {
        gname <- paste0(global_config$dominant_unit, ".",
                        global_config$var_names[g])
        if (gname %in% names(global_idx)) {
          Wi[k_i + k_s + g, global_idx[gname]] <- 1
        }
      }
    }

    # ---- A_{i,0}  =  [ I_{k_i} , -B_{i,0} , -D_{i,0} ]  ----
    Ai0 <- cbind(diag(k_i), -t(fit$B[[1]]))
    if (k_g > 0 && length(fit$D) > 0) {
      Ai0 <- cbind(Ai0, -t(fit$D[[1]]))
    }

    # G_0 contribution:  A_{i,0} * W_i
    G0[rows_i, ] <- Ai0 %*% Wi

    # Intercept and trend
    a0_vec[rows_i] <- fit$intercept
    if (!is.null(fit$trend)) {
      a_trend_vec[rows_i] <- fit$trend
    }

    # ---- Build A_{i,l} matrices for l = 1,...,s_global ----
    for (l in seq_len(s_global)) {
      # Domestic part: A_{i,l} if l <= p_i, else 0
      A_il <- if (l <= fit$p_lag) fit$A[[l]] else matrix(0, k_i, k_i)

      # Foreign part: B_{i,l} if l <= q_i, else 0
      B_il <- if (l <= fit$q_lag && (l + 1) <= length(fit$B)) {
        fit$B[[l + 1]]
      } else {
        matrix(0, k_s, k_i)
      }

      # Combined: [A_il , B_il']
      AB_il <- cbind(A_il, t(B_il))

      # Global exogenous part: D_{i,l} if l <= d_i, else 0
      if (k_g > 0) {
        d_lag_i <- if (!is.null(fit$d_lag)) fit$d_lag else 0
        D_il <- if (l <= d_lag_i && (l + 1) <= length(fit$D)) {
          fit$D[[l + 1]]
        } else {
          matrix(0, k_g, k_i)
        }
        AB_il <- cbind(AB_il, t(D_il))
      }

      # G_l contribution:  AB_il * W_i
      G[[l]][rows_i, ] <- AB_il %*% Wi
    }

    # Stack residuals (take last T_eff rows for conformability)
    n_res <- nrow(fit$residuals)
    resid_global[, rows_i] <- fit$residuals[(n_res - T_eff + 1):n_res, ]

    row_offset <- row_offset + k_i
  }

  # ---- Solve for reduced form ----
  # x_t = G_0^{-1} a_0  +  G_0^{-1} G_1 x_{t-1}  + ...
  G0_inv <- solve(G0)

  F0       <- as.numeric(G0_inv %*% a0_vec)
  F0_trend <- as.numeric(G0_inv %*% a_trend_vec)
  FF <- lapply(G, function(Gl) G0_inv %*% Gl)

  # Global residual covariance
  Sigma_global <- crossprod(resid_global) / T_eff

  # ---- Companion matrix ----
  comp <- make_companion(FF)

  # Check stability (all eigenvalues inside unit circle)
  eig_mod <- Mod(eigen(comp, only.values = TRUE)$values)
  is_stable <- all(eig_mod < 1)
  if (!is_stable) {
    message("[GVAR] WARNING: The global model is NOT stable (max eigenvalue modulus = ",
            round(max(eig_mod), 4), ").")
  } else {
    message("[GVAR] Global model is stable (max eigenvalue modulus = ",
            round(max(eig_mod), 4), ").")
  }

  return(list(
    F0        = F0,
    F0_trend  = F0_trend,
    F         = FF,
    G0        = G0,
    G         = G,
    G0_inv    = G0_inv,
    Sigma     = Sigma_global,
    companion = comp,
    k_total   = k_total,
    p_global  = s_global,
    var_names = var_names,
    is_stable = is_stable,
    eig_mod   = eig_mod,
    residuals = resid_global
  ))
}


# ───────────────────────────────────────────────────────────────────────────────
# 3.  Master Estimation Function
# ───────────────────────────────────────────────────────────────────────────────

#' Estimate the full GVAR model: individual VARX* + global stacking.
#'
#' @param gvar_data  Output from prepare_gvar_dataset()
#' @return           A list containing everything from stack_to_gvar() plus
#'                   the individual unit fits.
estimate_gvar <- function(gvar_data, deterministic = "intercept") {

  print_banner("GVAR Model Estimation")

  unit_names <- gvar_data$unit_names
  N <- length(unit_names)

  # Stage 1: estimate individual VARX* models
  unit_fits <- setNames(vector("list", N), unit_names)

  for (u in unit_names) {
    cd <- gvar_data$unit_data[[u]]
    message(sprintf("  Estimating VARX*(%d,%d) for unit '%s' ...",
                    cd$p_lag, cd$q_lag, u))

    k_dom    <- cd$k_dom
    k_star   <- cd$k_star
    k_global <- if (!is.null(cd$k_global)) cd$k_global else 0
    p_lag    <- cd$p_lag
    q_lag    <- cd$q_lag
    d_lag    <- if (!is.null(cd$d_lag)) cd$d_lag else 0

    det_info <- build_deterministic_columns(nrow(cd$Y), deterministic)
    X_full <- if (!is.null(det_info$D_mat)) cbind(det_info$D_mat, cd$X) else cd$X
    fit <- ols_estimate(cd$Y, X_full, intercept = FALSE)

    beta <- fit$beta
    dp <- partition_deterministic(beta, det_info$n_det,
                                   det_info$has_intercept,
                                   det_info$has_trend, k_dom)
    intercept <- dp$intercept
    trend     <- dp$trend
    idx       <- dp$start_idx

    # B_0 (contemporaneous star)
    B <- list()
    B[[1]] <- beta[idx:(idx + k_star - 1), , drop = FALSE]
    idx <- idx + k_star

    # D_0 (contemporaneous global exogenous) – if present
    D <- list()
    if (k_global > 0) {
      D[[1]] <- beta[idx:(idx + k_global - 1), , drop = FALSE]
      idx <- idx + k_global
    }

    # A_1,...,A_p
    A <- list()
    for (l in seq_len(p_lag)) {
      A[[l]] <- beta[idx:(idx + k_dom - 1), , drop = FALSE]
      idx <- idx + k_dom
    }

    # B_1,...,B_q
    for (l in seq_len(q_lag)) {
      B[[l + 1]] <- beta[idx:(idx + k_star - 1), , drop = FALSE]
      idx <- idx + k_star
    }

    # D_1,...,D_d (lagged global exogenous)
    if (k_global > 0 && d_lag > 0) {
      for (l in seq_len(d_lag)) {
        D[[l + 1]] <- beta[idx:(idx + k_global - 1), , drop = FALSE]
        idx <- idx + k_global
      }
    }

    unit_fits[[u]] <- list(
      beta      = beta,
      intercept = intercept,
      trend     = trend,
      A         = A,
      B         = B,
      D         = D,
      residuals = fit$residuals,
      sigma     = fit$sigma,
      fitted    = fit$fitted,
      k_dom     = k_dom,
      k_star    = k_star,
      k_global  = k_global,
      p_lag     = p_lag,
      q_lag     = q_lag,
      d_lag     = d_lag
    )
  }

  # Stage 2: stack into global VAR
  message("\n  Stacking individual models into the global VAR ...")

  # We need the original (not lagged) data list to build the link matrices
  # Reconstruct from gvar_data: the original matrices are implicitly inside
  # unit_data, but we need the full T versions.  The caller must supply them.
  # WORKAROUND: we re-derive from the unit_data and star_list.

  # Actually, we need the raw data_list.  We'll add it to gvar_data in
  # prepare_gvar_dataset(). For now, let's build a simple version from
  # the Y matrices (not ideal but functional).

  # Better approach: store data_list inside gvar_data.  We assume the
  # caller has done so (see updated prepare_gvar_dataset).

  # We need to pass a data_list with column names for link-matrix construction.
  # Build it from the gvar_data structure.
  data_list_for_stack <- list()
  for (u in unit_names) {
    cd <- gvar_data$unit_data[[u]]
    # The Y matrix has the correct column names for domestic variables
    data_list_for_stack[[u]] <- cd$Y
  }

  global <- stack_to_gvar(unit_fits, data_list_for_stack, gvar_data$W,
                          global_config = gvar_data$global_config)

  # Merge everything into one output object
  result <- c(global, list(
    unit_fits     = unit_fits,
    unit_names    = unit_names,
    freq          = gvar_data$freq,
    W             = gvar_data$W,
    deterministic = deterministic
  ))

  class(result) <- "gvar_model"
  message("[GVAR] Estimation complete.")
  return(result)
}
