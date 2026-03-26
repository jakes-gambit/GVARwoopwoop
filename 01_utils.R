###############################################################################
#  01_utils.R  –  Utility / Helper Functions for the GVAR Toolkit
#
#  Contents
#  --------
#  1. check_and_install()   – silently install missing CRAN packages
#  2. load_packages()       – load all required libraries
#  3. make_companion()      – build VAR companion matrix from coefficient list
#  4. vec_to_mat()          – reshape a stacked vector into a T×k matrix
#  5. lag_matrix()          – create a matrix of lagged values
#  6. ols_estimate()        – simple OLS with optional intercept
#  7. information_criteria() – compute AIC / BIC for a given residual matrix
#  8. print_banner()        – pretty-print section headers to the console
#  9. build_deterministic_columns() – create intercept / trend columns
# 10. partition_deterministic()     – extract deterministic coefficients from beta
###############################################################################

# ───────────────────────────────────────────────────────────────────────────────
# 1.  Package Management
# ───────────────────────────────────────────────────────────────────────────────

#' Silently install any packages that are not yet available
#' @param pkgs Character vector of CRAN package names
check_and_install <- function(pkgs) {
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    message("[GVAR] Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}

#' Load all packages needed by the toolkit.
#' Installs anything that is missing first.
load_packages <- function() {
  required <- c(
    "vars",        # VAR estimation, IRF
    "urca",        # Unit-root / cointegration tests
    "tseries",     # Additional time-series tests (ADF, KPSS)
    "MASS",        # ginv (generalised inverse)
    "Matrix",      # Sparse / efficient matrix operations
    "ggplot2",     # Plotting
    "reshape2",    # Data reshaping for ggplot
    "stats"        # Base stats (Kalman filter helpers)
  )
  check_and_install(required)
  invisible(lapply(required, library, character.only = TRUE, quietly = TRUE))
  message("[GVAR] All packages loaded successfully.")
}


# ───────────────────────────────────────────────────────────────────────────────
# 2.  Companion-Form Matrix
# ───────────────────────────────────────────────────────────────────────────────

#' Build the companion matrix for a VAR(p) or VARX(p,q) system.
#'
#' Given a list of k×k coefficient matrices  A1, A2, ..., Ap  the companion
#' matrix is the (k*p) × (k*p) block matrix
#'
#'   | A1  A2  ... Ap |
#'   | I   0   ... 0  |
#'   | 0   I   ... 0  |
#'   | :           :  |
#'   | 0   0   I   0  |
#'
#' @param coef_list  List of coefficient matrices (each k×k), in order A1..Ap
#' @return           Companion matrix of dimension (k*p) × (k*p)
make_companion <- function(coef_list) {
  p <- length(coef_list)
  k <- nrow(coef_list[[1]])
  kp <- k * p

  comp <- matrix(0, nrow = kp, ncol = kp)

  # First block-row: [A1 | A2 | ... | Ap]
  for (i in seq_len(p)) {
    cols <- ((i - 1) * k + 1):(i * k)
    comp[1:k, cols] <- coef_list[[i]]
  }

  # Identity sub-diagonal blocks
  if (p > 1) {
    comp[(k + 1):kp, 1:((p - 1) * k)] <- diag((p - 1) * k)
  }

  return(comp)
}


# ───────────────────────────────────────────────────────────────────────────────
# 3.  Lag Matrix Constructor
# ───────────────────────────────────────────────────────────────────────────────

#' Create a matrix of lagged values for each column of Y.
#'
#' @param Y      T × k matrix (or data frame) of time-series data
#' @param max_lag  Integer; maximum number of lags to include
#' @return       A (T - max_lag) × (k * max_lag) matrix, with columns ordered
#'               as  Y1_lag1, Y2_lag1, ..., Yk_lag1, Y1_lag2, ...
lag_matrix <- function(Y, max_lag) {
  Y <- as.matrix(Y)
  n <- nrow(Y)
  k <- ncol(Y)

  if (max_lag >= n) stop("max_lag must be less than the number of observations.")

  lags <- matrix(NA, nrow = n - max_lag, ncol = k * max_lag)
  col_names <- c()
  idx <- 1
  for (lag in seq_len(max_lag)) {
    rows <- (max_lag - lag + 1):(n - lag)
    lags[, idx:(idx + k - 1)] <- Y[rows, , drop = FALSE]
    col_names <- c(col_names, paste0(colnames(Y), "_lag", lag))
    idx <- idx + k
  }
  colnames(lags) <- col_names
  return(lags)
}


# ───────────────────────────────────────────────────────────────────────────────
# 4.  Simple OLS Estimator
# ───────────────────────────────────────────────────────────────────────────────

#' OLS regression  Y = X %*% beta + e   with optional intercept.
#'
#' @param Y          T × k  matrix of dependent variables
#' @param X          T × m  matrix of regressors (already lagged / stacked)
#' @param intercept  Logical; if TRUE a column of ones is prepended to X
#' @return A list with elements:
#'   \item{beta}{Coefficient matrix (m+1) × k  (or m × k if no intercept)}
#'   \item{residuals}{T × k residual matrix}
#'   \item{fitted}{T × k fitted values}
#'   \item{sigma}{k × k residual covariance matrix (1/T scaling)}
ols_estimate <- function(Y, X, intercept = TRUE, lambda = 0) {
  Y <- as.matrix(Y)
  X <- as.matrix(X)

  if (intercept) {
    X <- cbind(1, X)
  }

  p <- ncol(X)

  # Ridge regression: beta = (X'X + lambda * D)^{-1} X'Y
  # D is the penalty matrix: identity but with a zero for the intercept column
  # so the intercept is never penalised.
  if (lambda > 0) {
    D <- diag(p)
    if (intercept) D[1, 1] <- 0   # do not penalise the intercept
    beta <- solve(crossprod(X) + lambda * D, crossprod(X, Y))
  } else {
    beta <- solve(crossprod(X), crossprod(X, Y))
  }

  fitted    <- X %*% beta
  residuals <- Y - fitted

  TT <- nrow(Y)
  sigma <- crossprod(residuals) / TT   # MLE-type (divide by T)

  return(list(
    beta      = beta,
    residuals = residuals,
    fitted    = fitted,
    sigma     = sigma
  ))
}


# ───────────────────────────────────────────────────────────────────────────────
# 5.  Information Criteria (AIC / BIC)
# ───────────────────────────────────────────────────────────────────────────────

#' Compute AIC and BIC for a VAR-type model given residuals.
#'
#' AIC = log|Sigma| + 2 * n_params / T
#' BIC = log|Sigma| + log(T) * n_params / T
#'
#' where Sigma is the MLE covariance of residuals and n_params is the total
#' number of freely estimated parameters across all equations.
#'
#' @param residuals  T × k residual matrix
#' @param n_params   Integer; total number of estimated parameters
#' @return Named list with elements `aic` and `bic`
information_criteria <- function(residuals, n_params) {
  TT <- nrow(residuals)
  k  <- ncol(residuals)
  Sigma <- crossprod(residuals) / TT
  log_det <- log(det(Sigma))

  aic <- log_det + (2 * n_params) / TT
  bic <- log_det + (log(TT) * n_params) / TT

  return(list(aic = aic, bic = bic))
}


# ───────────────────────────────────────────────────────────────────────────────
# 6.  Console Banner
# ───────────────────────────────────────────────────────────────────────────────

#' Print a formatted section banner to the console.
#' @param text  Character string; the title text
print_banner <- function(text) {
  width <- max(nchar(text) + 6, 60)
  line  <- paste(rep("=", width), collapse = "")
  cat("\n", line, "\n", "   ", text, "\n", line, "\n\n", sep = "")
}


# ───────────────────────────────────────────────────────────────────────────────
# 9.  Deterministic Column Builder
# ───────────────────────────────────────────────────────────────────────────────

#' Build deterministic regressor columns (intercept and/or trend).
#'
#' @param n_obs         Integer; number of observations (rows)
#' @param deterministic Character; one of "none", "intercept", "trend", "both"
#' @param t_offset      Integer; starting value for the trend counter (default 1)
#' @return A list with:
#'   \item{D_mat}{n_obs x d matrix of deterministic columns (NULL if d = 0)}
#'   \item{n_det}{Integer; number of deterministic terms (0, 1, or 2)}
#'   \item{has_intercept}{Logical}
#'   \item{has_trend}{Logical}
build_deterministic_columns <- function(n_obs, deterministic = "intercept",
                                         t_offset = 1) {
  deterministic <- match.arg(deterministic,
                              c("none", "intercept", "trend", "both"))
  cols <- list()
  has_intercept <- FALSE
  has_trend     <- FALSE

  if (deterministic %in% c("intercept", "both")) {
    cols[["const"]] <- rep(1, n_obs)
    has_intercept <- TRUE
  }
  if (deterministic %in% c("trend", "both")) {
    cols[["trend"]] <- seq(t_offset, length.out = n_obs)
    has_trend <- TRUE
  }

  n_det <- length(cols)
  if (n_det == 0) {
    D_mat <- NULL
  } else {
    D_mat <- do.call(cbind, cols)
  }

  list(D_mat = D_mat, n_det = n_det,
       has_intercept = has_intercept, has_trend = has_trend)
}


# ───────────────────────────────────────────────────────────────────────────────
# 10.  Deterministic Coefficient Partitioner
# ───────────────────────────────────────────────────────────────────────────────

#' Extract deterministic-term coefficients from the top rows of a beta matrix.
#'
#' @param beta          Full coefficient matrix (rows x k_dom)
#' @param n_det         Number of deterministic rows (0, 1, or 2)
#' @param has_intercept Logical
#' @param has_trend     Logical
#' @param k_dom         Number of domestic (dependent) variables
#' @return A list with:
#'   \item{intercept}{k_dom-vector (zeros if no intercept)}
#'   \item{trend}{k_dom-vector (zeros if no trend)}
#'   \item{start_idx}{Row index where stochastic coefficients begin}
partition_deterministic <- function(beta, n_det, has_intercept, has_trend,
                                     k_dom) {
  idx <- 1
  intercept  <- rep(0, k_dom)
  trend_coef <- rep(0, k_dom)

  if (has_intercept) {
    intercept <- beta[idx, ]
    idx <- idx + 1
  }
  if (has_trend) {
    trend_coef <- beta[idx, ]
    idx <- idx + 1
  }

  list(intercept = intercept, trend = trend_coef, start_idx = idx)
}


# ───────────────────────────────────────────────────────────────────────────────
# 11.  Variable Recovery  (model space → economic variables)
# ───────────────────────────────────────────────────────────────────────────────

#' Convert GVAR forecast / history matrices from model-space transformations
#' back to interpretable economic variables.
#'
#' The model stores variables in three possible forms:
#'   - "*_logdiff"  : quarterly log-change  → annualised % growth = × 400
#'   - "*_log"      : log level             → exponentiate to recover level
#'   - "*_level" / raw levels (rate, lt_rate, oil_level, …)  → pass through
#'
#' @param mat       T × K matrix of model-space values.  Columns must be named
#'                  with the convention "UNIT.varname" (as in gvar_model$var_names).
#' @param freq      "quarterly" (default) or "annual".  Controls annualisation.
#' @param log_base_list  Named list: for "*_log" variables, supply the last
#'                  observed level so that the forecast log can be exponentiated.
#'                  Element name must match the full "UNIT.varname" string.
#'                  Pass NULL to skip level recovery for log variables.
#' @return  A data frame with the same row structure as mat and columns:
#'   variable, period (1…T), model_value, econ_value, econ_unit
recover_economic_variables <- function(mat, freq = "quarterly",
                                       log_base_list = NULL) {

  mat   <- as.matrix(mat)
  vnames <- colnames(mat)
  if (is.null(vnames)) stop("mat must have named columns (UNIT.varname).")

  annualise <- if (freq == "quarterly") 400 else 100

  rows <- list()
  for (j in seq_along(vnames)) {
    vfull <- vnames[j]
    # Strip the "UNIT." prefix to get the raw variable name
    dot    <- regexpr("\\.", vfull)[1]
    varstr <- substr(vfull, dot + 1, nchar(vfull))

    for (t in seq_len(nrow(mat))) {
      mv <- mat[t, j]
      ev <- NA_real_
      eu <- "raw"

      if (grepl("_logdiff$", varstr)) {
        # Annualised percentage growth rate
        ev <- mv * annualise
        eu <- "% growth (ann.)"

      } else if (grepl("_log$", varstr)) {
        # Log level: convert to index level if base supplied
        if (!is.null(log_base_list) && vfull %in% names(log_base_list)) {
          ev <- exp(mv) * log_base_list[[vfull]]
        } else {
          ev <- exp(mv)
        }
        eu <- "index level"

      } else {
        # Rates, levels, oil — pass through unchanged
        ev <- mv
        eu <- if (grepl("rate|lt_rate", varstr)) "% p.a." else "level"
      }

      rows[[length(rows) + 1]] <- data.frame(
        variable    = vfull,
        period      = t,
        model_value = mv,
        econ_value  = ev,
        econ_unit   = eu,
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, rows)
}


#' Quick helper: extract the last observed value (in model space) for each
#' variable in a sim_data list, returning a named list suitable for
#' log_base_list in recover_economic_variables().
#'
#' @param sim_data  Named list of T × k matrices (as built by to_gvar_list)
#' @return          Named list of scalars: "UNIT.varname" → last_value
last_observed_levels <- function(sim_data) {
  out <- list()
  for (u in names(sim_data)) {
    mat <- sim_data[[u]]
    for (v in colnames(mat)) {
      out[[paste0(u, ".", v)]] <- mat[nrow(mat), v]
    }
  }
  out
}
