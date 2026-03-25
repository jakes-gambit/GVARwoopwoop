###############################################################################
#  03_unit_root_tests.R  –  Unit-Root Testing for GVAR Variables
#
#  Before estimating a GVAR it is important to check the integration order
#  of every variable in every unit.  This module provides:
#
#  1. adf_test_single()     – ADF test on a single series
#  2. kpss_test_single()    – KPSS stationarity test on a single series
#  3. test_unit_root()      – run ADF + KPSS for one unit's variables
#  4. panel_unit_root()     – loop over all units and build a summary table
#  5. auto_difference()     – difference any I(1) series and flag them
###############################################################################


# ───────────────────────────────────────────────────────────────────────────────
# 1.  ADF Test Wrapper
# ───────────────────────────────────────────────────────────────────────────────

#' Augmented Dickey-Fuller test on a single numeric vector.
#'
#' Uses urca::ur.df() which allows selection of "drift" (constant) or "trend"
#' (constant + linear trend) deterministic components, plus the number of
#' augmented lags.
#'
#' @param x        Numeric vector (time series)
#' @param type     Character; "drift" or "trend"
#' @param lags     Integer; number of augmented lags (if NULL, chosen by AIC)
#' @param sig_level  Significance level for the decision (default 0.05)
#' @return         Named list:  statistic, critical values, p_approx,
#'                 is_stationary (logical), lags_used
adf_test_single <- function(x, type = "drift", lags = NULL, sig_level = 0.05) {

  x <- as.numeric(na.omit(x))

  # Default lag selection: floor(12*(T/100)^{1/4})  (Schwert rule)
  if (is.null(lags)) {
    lags <- max(1, floor(12 * (length(x) / 100)^0.25))
  }
  # Ensure lags don't exceed feasible range

  lags <- min(lags, floor(length(x) / 3))

  test <- urca::ur.df(x, type = type, lags = lags, selectlags = "AIC")

  # Extract test statistic (tau) and critical values
  # Handle both matrix and vector returns from urca (depends on test type)
  tau_stat <- if (is.matrix(test@teststat)) test@teststat[1, 1] else test@teststat[1]
  crit     <- if (is.matrix(test@cval)) test@cval[1, ] else test@cval

  # Decision at chosen significance level
  crit_val <- crit[paste0(sig_level * 100, "pct")]
  if (length(crit_val) == 0) crit_val <- crit["5pct"]   # fallback
  is_stat  <- tau_stat < crit_val                       # reject H0 => stationary

  return(list(
    test_name      = "ADF",
    statistic      = tau_stat,
    critical       = crit,
    is_stationary  = is_stat,
    lags_used      = test@lags
  ))
}


# ───────────────────────────────────────────────────────────────────────────────
# 2.  KPSS Test Wrapper
# ───────────────────────────────────────────────────────────────────────────────

#' KPSS stationarity test on a single series.
#'
#' H0 for KPSS: the series is (trend-) stationary.
#' If the test statistic exceeds the critical value we REJECT stationarity.
#'
#' @param x         Numeric vector
#' @param type      "mu" (level stationarity) or "tau" (trend stationarity)
#' @param sig_level Significance level (default 0.05)
#' @return          Named list analogous to adf_test_single()
kpss_test_single <- function(x, type = "mu", sig_level = 0.05) {

  x <- as.numeric(na.omit(x))
  test <- urca::ur.kpss(x, type = type)

  # Handle both matrix and vector returns from urca (depends on test type)
  stat <- if (is.matrix(test@teststat)) test@teststat[1, 1] else test@teststat[1]
  crit <- if (is.matrix(test@cval)) test@cval[1, ] else test@cval

  # For KPSS, REJECT stationarity if stat > critical value
  crit_val <- crit[paste0(sig_level * 100, "pct")]
  if (length(crit_val) == 0) crit_val <- crit["5pct"]
  is_stat <- stat < crit_val   # fail to reject => stationary

  return(list(
    test_name     = "KPSS",
    statistic     = stat,
    critical      = crit,
    is_stationary = is_stat
  ))
}


# ───────────────────────────────────────────────────────────────────────────────
# 3.  Test All Variables of One Unit
# ───────────────────────────────────────────────────────────────────────────────

#' Run ADF and KPSS tests on every column of a unit's data matrix.
#'
#' @param data_mat  T × k matrix of domestic variables
#' @param unit_name Character label for the unit
#' @return          Data frame with one row per variable, columns:
#'                  unit, variable, adf_stat, adf_stationary,
#'                  kpss_stat, kpss_stationary, consensus
test_unit_root <- function(data_mat, unit_name = "Unit") {

  data_mat <- as.matrix(data_mat)
  var_names <- colnames(data_mat)
  if (is.null(var_names)) var_names <- paste0("V", seq_len(ncol(data_mat)))

  results <- data.frame(
    unit            = character(),
    variable        = character(),
    adf_stat        = numeric(),
    adf_stationary  = logical(),
    kpss_stat       = numeric(),
    kpss_stationary = logical(),
    consensus       = character(),
    stringsAsFactors = FALSE
  )

  for (v in seq_along(var_names)) {
    x <- data_mat[, v]

    adf_res  <- adf_test_single(x)
    kpss_res <- kpss_test_single(x)

    # Consensus logic:
    #   Both say stationary      => "I(0)"
    #   Both say non-stationary  => "I(1)"
    #   Conflicting              => "Inconclusive"
    if (adf_res$is_stationary & kpss_res$is_stationary) {
      consensus <- "I(0)"
    } else if (!adf_res$is_stationary & !kpss_res$is_stationary) {
      consensus <- "I(1)"
    } else {
      consensus <- "Inconclusive"
    }

    results <- rbind(results, data.frame(
      unit            = unit_name,
      variable        = var_names[v],
      adf_stat        = round(adf_res$statistic, 4),
      adf_stationary  = adf_res$is_stationary,
      kpss_stat       = round(kpss_res$statistic, 4),
      kpss_stationary = kpss_res$is_stationary,
      consensus       = consensus,
      stringsAsFactors = FALSE
    ))
  }

  return(results)
}


# ───────────────────────────────────────────────────────────────────────────────
# 4.  Panel-Wide Unit-Root Summary
# ───────────────────────────────────────────────────────────────────────────────

#' Run unit-root tests for ALL units and return a consolidated table.
#'
#' @param data_list  Named list of T × k_i matrices (raw data)
#' @return           Data frame with columns unit, variable, adf_stat, etc.
panel_unit_root <- function(data_list) {

  print_banner("Unit-Root Tests (ADF + KPSS)")

  all_results <- data.frame()

  for (u in names(data_list)) {
    res <- test_unit_root(data_list[[u]], unit_name = u)
    all_results <- rbind(all_results, res)
  }

  # Pretty-print
  message("  Summary of integration order across all units/variables:\n")
  print(all_results, row.names = FALSE)

  n_i1 <- sum(all_results$consensus == "I(1)")
  n_i0 <- sum(all_results$consensus == "I(0)")
  n_inc <- sum(all_results$consensus == "Inconclusive")

  message(sprintf(
    "\n  => I(0): %d | I(1): %d | Inconclusive: %d\n", n_i0, n_i1, n_inc))

  return(all_results)
}


# ───────────────────────────────────────────────────────────────────────────────
# 5.  Automatic Differencing
# ───────────────────────────────────────────────────────────────────────────────

#' Difference any series diagnosed as I(1) and return a cleaned data list.
#'
#' Variables labelled "I(1)" or "Inconclusive" are first-differenced.
#' Variables labelled "I(0)" are left in levels.
#'
#' @param data_list   Named list of T × k_i matrices
#' @param ur_results  Data frame from panel_unit_root()
#' @param diff_inconclusive  Logical; if TRUE, treat inconclusive as I(1)
#' @return            Named list of (T-1) × k_i matrices of transformed data
auto_difference <- function(data_list, ur_results, diff_inconclusive = TRUE) {

  message("[GVAR] Applying automatic differencing where needed...")

  out <- list()

  for (u in names(data_list)) {
    mat <- as.matrix(data_list[[u]])
    new_mat <- matrix(NA, nrow = nrow(mat) - 1, ncol = ncol(mat))
    colnames(new_mat) <- colnames(mat)

    for (v in colnames(mat)) {
      row <- ur_results[ur_results$unit == u & ur_results$variable == v, ]

      needs_diff <- FALSE
      if (nrow(row) > 0) {
        if (row$consensus == "I(1)") needs_diff <- TRUE
        if (row$consensus == "Inconclusive" & diff_inconclusive) needs_diff <- TRUE
      }

      if (needs_diff) {
        new_mat[, v] <- diff(mat[, v])
        message(sprintf("  [%s] %s  =>  differenced", u, v))
      } else {
        new_mat[, v] <- mat[-1, v]   # drop first obs to keep conformable
        message(sprintf("  [%s] %s  =>  kept in levels", u, v))
      }
    }

    out[[u]] <- new_mat
  }

  return(out)
}


# ───────────────────────────────────────────────────────────────────────────────
# 6.  Johansen Cointegration Test for a Single Unit
# ───────────────────────────────────────────────────────────────────────────────

#' Johansen cointegration test on the combined (x_it, x*_it) system for one unit.
#'
#' Applies the Johansen trace and maximum-eigenvalue tests to determine the
#' cointegrating rank r_i.  The rank is selected by sequential testing from
#' r = 0 upward: reject H0(r <= j) if test stat > critical value, then
#' increment j.  Stop when the null is not rejected.
#'
#' @param z_mat     T x m numeric matrix where m = k_i + k*_i
#'                  (domestic + star variables stacked column-wise)
#' @param type      Character; which test statistic to use for rank decision:
#'                  "trace" (default) or "eigen" (maximum eigenvalue)
#' @param ecdet     Deterministic specification: "const" (restricted constant,
#'                  Case III – most common), "trend", or "none"
#' @param K         Integer; number of VAR lags in levels for the Johansen
#'                  procedure (VECM has K-1 lagged differences).  Default 2.
#' @param sig_level Significance level for rank decision (default 0.05)
#' @return A list with:
#'   \item{rank}{Selected cointegrating rank r_i}
#'   \item{trace_stats}{Trace test statistics for r = 0, ..., m-1}
#'   \item{trace_crit}{Corresponding critical values at sig_level}
#'   \item{eigen_stats}{Max-eigenvalue test statistics}
#'   \item{eigen_crit}{Corresponding critical values}
#'   \item{eigenvalues}{Ordered eigenvalues from the Johansen procedure}
#'   \item{beta}{m x rank matrix of cointegrating vectors (NULL if rank = 0)}
#'   \item{alpha}{m x rank loading matrix (NULL if rank = 0)}
#'   \item{ca_jo_trace}{Raw urca::ca.jo object (trace test)}
#'   \item{ca_jo_eigen}{Raw urca::ca.jo object (eigenvalue test)}
johansen_test_single <- function(z_mat,
                                  type = c("trace", "eigen"),
                                  ecdet = c("const", "trend", "none"),
                                  K = 2,
                                  sig_level = 0.05) {

  type  <- match.arg(type)
  ecdet <- match.arg(ecdet)

  z_mat <- as.matrix(z_mat)
  z_mat <- z_mat[complete.cases(z_mat), , drop = FALSE]
  m     <- ncol(z_mat)
  TT    <- nrow(z_mat)

  # Ensure K is feasible
  K <- max(2, min(K, floor(TT / (m + 2))))

  # Run both trace and eigenvalue tests
  ca_trace <- urca::ca.jo(z_mat, type = "trace", ecdet = ecdet, K = K,
                           spec = "transitory")
  ca_eigen <- urca::ca.jo(z_mat, type = "eigen", ecdet = ecdet, K = K,
                           spec = "transitory")

  # Extract test statistics and critical values
  # ca.jo stores statistics in reverse order (r = m-1, ..., 0)
  trace_stats <- as.numeric(ca_trace@teststat)
  trace_cval  <- ca_trace@cval
  eigen_stats <- as.numeric(ca_eigen@teststat)
  eigen_cval  <- ca_eigen@cval

  # Identify the critical value column for chosen significance
  crit_col <- paste0(sig_level * 100, "pct")
  if (!(crit_col %in% colnames(trace_cval))) crit_col <- "5pct"

  trace_crit <- as.numeric(trace_cval[, crit_col])
  eigen_crit <- as.numeric(eigen_cval[, crit_col])

  # The statistics are ordered from r = 0 at the END to r = m-1 at the START
  # in urca::ca.jo.  Reverse to get natural ordering r = 0, 1, ..., m-1.
  trace_stats <- rev(trace_stats)
  trace_crit  <- rev(trace_crit)
  eigen_stats <- rev(eigen_stats)
  eigen_crit  <- rev(eigen_crit)

  # Sequential rank selection
  select_rank <- function(stats, crits) {
    r <- 0
    for (j in seq_along(stats)) {
      if (stats[j] > crits[j]) {
        r <- j
      } else {
        break
      }
    }
    return(r)
  }

  rank_trace <- min(select_rank(trace_stats, trace_crit), m - 1)
  rank_eigen <- min(select_rank(eigen_stats, eigen_crit), m - 1)

  # Use the rank from the chosen test type
  rank_used <- if (type == "trace") rank_trace else rank_eigen

  # Extract cointegrating vectors (beta) and loading matrix (alpha)
  beta  <- NULL
  alpha <- NULL
  if (rank_used > 0) {
    # ca.jo stores beta in @V and alpha in @W
    # When ecdet = "const" or "trend", @V has (m+1) rows; the last row is
    # the restricted deterministic coefficient.  We keep only the first m rows
    # so that beta is m × rank (matching the variable space).
    V_full <- ca_trace@V
    W_full <- ca_trace@W
    beta  <- V_full[1:m, 1:rank_used, drop = FALSE]
    alpha <- W_full[1:m, 1:rank_used, drop = FALSE]
  }

  # Eigenvalues
  eigenvalues <- ca_trace@lambda

  return(list(
    rank        = rank_used,
    rank_trace  = rank_trace,
    rank_eigen  = rank_eigen,
    trace_stats = trace_stats,
    trace_crit  = trace_crit,
    eigen_stats = eigen_stats,
    eigen_crit  = eigen_crit,
    eigenvalues = eigenvalues,
    beta        = beta,
    alpha       = alpha,
    ca_jo_trace = ca_trace,
    ca_jo_eigen = ca_eigen
  ))
}


# ───────────────────────────────────────────────────────────────────────────────
# 7.  Panel-Wide Cointegration Summary
# ───────────────────────────────────────────────────────────────────────────────

#' Run Johansen cointegration tests for ALL units and return a summary.
#'
#' For each unit, the test is applied to the combined system z_it = (x_it, x*_it).
#' The resulting cointegrating rank r_i determines whether a VECM or VAR in
#' levels should be estimated for that unit.
#'
#' @param data_list  Named list of T x k_i matrices (domestic variables)
#' @param star_list  Named list of T x k*_i matrices (star variables, from
#'                   compute_star_vars())
#' @param K          Number of VAR lags for the Johansen procedure (default 2)
#' @param ecdet      Deterministic specification: "const", "trend", or "none"
#' @param sig_level  Significance level (default 0.05)
#' @param type       Which test statistic for rank decision: "trace" or "eigen"
#' @return A list with:
#'   \item{rank_vec}{Named integer vector of cointegrating ranks per unit}
#'   \item{summary_table}{Data frame with unit, m, rank_trace, rank_eigen, rank_used}
#'   \item{details}{Named list of full johansen_test_single() results per unit}
panel_cointegration <- function(data_list, star_list,
                                 K = 2,
                                 ecdet = "const",
                                 sig_level = 0.05,
                                 type = "trace") {

  print_banner("Johansen Cointegration Tests")

  unit_names <- names(data_list)
  N <- length(unit_names)

  rank_vec <- setNames(integer(N), unit_names)
  details  <- setNames(vector("list", N), unit_names)

  summary_rows <- list()

  for (u in unit_names) {
    dom  <- as.matrix(data_list[[u]])
    star <- as.matrix(star_list[[u]])
    z    <- cbind(dom, star)
    m    <- ncol(z)

    message(sprintf("  Testing unit '%s' (m = %d variables) ...", u, m))

    res <- johansen_test_single(z, type = type, ecdet = ecdet,
                                 K = K, sig_level = sig_level)

    rank_vec[u]  <- res$rank
    details[[u]] <- res

    summary_rows[[u]] <- data.frame(
      unit       = u,
      m          = m,
      rank_trace = res$rank_trace,
      rank_eigen = res$rank_eigen,
      rank_used  = res$rank,
      stringsAsFactors = FALSE
    )
  }

  summary_table <- do.call(rbind, summary_rows)
  rownames(summary_table) <- NULL

  message("\n  Cointegration rank summary:\n")
  print(summary_table, row.names = FALSE)

  n_coint <- sum(rank_vec > 0)
  n_none  <- sum(rank_vec == 0)
  message(sprintf(
    "\n  => Cointegrated: %d | No cointegration: %d\n", n_coint, n_none))

  return(list(
    rank_vec      = rank_vec,
    summary_table = summary_table,
    details       = details
  ))
}
