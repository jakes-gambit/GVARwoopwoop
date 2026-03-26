###############################################################################
#  02_data_preparation.R  –  Data Preparation for GVAR
#
#  A GVAR model treats a panel of N countries (or units), each observed over
#  T periods on k_i domestic variables.  Foreign ("star") variables for each
#  unit are trade-weighted (or otherwise weighted) averages of the other
#  units' domestic variables.
#
#  Contents
#  --------
#  1. validate_gvar_data()    – sanity checks on inputs
#  1b. align_to_common_dates() – align all units to common time period
#  2. build_weight_matrix()   – construct & row-normalise a weight matrix
#  3. compute_star_vars()     – create foreign-star variables  x*_it
#  4. prepare_country_data()  – assemble (x_it, x*_it) for one unit
#  5. prepare_gvar_dataset()  – master function that loops over all units
###############################################################################


# ───────────────────────────────────────────────────────────────────────────────
# 1.  Input Validation
# ───────────────────────────────────────────────────────────────────────────────

#' Validate the raw inputs supplied by the user before estimation.
#'
#' @param data_list   Named list of length N.  Each element is a T × k_i
#'                    matrix (or data.frame) of domestic variables for unit i.
#'                    Row names or a "date" column may carry the time index.
#' @param weights     N × N numeric matrix of bilateral linkage weights.
#'                    Diagonal elements should be zero; the function will
#'                    row-normalise if needed.
#' @param freq        Character; one of "quarterly" or "yearly".
#' @return            Invisible TRUE if all checks pass; stops with an
#'                    informative error otherwise.
validate_gvar_data <- function(data_list, weights, freq = c("quarterly", "yearly")) {
  
  freq <- match.arg(freq)
  N <- length(data_list)
  
  # --- basic structure ---
  if (!is.list(data_list) || is.null(names(data_list)))
    stop("data_list must be a *named* list of matrices / data.frames.")
  
  if (!is.matrix(weights) || nrow(weights) != N || ncol(weights) != N)
    stop("weights must be an N x N matrix (N = ", N, ").")
  
  # --- matching unit names ---
  unit_names <- names(data_list)
  if (!all(unit_names %in% rownames(weights)) ||
      !all(unit_names %in% colnames(weights))) {
    stop("Row/column names of the weight matrix must match names(data_list).")
  }
  
  # --- equal time dimensions ---
  TT_vec <- sapply(data_list, nrow)
  if (length(unique(TT_vec)) != 1)
    stop("All units must have the same number of time-series observations (T).",
         "\n  Found: ", paste(TT_vec, collapse = ", "))
  
  message(sprintf("[GVAR] Data validated: N=%d units, T=%d obs, freq=%s",
                  N, TT_vec[1], freq))
  invisible(TRUE)
}


# ───────────────────────────────────────────────────────────────────────────────
# 1b.  Align Units to Common Time Period
# ───────────────────────────────────────────────────────────────────────────────

#' Align all units to a common time period based on row names (dates).
#'
#' When units have different time coverage, this function identifies the
#' overlapping period and subsets all units to those common dates.
#'
#' @param data_list   Named list of data.frames/matrices with date row names
#' @param date_col    Optional: name of a date column (if dates aren't row names)
#' @return            Named list with all units aligned to common dates
align_to_common_dates <- function(data_list, date_col = NULL) {
  
  unit_names <- names(data_list)
  N <- length(unit_names)
  
  # Report current state
  row_counts <- sapply(data_list, nrow)
  message("[GVAR] Current row counts per unit:")
  for (u in unit_names) {
    message(sprintf("  %s: %d rows", u, nrow(data_list[[u]])))
  }
  
  # Extract dates from each unit
  get_dates <- function(df) {
    if (!is.null(date_col) && date_col %in% colnames(df)) {
      return(as.character(df[[date_col]]))
    } else if (!is.null(rownames(df))) {
      return(rownames(df))
    } else {
      return(as.character(seq_len(nrow(df))))
    }
  }
  
  all_dates <- lapply(data_list, get_dates)
  
  # Find common dates across all units
  common_dates <- Reduce(intersect, all_dates)
  
  if (length(common_dates) == 0) {
    stop("No common dates found across all units. Check your date formats.")
  }
  
  # Sort dates (assumes dates are sortable as strings, e.g., "2000Q1" or "2000-01-01")
  common_dates <- sort(common_dates)
  
  message(sprintf("[GVAR] Common period: %s to %s (%d observations)",
                  common_dates[1], common_dates[length(common_dates)],
                  length(common_dates)))
  
  # Subset each unit to common dates
  aligned_list <- setNames(vector("list", N), unit_names)
  
  for (u in unit_names) {
    df <- data_list[[u]]
    dates_u <- get_dates(df)
    
    # Find indices of common dates in this unit
    idx <- which(dates_u %in% common_dates)
    
    # Subset and preserve structure
    aligned_list[[u]] <- df[idx, , drop = FALSE]
    
    # Ensure row names are set
    if (is.null(rownames(aligned_list[[u]]))) {
      rownames(aligned_list[[u]]) <- common_dates
    }
  }
  
  # Verify alignment
  new_counts <- sapply(aligned_list, nrow)
  if (length(unique(new_counts)) != 1) {
    stop("Alignment failed. Row counts still differ: ",
         paste(new_counts, collapse = ", "))
  }
  
  message(sprintf("[GVAR] All %d units aligned to %d common observations.",
                  N, new_counts[1]))
  
  return(aligned_list)
}


# ───────────────────────────────────────────────────────────────────────────────
# 2.  Weight Matrix Construction
# ───────────────────────────────────────────────────────────────────────────────

#' Build and row-normalise a weight matrix.
#'
#' The user may supply raw bilateral flows (e.g. trade volumes).  This
#' function zeroes the diagonal and row-normalises so each row sums to 1.
#'
#' @param raw_weights  N × N numeric matrix (e.g. bilateral trade flows)
#' @return             N × N row-normalised weight matrix W with w_ii = 0
build_weight_matrix <- function(raw_weights) {
  W <- as.matrix(raw_weights)
  
  # Zero-out the diagonal (a unit cannot "trade with itself")
  diag(W) <- 0
  
  # Row-normalise: each row sums to 1
  row_sums <- rowSums(W)
  if (any(row_sums == 0))
    stop("At least one unit has zero total weight to all others.")
  
  W <- W / row_sums
  
  message("[GVAR] Weight matrix constructed and row-normalised (",
          nrow(W), " x ", ncol(W), ").")
  return(W)
}


# ───────────────────────────────────────────────────────────────────────────────
# 3.  Compute Foreign ("Star") Variables
# ───────────────────────────────────────────────────────────────────────────────

#' Create the foreign-star variables x*_{it} for every unit i.
#'
#' For unit i with k domestic variables, x*_{it} is a k-vector where each
#' element is the weighted average of the same variable across all OTHER units:
#'
#'    x*_{it,v} = sum_{j ≠ i}  w_{ij} * x_{jt,v}
#'
#' We restrict star variables to the *common* set of variable names that
#' appear in ALL units.  If a unit has extra domestic variables they are kept
#' domestically but no star counterpart is created.
#'
#' @param data_list  Named list of T × k_i matrices (one per unit)
#' @param W          N × N row-normalised weight matrix
#' @return           Named list of T × k_common matrices of star variables
compute_star_vars <- function(data_list, W, global_var_names = NULL) {
  
  unit_names <- names(data_list)
  N <- length(unit_names)
  TT <- nrow(data_list[[1]])
  
  # Identify the common variable names across all units
  all_vars <- lapply(data_list, colnames)
  common_vars <- Reduce(intersect, all_vars)
  
  # Exclude global variables from star construction (they are handled separately)
  if (!is.null(global_var_names)) {
    common_vars <- setdiff(common_vars, global_var_names)
    message("[GVAR] Excluded from star construction (global): ",
            paste(global_var_names, collapse = ", "))
  }
  
  if (length(common_vars) == 0)
    stop("No common variable names found across units.")
  
  message("[GVAR] Common variables for star construction: ",
          paste(common_vars, collapse = ", "))
  
  # For each common variable, stack all units into a T × N matrix,
  # then multiply by the weight matrix to get star values.
  star_list <- setNames(vector("list", N), unit_names)
  
  # Ensure W is a numeric matrix with proper names
  W <- as.matrix(W)
  storage.mode(W) <- "numeric"
  
  for (i in seq_len(N)) {
    star_mat <- matrix(0, nrow = TT, ncol = length(common_vars))
    colnames(star_mat) <- paste0(common_vars, "_star")
    
    for (v in seq_along(common_vars)) {
      vname <- common_vars[v]
      
      # Build T × N matrix for variable v explicitly (more robust than sapply)
      X_v <- matrix(NA_real_, nrow = TT, ncol = N)
      colnames(X_v) <- unit_names
      
      for (j in seq_len(N)) {
        unit_df <- data_list[[unit_names[j]]]
        
        # Check if variable exists
        if (!(vname %in% colnames(unit_df))) {
          stop("Variable '", vname, "' not found in unit '", unit_names[j], "'")
        }
        
        col_data <- unit_df[, vname, drop = TRUE]
        
        # Handle factors and ensure numeric
        if (is.factor(col_data)) {
          col_data <- as.numeric(as.character(col_data))
        } else {
          col_data <- as.numeric(col_data)
        }
        
        # Verify length matches
        if (length(col_data) != TT) {
          stop("Length mismatch for variable '", vname, "' in unit '", unit_names[j],
               "': expected ", TT, " rows but got ", length(col_data),
               ".\nHint: Use align_to_common_dates() to align all units first.")
        }
        
        X_v[, j] <- col_data
      }
      
      # Get weight vector for unit i (ensure numeric vector)
      w_i <- as.numeric(W[unit_names[i], unit_names]) 
      
      # Star variable for unit i:  x*_{it,v} = sum_j w_{ij} x_{jt,v}
      star_mat[, v] <- X_v %*% w_i
    }
    star_list[[unit_names[i]]] <- star_mat
  }
  
  return(star_list)
}


# ───────────────────────────────────────────────────────────────────────────────
# 4.  Assemble Data for a Single Unit
# ───────────────────────────────────────────────────────────────────────────────

#' Combine domestic and star variables for one unit, applying the chosen
#' number of lags for domestic (p) and foreign (q) variables.
#'
#' The individual VARX* model for unit i is:
#'
#'    x_{it} = c_i  +  sum_{l=1}^{p} A_{il} x_{i,t-l}
#'                   +  sum_{l=0}^{q} B_{il} x*_{i,t-l}  +  u_{it}
#'
#' Note: the contemporaneous star variable (l=0) is included by default.
#'
#' @param domestic  T × k_i matrix of domestic variables
#' @param star      T × k*  matrix of star variables
#' @param p_lag     Integer; number of own lags
#' @param q_lag     Integer; number of foreign lags (excluding contemporaneous)
#' @return          List with Y (dependent), X (regressors), and metadata
prepare_country_data <- function(domestic, star, p_lag, q_lag,
                                 global_exog   = NULL, d_lag = NULL,
                                 covid_dummies = NULL) {  # [NEW] crisis dummies
  
  domestic <- as.matrix(domestic)
  star     <- as.matrix(star)
  TT       <- nrow(domestic)
  k_dom    <- ncol(domestic)
  k_star   <- ncol(star)
  
  # Global exogenous handling
  k_global <- 0
  if (!is.null(global_exog)) {
    global_exog <- as.matrix(global_exog)
    k_global <- ncol(global_exog)
    if (is.null(d_lag)) d_lag <- q_lag
  } else {
    d_lag <- 0
  }
  
  max_lag <- max(p_lag, q_lag, d_lag)
  
  # Dependent variable: x_{it} from t = max_lag+1 ... T
  Y <- domestic[(max_lag + 1):TT, , drop = FALSE]
  n_eff <- nrow(Y)   # effective sample size
  
  # --- Regressors ---
  X_parts <- list()
  
  # (a) Contemporaneous star: x*_{it}
  X_parts[["star_contemp"]] <- star[(max_lag + 1):TT, , drop = FALSE]
  
  # (a2) Contemporaneous global exogenous: d_{it}
  if (k_global > 0) {
    X_parts[["global_contemp"]] <- global_exog[(max_lag + 1):TT, , drop = FALSE]
  }
  
  # (b) Lagged domestic variables
  for (l in seq_len(p_lag)) {
    rows <- (max_lag + 1 - l):(TT - l)
    block <- domestic[rows, , drop = FALSE]
    colnames(block) <- paste0(colnames(domestic), "_lag", l)
    X_parts[[paste0("dom_lag", l)]] <- block
  }
  
  # (c) Lagged star variables
  for (l in seq_len(q_lag)) {
    rows <- (max_lag + 1 - l):(TT - l)
    block <- star[rows, , drop = FALSE]
    colnames(block) <- paste0(colnames(star), "_lag", l)
    X_parts[[paste0("star_lag", l)]] <- block
  }
  
  # (d) Lagged global exogenous variables
  if (k_global > 0 && d_lag > 0) {
    for (l in seq_len(d_lag)) {
      rows <- (max_lag + 1 - l):(TT - l)
      block <- global_exog[rows, , drop = FALSE]
      colnames(block) <- paste0(colnames(global_exog), "_lag", l)
      X_parts[[paste0("global_lag", l)]] <- block
    }
  }

  # [NEW] (e) Contemporaneous COVID / crisis dummy variables
  # Dummies enter without lagging — they are outlier/pulse controls, not
  # dynamic regressors.  They are added last so existing coefficient indexing
  # (for stars and globals) is unchanged.
  k_dummy <- 0
  if (!is.null(covid_dummies)) {
    covid_dummies <- as.matrix(covid_dummies)
    if (nrow(covid_dummies) != TT)
      stop("covid_dummies row count (", nrow(covid_dummies), ") != TT (", TT, ").")
    dummy_block <- covid_dummies[(max_lag + 1):TT, , drop = FALSE]
    X_parts[["covid_dummies"]] <- dummy_block
    k_dummy <- ncol(covid_dummies)
  }

  X <- do.call(cbind, X_parts)

  return(list(
    Y        = Y,
    X        = X,
    k_dom    = k_dom,
    k_star   = k_star,
    k_global = k_global,
    k_dummy  = k_dummy,    # [NEW] number of dummy columns
    p_lag    = p_lag,
    q_lag    = q_lag,
    d_lag    = d_lag,
    n_eff    = n_eff
  ))
}


# ───────────────────────────────────────────────────────────────────────────────
# 4b.  [NEW] COVID / Crisis Dummy Construction
#      Creates impulse and step dummies for the COVID-19 episode (or any crisis)
#      that can be added as exogenous regressors to control for extreme outliers.
#      Dummies are aligned to the same row names as data_list matrices.
# ───────────────────────────────────────────────────────────────────────────────

#' Build a matrix of crisis dummy variables.
#'
#' Two types of dummies are generated:
#'   - "pulse"     : 1 in a specific quarter, 0 elsewhere  (outlier correction)
#'   - "step"      : 1 from a quarter onward               (permanent shift)
#'
#' The default specification matches COVID-19: a pulse for the two worst
#' quarters (2020-Q1, 2020-Q2) and a recovery pulse for 2020-Q3.
#' All periods are matched against the row names of the data matrices,
#' which follow the "YYYY-QN" format produced by 00_fetch_fred_data.R.
#'
#' @param date_labels   Character vector of "YYYY-QN" labels (= rownames of a
#'                      data_list matrix).  Must cover all crisis periods.
#' @param pulse_periods Character vector of "YYYY-QN" quarters that receive a
#'                      pulse dummy (default: 2020-Q1 and 2020-Q2).
#' @param step_periods  Named list; each element is a "YYYY-QN" start date for
#'                      a step dummy.  Name becomes the column name.
#'                      (default: NULL – no step dummies)
#' @return  T × (n_pulse + n_step) numeric matrix; colnames indicate the type
make_covid_dummies <- function(date_labels,
                               pulse_periods = c("2020-Q1", "2020-Q2"),
                               step_periods  = NULL) {

  T <- length(date_labels)
  dummies <- list()

  # ── Pulse dummies: 1 exactly at the specified quarter ──────────────────────
  for (qtr in pulse_periods) {
    d <- as.integer(date_labels == qtr)     # 1 if match, 0 otherwise
    if (sum(d) == 0) {
      message(sprintf("[COVID dummy] Pulse period '%s' not found in date labels – skipped.", qtr))
    }
    col_nm <- paste0("d_", gsub("-", "_", qtr))   # e.g. "d_2020_Q1"
    dummies[[col_nm]] <- d
  }

  # ── Step dummies: 1 from the start quarter onward ──────────────────────────
  if (!is.null(step_periods)) {
    for (nm in names(step_periods)) {
      start_qtr <- step_periods[[nm]]
      d <- as.integer(date_labels >= start_qtr)
      if (sum(d) == 0) {
        message(sprintf("[COVID dummy] Step period '%s' not found – skipped.", start_qtr))
      }
      dummies[[nm]] <- d
    }
  }

  if (length(dummies) == 0) {
    message("[COVID dummy] No valid dummy periods – returning NULL.")
    return(NULL)
  }

  D <- do.call(cbind, dummies)
  rownames(D) <- date_labels
  message(sprintf("[COVID dummy] Created %d dummy column(s): %s",
                  ncol(D), paste(colnames(D), collapse = ", ")))
  return(D)
}


# ───────────────────────────────────────────────────────────────────────────────
# 5.  Master Preparation Function
# ───────────────────────────────────────────────────────────────────────────────

#' Prepare the full GVAR dataset: validate, build stars, assemble per-unit
#' regression data.
#'
#' @param data_list       Named list of T × k_i matrices
#' @param weights         N × N raw weight matrix
#' @param p_lag           Integer or named integer vector of domestic lag orders
#' @param q_lag           Integer or named integer vector of foreign lag orders
#' @param freq            "quarterly" or "yearly"
#' @param global_vars     Optional list(var_names, dominant_unit, d_lag)
#' @param deterministic   "none" | "intercept" | "trend" | "both"
#' @param covid_dummies   Optional T × d matrix of crisis dummies.
#' @param star_restrict   Controls which star variables enter each unit's model.
#'   NULL (default): every unit receives all star variables (standard GVAR).
#'   Character vector: only these variable names (without _star suffix) enter
#'     as stars in every unit, e.g. c("gdp_log", "cpi_logdiff").
#'   Named list: per-unit overrides.  Each element is a character vector for
#'     that unit; units not listed receive the full star set.
#'     e.g. list(USA = c("gdp_log","cpi_logdiff"), CHN = c("gdp_log")).
#' @return                A list with unit_data, star_list, W, unit_names, etc.
prepare_gvar_dataset <- function(data_list, weights, p_lag = 1, q_lag = 1,
                                 freq = "quarterly",
                                 global_vars    = NULL,
                                 deterministic  = "intercept",
                                 covid_dummies  = NULL,
                                 star_restrict  = NULL) {
  
  # Step 1: validate
  validate_gvar_data(data_list, weights, freq)
  
  # Step 2: weight matrix
  W <- build_weight_matrix(weights)
  unit_names <- names(data_list)
  N <- length(unit_names)
  
  # Step 3: global variable configuration
  global_config <- NULL
  global_var_names <- NULL
  global_series <- NULL
  
  if (!is.null(global_vars)) {
    global_var_names <- global_vars$var_names
    dominant_unit    <- global_vars$dominant_unit
    d_lag            <- if (!is.null(global_vars$d_lag)) global_vars$d_lag else 1
    
    if (!(dominant_unit %in% unit_names))
      stop("Dominant unit '", dominant_unit, "' not found in data_list.")
    
    # Extract global variable series from the dominant unit's data
    dom_data <- as.matrix(data_list[[dominant_unit]])
    missing_gv <- setdiff(global_var_names, colnames(dom_data))
    if (length(missing_gv) > 0)
      stop("Global variable(s) not found in dominant unit: ",
           paste(missing_gv, collapse = ", "))
    
    global_series <- dom_data[, global_var_names, drop = FALSE]
    
    message(sprintf("[GVAR] Global variable(s): %s | Dominant unit: %s | d_lag: %d",
                    paste(global_var_names, collapse = ", "), dominant_unit, d_lag))
    
    global_config <- list(
      var_names     = global_var_names,
      dominant_unit = dominant_unit,
      d_lag         = d_lag,
      global_series = global_series
    )
  }
  
  # Step 4: star variables (excluding global variables from star construction)
  star_list <- compute_star_vars(data_list, W,
                                 global_var_names = global_var_names)
  
  # Allow per-unit lag orders (scalar is broadcast to all units)
  if (length(p_lag) == 1) p_lag <- setNames(rep(p_lag, N), unit_names)
  if (length(q_lag) == 1) q_lag <- setNames(rep(q_lag, N), unit_names)
  
  # [NEW] Step 4b: validate and align COVID dummies if provided
  if (!is.null(covid_dummies)) {
    covid_dummies <- as.matrix(covid_dummies)
    T_first <- nrow(data_list[[unit_names[1]]])
    if (nrow(covid_dummies) != T_first)
      stop("covid_dummies must have the same number of rows (T) as the data matrices. ",
           "Got ", nrow(covid_dummies), " vs ", T_first, ".")
    message(sprintf("[GVAR] COVID dummies included: %s",
                    paste(colnames(covid_dummies), collapse = ", ")))
  }

  # Step 5: assemble per-unit regression datasets
  unit_data <- setNames(vector("list", N), unit_names)
  for (u in unit_names) {
    # Determine global exogenous for this unit
    g_exog <- NULL
    g_d_lag <- NULL
    dom_u <- data_list[[u]]

    if (!is.null(global_config) && u != global_config$dominant_unit) {
      g_exog  <- global_config$global_series
      g_d_lag <- global_config$d_lag

      # Remove global variables from domestic set for non-dominant units.
      # The global variable enters only through the exogenous channel;
      # keeping it in domestic would create (near-)perfect collinearity
      # (e.g., oil is a world price — same series in both channels).
      gv_in_dom <- intersect(global_var_names, colnames(dom_u))
      if (length(gv_in_dom) > 0) {
        message(sprintf("  [%s] Dropping '%s' from domestic (enters as global exogenous).",
                        u, paste(gv_in_dom, collapse = "', '")))
        dom_u <- dom_u[, !(colnames(dom_u) %in% gv_in_dom), drop = FALSE]
      }
    }

    # Apply star variable restriction for this unit
    star_u <- star_list[[u]]
    if (!is.null(star_restrict)) {
      allowed_vars <- if (is.list(star_restrict)) star_restrict[[u]] else star_restrict
      if (!is.null(allowed_vars)) {
        allowed_cols <- paste0(allowed_vars, "_star")
        keep         <- intersect(allowed_cols, colnames(star_u))
        dropped      <- setdiff(colnames(star_u), keep)
        if (length(dropped) > 0)
          message(sprintf("  [%s] Star restricted to: %s (dropped: %s)",
                          u, paste(keep, collapse = ", "), paste(dropped, collapse = ", ")))
        star_u <- star_u[, keep, drop = FALSE]
      }
    }

    unit_data[[u]] <- prepare_country_data(
      domestic      = dom_u,
      star          = star_u,
      p_lag         = p_lag[u],
      q_lag         = q_lag[u],
      global_exog   = g_exog,
      d_lag         = g_d_lag,
      covid_dummies = covid_dummies   # [NEW] pass dummies into each unit
    )
  }

  message("[GVAR] Dataset preparation complete.")
  return(list(
    unit_data     = unit_data,
    star_list     = star_list,
    W             = W,
    unit_names    = unit_names,
    freq          = freq,
    global_config = global_config,
    deterministic = deterministic,
    covid_dummies = covid_dummies     # [NEW] stored for downstream reference
  ))
}
