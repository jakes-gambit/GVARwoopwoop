###############################################################################
#  main.R  –  GVAR Toolkit: Master Estimation Script
#
#  Pipeline:
#    0.  Source modules
#    1.  Configuration   ← edit here
#    2.  Data loading    (fred_data + bilateral trade weights)
#    3.  Unit-root & cointegration tests
#    4.  Lag selection
#    5.  Prepare GVAR dataset
#    6.  Estimation  (frequentist, with optional ridge regularisation)
#    7.  Diagnostics (in-sample + OOS)
#    8.  Impulse responses
#    9.  Conditional forecasting
#   10.  Local Projections (Jordà 2005)
###############################################################################

# ── 0. Source all modules ────────────────────────────────────────────────────

# setwd("C:/Users/WZHJEAB/Desktop/Conditional Forecasting 2/")   # adjust path

source("01_utils.R")
source("02_data_preparation.R")
source("03_unit_root_tests.R")
source("04_lag_selection.R")
source("05_gvar_estimation.R")
source("06_diagnostics.R")
source("07_impulse_response.R")
source("08_kalman_filter.R")
source("09_waggoner_zha.R")
source("10_vecm_estimation.R")
source("11_local_projections.R")
source("00b_fetch_weight_data.R")   # bilateral trade weight functions

load_packages()
library(dplyr)


###############################################################################
#  1.  CONFIGURATION  ─  edit this section
###############################################################################

# Countries to include (must be present in fred_data)
COUNTRIES <- c("USA", "DEU", "JPN", "FRA", "ITA")

# Domestic variables to use in each country's VARX*.
# Pick column names from fred_data (level / log / logdiff transformations).
# Stationary variables preferred: logdiff for GDP/CPI/REER/EQ, levels for rates.
DOMESTIC_VARS <- c("gdp_log", "cpi_logdiff", "rho_lms", "reer_log", "eq_log")

# Global variable(s): enter as exogenous for all non-dominant units.
# All variables listed in GLOBAL_VAR_NAMES must exist as columns in sim_data.
#
# OIL_VAR options:  "oil_level"   – price in $/barrel (I(1), use with caution)
#                   "oil_log"     – log price level    (I(1))
#                   "oil_logdiff" – quarterly log-change (stationary)
#                   "oil_diff"    – custom name used here
#
# To add a second global variable, extend GLOBAL_VAR_NAMES:
#   e.g. c("oil", "comm_logdiff")  — both must be columns in sim_data
OIL_VAR          <- "oil_diff"
OIL_DOM_UNIT     <- "USA"       # country where oil is endogenous
GLOBAL_VAR_NAMES <- c("oil")    # extend to c("oil","other_global") for multiples
GLOBAL_DOM_UNIT  <- OIL_DOM_UNIT
GLOBAL_D_LAG     <- 1           # lags for global exogenous channel

# Lags
MAX_P         <- 4        # max domestic lags for AIC/BIC search
MAX_Q         <- 2        # max foreign lags
LAG_CRITERION <- "bic"

# Deterministic terms: "none" | "intercept" | "trend" | "both"
DETERMINISTIC <- "intercept"

# COVID dummies: pulse dummies for specified quarters (set to NULL to disable)
COVID_DUMMIES <-c("2020-Q1", "2020-Q2")   # set NULL to turn off

# ── Ridge regularisation ──────────────────────────────────────────────────────
# Set RIDGE_LAMBDA > 0 to apply ridge (L2) shrinkage to the frequentist VAR.
# A value of 0 (default) means plain OLS – no regularisation.
# Typical starting range: 0.01 – 1.0.  Larger values impose more shrinkage.
# The intercept is never penalised.
RIDGE_LAMBDA  <- 0       # 0 = plain OLS; e.g. 0.1 for mild ridge

# OOS evaluation
OOS_H      <- 1
OOS_T0_FRAC <- 0.70

# IRF
IRF_HORIZON <- 200
IRF_BOOT    <- 0
IRF_CI      <- 0.95

# Cointegration (Johansen)
COINT_K    <- 2
COINT_DET  <- "const"
COINT_SIG  <- 0.05

# Weight matrix: years to average bilateral trade data over
WEIGHT_YEARS     <- 2014:2016
WEIGHT_FIN_ALPHA <- 0        # 0 = trade only; >0 blends BIS financial weights


###############################################################################
#  2.  DATA LOADING
###############################################################################

# Load pre-fetched FRED data (run 00_fetch_fred_data.R once to produce this)
load("gvar_fred_data.RData")   # loads: fred_data, to_gvar_list()

# Build GVAR list: domestic vars + selected oil transformation
all_vars <- c(DOMESTIC_VARS, OIL_VAR)
sim_data <- to_gvar_list(fred_data,
                          var_cols  = all_vars,
                          countries = COUNTRIES)

# Rename oil column to "oil" so the pipeline refers to it consistently
sim_data <- lapply(sim_data, function(m) {
  colnames(m)[colnames(m) == OIL_VAR] <- "oil"
  m
})

sim_data <- align_to_common_dates(sim_data)


###############################################################################
#  2b.  VARIABLE PLOTS
#        One panel per variable in DOMESTIC_VARS + OIL_VAR.
#        All countries overlaid as coloured lines.  Saved to PDF.
###############################################################################

plot_gvar_series(sim_data,
                 var_list  = c(DOMESTIC_VARS, "oil"),
                 file_out  = "gvar_variable_plots.pdf")


# COVID pulse dummies (appended as exogenous regressors, one column per quarter)
covid_dummies <- 
  if (!is.null(COVID_DUMMIES)) {
    make_covid_dummies(
      rownames(sim_data[[1]]),
      pulse_periods = COVID_DUMMIES
    )
  } else {
    NULL
  }


# # Bilateral trade weight matrix from empirical data
W_raw <- fetch_trade_weights(COUNTRIES, avg_years = WEIGHT_YEARS,
                               finance_alpha = WEIGHT_FIN_ALPHA)

# Global variable config (supports multiple globals via GLOBAL_VAR_NAMES vector)
GLOBAL_VARS <- list(var_names     = GLOBAL_VAR_NAMES,
                    dominant_unit = GLOBAL_DOM_UNIT,
                    d_lag         = GLOBAL_D_LAG)

message(sprintf("[Data] %d countries, %d observations each, variables: %s",
                length(sim_data), nrow(sim_data[[1]]),
                paste(colnames(sim_data[[1]]), collapse=", ")))


###############################################################################
#  3.  PRE-ESTIMATION TESTS
###############################################################################

ur_results   <- panel_unit_root(sim_data)

W        <- build_weight_matrix(W_raw)
star_list <- compute_star_vars(sim_data, W, global_var_names = GLOBAL_VARS$var_names)

USE_VECM     <- FALSE
coint_results <- NULL

if (any(ur_results$consensus == "I(1)")) {
  coint_results <- panel_cointegration(sim_data, star_list,
                                       K = COINT_K, ecdet = COINT_DET,
                                       sig_level = COINT_SIG)
  USE_VECM <- any(coint_results$rank_vec > 0)
}
message(if (USE_VECM) "[GVAR] Cointegration detected – using GVECM."
        else          "[GVAR] No cointegration – using standard GVAR.")


###############################################################################
#  4.  LAG SELECTION
###############################################################################

lag_sel <- select_lags_all(sim_data, star_list,
                            max_p = MAX_P, max_q = MAX_Q,
                            criterion = LAG_CRITERION,
                            deterministic = DETERMINISTIC)
p_vec <- lag_sel$p_vec
q_vec <- lag_sel$q_vec


###############################################################################
#  5.  PREPARE GVAR DATASET
###############################################################################

gvar_data <- prepare_gvar_dataset(
  data_list     = sim_data,
  weights       = W_raw,
  p_lag         = p_vec,
  q_lag         = q_vec,
  freq          = "quarterly",
  global_vars   = GLOBAL_VARS,
  deterministic = DETERMINISTIC,
  covid_dummies = covid_dummies
)


###############################################################################
#  6.  ESTIMATION
###############################################################################

if (USE_VECM) {
  gvar_model <- estimate_gvecm(gvar_data, rank_vec = coint_results$rank_vec,
                                ecdet = COINT_DET, deterministic = DETERMINISTIC,
                                ridge_lambda = RIDGE_LAMBDA)
} else {
  gvar_model <- estimate_gvar(gvar_data,
                               deterministic = DETERMINISTIC,
                               ridge_lambda  = RIDGE_LAMBDA)
}

if (RIDGE_LAMBDA > 0) {
  message(sprintf("[GVAR] Ridge regularisation applied (lambda = %.4f).", RIDGE_LAMBDA))
}


###############################################################################
#  6b.  EIGENVALUE DIAGNOSTICS
#        Individual country companion eigenvalues + global companion eigenvalues.
#        Max modulus < 1 ↔ model stability.  Prints top N sorted by modulus.
###############################################################################

cat("\n", strrep("=", 70), "\n", sep = "")
cat("  INDIVIDUAL UNIT COMPANION EIGENVALUES\n")
cat(strrep("=", 70), "\n", sep = "")
for (u in gvar_model$unit_names) {
  fit <- gvar_model$unit_fits[[u]]
  if (length(fit$A) > 0 && fit$k_dom > 0) {
    local_comp <- make_companion(fit$A)
    print_eigenvalues(local_comp, label = u, top_n = fit$k_dom * fit$p_lag)
  }
}

cat("\n", strrep("=", 70), "\n", sep = "")
cat("  GLOBAL COMPANION EIGENVALUES\n")
cat(strrep("=", 70), "\n", sep = "")
print_eigenvalues(gvar_model$companion,
                  label = "Global GVAR",
                  top_n = min(30, nrow(gvar_model$companion)))


###############################################################################
#  7.  DIAGNOSTICS
###############################################################################

diag_results <- run_all_diagnostics(
  gvar_model, sim_data, star_list, W, h = OOS_H, t0_frac = OOS_T0_FRAC)


###############################################################################
#  8.  IMPULSE RESPONSES  (VAR-based GIRF)
###############################################################################

shock_var <- paste0(COUNTRIES[1], ".gdp_log")

irf_freq <- bootstrap_girf(gvar_model, shock_var,
                            horizon = 20, n_boot = IRF_BOOT,
                            ci_level = IRF_CI)
plot_irf(irf_freq, shock_label = shock_var)


###############################################################################
#  9.  CONDITIONAL FORECASTING
###############################################################################

conditions <- data.frame(
  variable  = rep(paste0(COUNTRIES[1], ".oil"), 20),
  horizon   = 1:20,
  value     = rep(12,20),
  type      = "hard",
  tolerance = 0
)

# ── 9a. Kalman filter ────────────────────────────────────────────────────────
cf_kf <- conditional_forecast(gvar_model, conditions, max_h = 20,
                               data_list = sim_data, ci_level = 0)
plot_conditional_forecast(cf_kf)

# ── 9b. Waggoner & Zha ───────────────────────────────────────────────────────
cf_wz <- conditional_forecast_wz(gvar_model, conditions, max_h = 20,
                                  data_list = sim_data, n_draws = 2000,
                                  ci_level = 0.95)
plot_conditional_forecast_wz(cf_wz, show_uncond = TRUE)

compare_kf_wz(cf_kf, cf_wz)


###############################################################################
# 10.  VARIABLE RECOVERY
#
#  Converts model-space matrices to interpretable economic units.
#  Two contexts are distinguished:
#    "forecast" – level recovery for *_log, exact rho→% p.a., ×400 for logdiffs
#    "irf"      – % change (×100) for *_log, Δ% p.a. (×400) for rho, ×100 for logdiffs
#
#  Covered: IRF point responses, KF conditional forecast, WZ conditional forecast.
###############################################################################

# Last observed levels — needed for absolute *_log level recovery in forecasts
log_bases <- last_observed_levels(sim_data)

# Pattern to filter the print output to variables of interest
ALL_VARS_PATTERN <- paste(c(DOMESTIC_VARS, GLOBAL_VARS$var_names), collapse = "|")

# ── 10a. IRF recovery ────────────────────────────────────────────────────────
econ_irf <- recover_economic_variables(
  irf_freq$point,
  freq          = "quarterly",
  log_base_list = NULL,    # no base needed for % change IRF responses
  context       = "irf"
)
message("\n[Recovery] IRF responses (shock: ", shock_var, ") in economic units:")
print(subset(econ_irf, grepl(ALL_VARS_PATTERN, variable)))

# ── 10b. Conditional forecast – Kalman filter ────────────────────────────────
econ_cf_kf <- recover_economic_variables(
  cf_kf$forecast_mean,
  freq          = "quarterly",
  log_base_list = log_bases,
  context       = "forecast"
)
message("\n[Recovery] Conditional forecast (KF) in economic units:")
print(subset(econ_cf_kf, grepl(ALL_VARS_PATTERN, variable)))

# ── 10c. Conditional forecast – Waggoner-Zha ────────────────────────────────
econ_cf_wz <- recover_economic_variables(
  cf_wz$gibbs$mean,
  freq          = "quarterly",
  log_base_list = log_bases,
  context       = "forecast"
)
message("\n[Recovery] Conditional forecast (WZ) in economic units:")
print(subset(econ_cf_wz, grepl(ALL_VARS_PATTERN, variable)))

# Wide-format convenience tables (one row per horizon, one column per unit.var)
spread_econ <- function(df) {
  stats::reshape(df[, c("variable", "period", "econ_value")],
                 idvar = "period", timevar = "variable", direction = "wide")
}
econ_irf_wide   <- spread_econ(econ_irf)
econ_cf_kf_wide <- spread_econ(econ_cf_kf)
econ_cf_wz_wide <- spread_econ(econ_cf_wz)
