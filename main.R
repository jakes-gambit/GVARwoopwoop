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

setwd("C:/Users/WZHJEAB/Desktop/Conditional Forecasting 2/")   # adjust path

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
COUNTRIES <- c("USA", "DEU", "FRA")

# Domestic variables to use in each country's VARX*.
# Pick column names from fred_data (level / log / logdiff transformations).
# Stationary variables preferred: logdiff for GDP/CPI/REER/EQ, levels for rates.
DOMESTIC_VARS <- c("gdp_logdiff", "cpi_logdiff", "lt_rate", "reer_log",'eq_log')

# Global variable: one oil transformation enters via the dominant-unit approach.
# Options: "oil_level"  – price in $/barrel (I(1), use with unit-root caution)
#          "oil_log"    – log price level   (I(1))
#          "oil_logdiff"– quarterly log-change (stationary, lowest instability)
OIL_VAR       <- "oil_logdiff"
OIL_DOM_UNIT  <- "USA"    # country where oil is endogenous

# Lags
MAX_P         <- 4        # max domestic lags for AIC/BIC search
MAX_Q         <- 2        # max foreign lags
LAG_CRITERION <- "bic"

# Deterministic terms: "none" | "intercept" | "trend" | "both"
DETERMINISTIC <- "intercept"

# COVID dummies: pulse dummies for specified quarters (set to NULL to disable)
COVID_DUMMIES <- c("2020-Q1", "2020-Q2")   # set NULL to turn off

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
IRF_HORIZON <- 20
IRF_BOOT    <- 2000
IRF_CI      <- 0.95

# Local Projections
LP_HORIZON  <- 20         # max horizon for LP-IRFs
LP_N_LAGS   <- 1          # number of control lags in each LP regression
LP_CI       <- 0.95       # confidence level (Newey-West HAC bands)
LP_BOOT     <- 500        # bootstrap reps (set 0 to skip bootstrap)

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


# Bilateral trade weight matrix from empirical data
W_raw <- fetch_trade_weights(COUNTRIES, avg_years = WEIGHT_YEARS,
                              finance_alpha = WEIGHT_FIN_ALPHA)

# Global variable config
GLOBAL_VARS <- list(var_names = "oil", dominant_unit = OIL_DOM_UNIT, d_lag = 1)

message(sprintf("[Data] %d countries, %d observations each, variables: %s",
                length(sim_data), nrow(sim_data[[1]]),
                paste(colnames(sim_data[[1]]), collapse=", ")))


###############################################################################
#  3.  PRE-ESTIMATION TESTS
###############################################################################

ur_results   <- panel_unit_root(sim_data)

W        <- build_weight_matrix(W_raw)
star_list <- compute_star_vars(sim_data, W, global_var_names = "oil")

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
                                ecdet = COINT_DET, deterministic = DETERMINISTIC)
} else {
  gvar_model <- estimate_gvar(gvar_data,
                               deterministic = DETERMINISTIC,
                               ridge_lambda  = RIDGE_LAMBDA)
}

if (RIDGE_LAMBDA > 0) {
  message(sprintf("[GVAR] Ridge regularisation applied (lambda = %.4f).", RIDGE_LAMBDA))
}


###############################################################################
#  7.  DIAGNOSTICS
###############################################################################

diag_results <- run_all_diagnostics(
  gvar_model, sim_data, star_list, W, h = OOS_H, t0_frac = OOS_T0_FRAC)


###############################################################################
#  8.  IMPULSE RESPONSES  (VAR-based GIRF)
###############################################################################

shock_var <- paste0(COUNTRIES[1], ".gdp_logdiff")

irf_freq <- bootstrap_girf(gvar_model, shock_var,
                            horizon = IRF_HORIZON, n_boot = IRF_BOOT,
                            ci_level = IRF_CI)
plot_irf(irf_freq, shock_label = shock_var)


###############################################################################
#  9.  CONDITIONAL FORECASTING
###############################################################################

conditions <- data.frame(
  variable  = rep(paste0(COUNTRIES[1], ".gdp_logdiff"), 3),
  horizon   = 1:3,
  value     = c(-0.03, -0.04, -0.01),
  type      = "hard",
  tolerance = 0
)

# ── 9a. Kalman filter ────────────────────────────────────────────────────────
cf_kf <- conditional_forecast(gvar_model, conditions, max_h = 8,
                               data_list = sim_data, ci_level = 0.95)
plot_conditional_forecast(cf_kf)

# ── 9b. Waggoner & Zha ───────────────────────────────────────────────────────
cf_wz <- conditional_forecast_wz(gvar_model, conditions, max_h = 8,
                                  data_list = sim_data, n_draws = 2000,
                                  ci_level = 0.95)
plot_conditional_forecast_wz(cf_wz, show_uncond = TRUE)

compare_kf_wz(cf_kf, cf_wz)


###############################################################################
# 10.  LOCAL PROJECTIONS  (Jordà 2005)
###############################################################################

# ── 10a. Analytic LP-IRFs with Newey-West HAC bands ──────────────────────────
lp_result <- lp_gvar(gvar_model, shock_var,
                      data_list = sim_data,
                      n_lags    = LP_N_LAGS,
                      horizon   = LP_HORIZON,
                      ci_level  = LP_CI)

plot_lp_irf(lp_result, shock_label = shock_var)

# ── 10b. Bootstrap LP-IRFs (wild bootstrap) ───────────────────────────────────
if (LP_BOOT > 0) {
  lp_boot <- bootstrap_lp(gvar_model, shock_var,
                           data_list = sim_data,
                           n_lags    = LP_N_LAGS,
                           horizon   = LP_HORIZON,
                           n_boot    = LP_BOOT,
                           ci_level  = LP_CI)
  plot_lp_irf(lp_boot, shock_label = paste0(shock_var, " (Bootstrap LP)"))

  # ── 10c. Compare LP vs VAR-GIRF ──────────────────────────────────────────
  compare_lp_girf(lp_boot, irf_freq, shock_label = shock_var)
}
