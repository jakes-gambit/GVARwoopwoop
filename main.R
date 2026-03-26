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
#    6.  Estimation  (frequentist + Bayesian)
#    7.  Diagnostics (in-sample + OOS)
#    8.  Impulse responses
#    9.  Conditional forecasting
###############################################################################

# ── 0. Source all modules ────────────────────────────────────────────────────

setwd("C:/Users/WZHJEAB/Desktop/Conditional Forecasting 1/")   # adjust path

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
source("11_bayesian_gvar.R")
source("12_bayesian_diagnostics.R")
source("00b_fetch_weight_data.R")   # bilateral trade weight functions

load_packages()
library(dplyr)


###############################################################################
#  1.  CONFIGURATION  ─  edit this section
###############################################################################

# Countries to include (must be present in fred_data)
COUNTRIES <- c("DEU", "ESP", "FRA", "ITA")

# Domestic variables to use in each country's VARX*.
# Pick column names from fred_data.  Log-levels are the GVAR standard (I(1) in
# levels, cointegration handled by GVECM if detected).
# Options per variable:  *_level (raw), *_log (log level), *_logdiff (stationary)
DOMESTIC_VARS <- c("gdp_log", "cpi_log", "lt_rate", "reer_log")

# Global variable: one oil transformation enters via the dominant-unit approach.
# Options: "oil_level"  – price in $/barrel (I(1), use with unit-root caution)
#          "oil_log"    – log price level   (I(1))
#          "oil_logdiff"– quarterly log-change (stationary, lowest instability)
OIL_VAR       <- "oil_logdiff"
OIL_DOM_UNIT  <- "DEU"    # country where oil is endogenous

# Lags
MAX_P         <- 4        # max domestic lags for AIC/BIC search
MAX_Q         <- 2        # max foreign lags
LAG_CRITERION <- "bic"

# Deterministic terms: "none" | "intercept" | "trend" | "both"
DETERMINISTIC <- "intercept"

# COVID dummies: pulse dummies for specified quarters (set to NULL to disable)
COVID_DUMMIES <- c("2020-Q1", "2020-Q2")   # set NULL to turn off

# Bayesian prior (Minnesota/Litterman)
BAYES_DRAWS          <- 1000
BAYES_LAMBDA1        <- 0.01    # overall tightness
BAYES_LAMBDA2        <- 0.01    # cross-variable tightness
BAYES_LAMBDA3        <- 0.01    # lag decay
BAYES_LAMBDA4        <- 0.01    # star-variable tightness
BAYES_RW             <- TRUE    # TRUE = random-walk prior (set FALSE for white-noise)
BAYES_REGULARISE     <- FALSE   # TRUE = reject explosive draws (rejection sampling)
BAYES_LAMBDA5        <- 0.05    # sum-of-coefficients prior (0 = off; 0.01–0.1 typical)
BAYES_PROJECT_STABLE <- FALSE   # TRUE = rescale explosive companion draws to sr=0.99

# OOS evaluation
OOS_H      <- 1
OOS_T0_FRAC <- 0.70

# IRF
IRF_HORIZON <- 20
IRF_BOOT    <- 2000
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

# COVID pulse dummies (appended as exogenous regressors, one column per quarter)
covid_dummies <- if (!is.null(COVID_DUMMIES))
  make_covid_dummies(rownames(sim_data[[1]]), pulse_periods = COVID_DUMMIES)
else NULL

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

# ── 6a. Frequentist ──────────────────────────────────────────────────────────

if (USE_VECM) {
  gvar_model <- estimate_gvecm(gvar_data, rank_vec = coint_results$rank_vec,
                                ecdet = COINT_DET, deterministic = DETERMINISTIC)
} else {
  gvar_model <- estimate_gvar(gvar_data, deterministic = DETERMINISTIC)
}

# ── 6b. Bayesian ─────────────────────────────────────────────────────────────

bayesian_model <- bayesian_estimate_gvar(
  gvar_data,
  n_draws        = BAYES_DRAWS,
  lambda_1       = BAYES_LAMBDA1, lambda_2 = BAYES_LAMBDA2,
  lambda_3       = BAYES_LAMBDA3, lambda_4 = BAYES_LAMBDA4,
  rw_prior       = BAYES_RW, seed = 42,
  deterministic  = DETERMINISTIC,
  regularise     = BAYES_REGULARISE,
  lambda_5       = BAYES_LAMBDA5,
  project_stable = BAYES_PROJECT_STABLE
)
message(sprintf("[Bayesian] %d / %d draws stable.",
                bayesian_model$n_stable, bayesian_model$n_draws))


###############################################################################
#  7.  DIAGNOSTICS
###############################################################################

diag_results <- run_all_diagnostics(
  gvar_model, sim_data, star_list, W, h = OOS_H, t0_frac = OOS_T0_FRAC)

bayesian_diag <- run_bayesian_all_diagnostics(
  bayesian_model, gvar_data, sim_data, star_list, W,
  h = OOS_H, t0_frac = OOS_T0_FRAC)


###############################################################################
#  8.  IMPULSE RESPONSES
###############################################################################

shock_var <- paste0(COUNTRIES[1], ".", DOMESTIC_VARS[1])

irf_freq <- bootstrap_girf(gvar_model, shock_var,
                            horizon = IRF_HORIZON, n_boot = IRF_BOOT,
                            ci_level = IRF_CI)
plot_irf(irf_freq, shock_label = paste(shock_var, "(Frequentist)"))

irf_bayes <- bayesian_girf(bayesian_model, shock_var,
                            horizon = IRF_HORIZON, ci_level = 0.90)
plot_bayesian_irf(irf_bayes, shock_label = paste(shock_var, "(Bayesian)"))

bayesian_model_comparison(irf_bayes, irf_freq)


###############################################################################
#  9.  CONDITIONAL FORECASTING
###############################################################################

# Conditional scenario: 3-quarter path for the first domestic variable of COUNTRIES[1].
# Values should be on the same scale as DOMESTIC_VARS[1] (log-levels by default).
conditions <- data.frame(
  variable  = rep(paste0(COUNTRIES[1], ".", DOMESTIC_VARS[1]), 3),
  horizon   = 1:3,
  value     = c(-0.01, -0.01, -0.005),   # illustrative log-level deviations; adjust as needed
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

# ── 9c. Bayesian conditional forecast ────────────────────────────────────────
cf_bayes <- bayesian_conditional_forecast(bayesian_model, conditions, max_h = 8,
                                           data_list = sim_data, ci_level = 0.90)
plot_bayesian_forecast(cf_bayes)
