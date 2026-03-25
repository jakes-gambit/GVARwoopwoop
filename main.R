###############################################################################
#  main.R  –  GVAR Toolkit: Master Script
#
#  Pipeline structure:
#
#    0.  Source modules & load packages
#    1.  Configuration
#    2.  Data loading / generation
#    3.  Pre-estimation: unit-root tests, cointegration tests
#    4.  Lag-order selection (AIC / BIC)
#    5.  Prepare the GVAR dataset
#    ─── ESTIMATION ───────────────────────────────────────────────
#    6a. Frequentist GVAR / GVECM
#    6b. Bayesian GVAR (Minnesota prior)
#    ─── DIAGNOSTICS (In-Sample + OOS) ────────────────────────────
#    7a. Frequentist diagnostics (R², LB, JB, OOS, RMSE, plots)
#    7b. Bayesian diagnostics   (DIC, ML, PPC + insample + OOS + plots)
#    ─── IMPULSE RESPONSE FUNCTIONS ───────────────────────────────
#    8a. Frequentist GIRFs (bootstrap)
#    8b. Bayesian GIRFs (posterior credible intervals)
#    8c. Bayesian vs Frequentist IRF comparison
#    ─── CONDITIONAL FORECASTING ──────────────────────────────────
#    9a. Frequentist: Kalman filter
#    9b. Frequentist: Waggoner & Zha (hard + soft constraints)
#    9c. Frequentist: KF vs WZ comparison
#    9d. Bayesian conditional forecasting
#   10.  Summary
#
#  To use with your own data, replace Section 2 with your data-loading code.
#  The only requirements are:
#    (a) A named list of T × k_i  matrices  (one per country / unit)
#    (b) An N × N  weight matrix (e.g. bilateral trade shares)
#
#  Author:  GVAR Toolkit
#  Date:    2026
###############################################################################


# ═══════════════════════════════════════════════════════════════════════════════
# 0.  Source All Module Files
# ═══════════════════════════════════════════════════════════════════════════════

# Set working directory to the project root (adjust if needed)
setwd("C:/Users/WZHJEAB/Desktop/Conditional Forecasting 1/")

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

# Load required packages
load_packages()

cat("\n")
message("=============================================================")
message("   GVAR Toolkit – Full Estimation Pipeline")
message("   (VAR / VECM / Bayesian GVAR)")
message("=============================================================\n")

library(tidyr)


# ═══════════════════════════════════════════════════════════════════════════════
# 1.  Configuration
# ═══════════════════════════════════════════════════════════════════════════════

# Frequency of the data: "quarterly" or "yearly"
FREQ <- "quarterly"

# Maximum lags to search over (for AIC/BIC selection)
MAX_P <- 8     # domestic lags
MAX_Q <- 4     # foreign lags

# Lag selection criterion: "aic" or "bic"
LAG_CRITERION <- "bic"

# Deterministic terms in the short-run dynamics (applies to VAR, VECM, Bayesian)
# Options: "none" | "intercept" | "trend" | "both"
DETERMINISTIC <- "intercept"

# IRF settings
IRF_HORIZON <- 5*4        # number of periods
IRF_BOOT    <- 2000      # bootstrap replications
IRF_CI      <- 0.95      # confidence level

# OOS settings
OOS_H       <- 1         # forecast horizon
OOS_T0_FRAC <- 0.70      # fraction used for initial estimation

# Global variable configuration (Pesaran et al. 2004 dominant-unit approach)
# Set to NULL to disable global variable treatment
GLOBAL_VARS <- list(
  var_names     = c("oil"),    # variables treated as global
  dominant_unit = "DEU",       # unit where oil is endogenous
  d_lag         = 1            # lags for global exogenous in non-dominant units
)

# Bayesian GVAR settings (Minnesota / Litterman prior)
BAYES_DRAWS   <- 1000      # number of posterior draws
BAYES_LAMBDA1 <- 0.01       # overall tightness
BAYES_LAMBDA2 <- 0.01       # cross-variable tightness
BAYES_LAMBDA3 <- 0.01       # lag decay rate
BAYES_LAMBDA4 <- 0.01      # star variable tightness
BAYES_RW      <- TRUE       # random-walk prior (TRUE) or white-noise (FALSE)

# Cointegration test settings
COINT_K       <- 2          # VAR lags for Johansen procedure
COINT_ECDET   <- "const"    # deterministic specification
COINT_SIG     <- 0.05       # significance level


# ═══════════════════════════════════════════════════════════════════════════════
# 2.  Data Loading  (Replace This Section With Your Own Data)
# ═══════════════════════════════════════════════════════════════════════════════

set.seed(123)
TT <- 120    # number of observations (quarters)
N  <- 3      # number of countries
k  <- 2      # domestic variables per country
selected_countries <- c("DEU", "ESP", 'FRA' ,"ITA")
var_names     <- c('gdp','inf', 'rate','lt_rate','reer', 'oil')

sim_data <- gvar_data_with_oil[names(gvar_data_with_oil) %in% selected_countries]
sim_data <- align_to_common_dates(sim_data)

# NOTE: With global variable support, we no longer need to exclude oil manually.
# Oil is handled by the dominant-unit approach in prepare_gvar_dataset().
# If you do NOT want global variable treatment, set GLOBAL_VARS <- NULL above
# and uncomment the exclusion below:
exclude_vars <- c('rate')
sim_data <- lapply(sim_data, function(df) {
  df[, !(colnames(df) %in% exclude_vars), drop = FALSE]
})

# Bilateral trade-weight matrix (illustrative)
W_raw <- matrix(c(
  0,    0.2, 0.1, 0.1,  # US exports to EU
  0.1, 0, 0.2, 0.2,
  0.3, 0.1, 0, 0.3,
  0.1,0.2,0.4,0# EU exports to US
), nrow = length(selected_countries), byrow = TRUE,
dimnames = list(selected_countries, selected_countries))

message("[Example] Data loaded: ", length(sim_data), " countries.\n")


# ═══════════════════════════════════════════════════════════════════════════════
# 3.  Pre-Estimation Tests
# ═══════════════════════════════════════════════════════════════════════════════

# ── 3a. Unit-Root Tests ──────────────────────────────────────────────────────

ur_results <- panel_unit_root(sim_data)

# Optional: automatically difference I(1) series
# sim_data <- auto_difference(sim_data, ur_results)

# ── 3b. Cointegration Tests (Johansen) ──────────────────────────────────────

# Build weight matrix and star variables first (needed for cointegration test)
W <- build_weight_matrix(W_raw)
star_list <- compute_star_vars(sim_data, W,
                                global_var_names = GLOBAL_VARS$var_names)

message("[Preview] Star variables for first unit (first 5 rows):")
print(head(star_list[[1]], 5))

# Run Johansen cointegration tests on each unit's (x_it, x*_it) system
USE_VECM <- FALSE
coint_results <- NULL

if (any(ur_results$consensus == "I(1)")) {
  coint_results <- panel_cointegration(
    data_list = sim_data,
    star_list = star_list,
    K         = COINT_K,
    ecdet     = COINT_ECDET,
    sig_level = COINT_SIG
  )
  USE_VECM <- any(coint_results$rank_vec > 0)

  if (USE_VECM) {
    message("[GVAR] Cointegration detected – will estimate GVECM.")
  } else {
    message("[GVAR] No cointegration found – will estimate standard GVAR.")
  }
} else {
  message("[GVAR] All variables appear I(0) – skipping cointegration tests.")
}


# ═══════════════════════════════════════════════════════════════════════════════
# 4.  Lag-Order Selection via AIC / BIC
# ═══════════════════════════════════════════════════════════════════════════════

lag_sel <- select_lags_all(
  data_list     = sim_data,
  star_list     = star_list,
  max_p         = MAX_P,
  max_q         = MAX_Q,
  criterion     = LAG_CRITERION,
  deterministic = DETERMINISTIC
)

# Use the recommended lags (or override manually)
p_vec <- lag_sel$p_vec
q_vec <- lag_sel$q_vec


# ═══════════════════════════════════════════════════════════════════════════════
# 5.  Prepare the GVAR Dataset
# ═══════════════════════════════════════════════════════════════════════════════

gvar_data <- prepare_gvar_dataset(
  data_list     = sim_data,
  weights       = W_raw,
  p_lag         = p_vec,
  q_lag         = q_vec,
  freq          = FREQ,
  global_vars   = GLOBAL_VARS,
  deterministic = DETERMINISTIC
)


###############################################################################
#                          E S T I M A T I O N
###############################################################################

# ═══════════════════════════════════════════════════════════════════════════════
# 6a. Frequentist GVAR / GVECM Estimation
# ═══════════════════════════════════════════════════════════════════════════════

if (USE_VECM) {
  message("\n[GVAR] Estimating Global VECM (error-correction form) ...")
  gvar_model <- estimate_gvecm(
    gvar_data     = gvar_data,
    rank_vec      = coint_results$rank_vec,
    ecdet         = COINT_ECDET,
    deterministic = DETERMINISTIC
  )
} else {
  message("\n[GVAR] Estimating standard GVAR (levels VAR) ...")
  gvar_model <- estimate_gvar(gvar_data, deterministic = DETERMINISTIC)
}

# Print eigenvalue summary (stability check)
cat("\n  Companion matrix eigenvalue moduli (top 5):\n  ")
top_eig <- sort(gvar_model$eig_mod, decreasing = TRUE)
cat(round(top_eig[1:min(5, length(top_eig))], 4), "\n")

# Check B_0 coefficients (contemporaneous star effects)
for (u in names(gvar_model$unit_fits)) {
  B0 <- gvar_model$unit_fits[[u]]$B[[1]]
  cat(u, "- max |B_0|:", round(max(abs(B0)), 3), "\n")
}


# ═══════════════════════════════════════════════════════════════════════════════
# 6b. Bayesian GVAR Estimation (Minnesota Prior)
# ═══════════════════════════════════════════════════════════════════════════════

bayesian_model <- bayesian_estimate_gvar(
  gvar_data     = gvar_data,
  n_draws       = BAYES_DRAWS,
  lambda_1      = BAYES_LAMBDA1,
  lambda_2      = BAYES_LAMBDA2,
  lambda_3      = BAYES_LAMBDA3,
  lambda_4      = BAYES_LAMBDA4,
  rw_prior      = BAYES_RW,
  seed          = 42,
  deterministic = DETERMINISTIC
)

message(sprintf("[Bayesian] %d / %d draws are stable.",
                bayesian_model$n_stable, bayesian_model$n_draws))


###############################################################################
#                D I A G N O S T I C S  (In-Sample + OOS)
###############################################################################

# ═══════════════════════════════════════════════════════════════════════════════
# 7a. Frequentist Diagnostics
#     R², Ljung-Box, Jarque-Bera, OOS recursive forecast, RMSE, MAE, Theil's U
#     Plots: residuals, fitted vs actual, eigenvalues, CUSUM
# ═══════════════════════════════════════════════════════════════════════════════

diag_results <- run_all_diagnostics(
  gvar_model = gvar_model,
  data_list  = sim_data,
  star_list  = star_list,
  W          = W,
  h          = OOS_H,
  t0_frac    = OOS_T0_FRAC
)


# ═══════════════════════════════════════════════════════════════════════════════
# 7b. Bayesian Diagnostics
#     Bayesian-specific: DIC, marginal likelihood, posterior predictive checks
#     + In-sample (R², LB, JB on posterior mean)
#     + OOS recursive forecast + evaluation (RMSE, MAE, Theil's U)
#     + Plots: residuals, fitted vs actual, eigenvalues (point + posterior
#       distribution), CUSUM, prior vs posterior densities
# ═══════════════════════════════════════════════════════════════════════════════

bayesian_full_diag <- run_bayesian_all_diagnostics(
  bayesian_gvar = bayesian_model,
  gvar_data     = gvar_data,
  data_list     = sim_data,
  star_list     = star_list,
  W             = W,
  h             = OOS_H,
  t0_frac       = OOS_T0_FRAC
)

# Prior vs posterior density plot
plot_prior_posterior(bayesian_model, unit_name = "USA")


###############################################################################
#             I M P U L S E   R E S P O N S E   F U N C T I O N S
###############################################################################

# ═══════════════════════════════════════════════════════════════════════════════
# 8a. Frequentist GIRFs (Bootstrap Confidence Bands)
# ═══════════════════════════════════════════════════════════════════════════════

irf_us_gdp <- bootstrap_girf(
  gvar_model = gvar_model,
  shock_var  = "USA.gdp",
  horizon    = IRF_HORIZON,
  n_boot     = IRF_BOOT,
  ci_level   = IRF_CI
)

plot_irf(irf_us_gdp, shock_label = "US GDP Shock (Frequentist)")


# ═══════════════════════════════════════════════════════════════════════════════
# 8b. Bayesian GIRFs (Posterior Credible Intervals)
# ═══════════════════════════════════════════════════════════════════════════════

bayesian_irf <- bayesian_girf(
  bayesian_gvar = bayesian_model,
  shock_var     = "USA.gdp",
  horizon       = IRF_HORIZON,
  ci_level      = 0.90
)

plot_bayesian_irf(bayesian_irf, shock_label = "US GDP Shock (Bayesian)")


# ═══════════════════════════════════════════════════════════════════════════════
# 8c. Bayesian vs Frequentist IRF Comparison
# ═══════════════════════════════════════════════════════════════════════════════

bayesian_model_comparison(bayesian_irf, irf_us_gdp)


###############################################################################
#              C O N D I T I O N A L   F O R E C A S T I N G
###############################################################################

# Define conditions (shared across methods)
conditions <- data.frame(
  variable  = c(rep("USA.gdp", 3)),
  horizon   = c(1, 2, 3),
  value     = c(-3.5, -4, -1),
  type      = "hard",
  tolerance = c(0, 0, 0)
)
cat("\n  Conditions imposed:\n")
print(conditions)


# ═══════════════════════════════════════════════════════════════════════════════
# 9a. Frequentist: Kalman Filter Conditional Forecast
# ═══════════════════════════════════════════════════════════════════════════════

cf <- conditional_forecast(
  gvar_model = gvar_model,
  conditions = conditions,
  max_h      = 8,
  data_list  = sim_data,
  ci_level   = 0.95
)

cat("\n  Conditional Forecast Means (KF):\n")
print(round(cf$forecast_mean, 4))
plot_conditional_forecast(cf)


# ═══════════════════════════════════════════════════════════════════════════════
# 9b. Frequentist: Waggoner & Zha (1999) Conditional Forecast
# ═══════════════════════════════════════════════════════════════════════════════

# ── Hard constraints ─────────────────────────────────────────────────────────

conditions_wz <- data.frame(
  variable  = c(rep("USA.gdp", 3)),
  horizon   = c(1, 2, 3),
  value     = c(-3.5, -4, -1),
  type      = "hard",
  tolerance = 0
)

wz_hard <- conditional_forecast_wz(
  gvar_model = gvar_model,
  conditions = conditions_wz,
  max_h      = 8,
  data_list  = sim_data,
  n_draws    = 2000,
  ci_level   = 0.95
)

cat("\n  WZ Conditional Forecast Means (Hard Constraints):\n")
print(round(wz_hard$forecast_mean, 4))
plot_conditional_forecast_wz(wz_hard, show_uncond = TRUE)

# ── Soft constraints ─────────────────────────────────────────────────────────

conditions_soft <- data.frame(
  variable  = c(rep("USA.gdp", 3)),
  horizon   = c(1, 2, 3),
  value     = c(-3.5, -4, -1),
  type      = "soft",
  tolerance = c(0.1, 0.1, 0.1)
)

wz_soft <- conditional_forecast_wz(
  gvar_model = gvar_model,
  conditions = conditions_soft,
  max_h      = 8,
  data_list  = sim_data,
  n_draws    = 2000,
  ci_level   = 0.90
)

cat("\n  WZ Soft Conditional Forecast Means:\n")
print(round(wz_soft$forecast_mean, 4))
plot_conditional_forecast_wz(wz_soft, show_uncond = TRUE)


# ═══════════════════════════════════════════════════════════════════════════════
# 9c. Frequentist: Kalman Filter vs Waggoner-Zha Comparison
# ═══════════════════════════════════════════════════════════════════════════════

compare_kf_wz(kf_result = cf, wz_result = wz_hard)

cat("\n  Max absolute discrepancy (KF mean vs WZ point):\n")
disc <- abs(cf$forecast_mean[1:8, ] - wz_hard$forecast_point[1:8, ])
cat("  ", round(max(disc), 6), "\n")


# ═══════════════════════════════════════════════════════════════════════════════
# 9d. Bayesian Conditional Forecasting
# ═══════════════════════════════════════════════════════════════════════════════

bayesian_cf <- bayesian_conditional_forecast(
  bayesian_gvar = bayesian_model,
  conditions    = conditions,
  max_h         = 8,
  data_list     = sim_data,
  method        = "kalman",
  ci_level      = 0.90,
  n_cf_draws    = 500
)

cat("\n  Bayesian Conditional Forecast Means:\n")
print(round(bayesian_cf$forecast_mean, 4))
plot_bayesian_forecast(bayesian_cf)


###############################################################################
#                           S U M M A R Y
###############################################################################

message("\n[GVAR] Pipeline complete.  All results stored in:")
message("")
message("  ESTIMATION:")
message("  - gvar_model         : frequentist GVAR / GVECM object")
message("  - bayesian_model     : Bayesian GVAR object (Minnesota prior)")
message("")
message("  DIAGNOSTICS:")
message("  - diag_results       : frequentist (insample + OOS + plots)")
message("  - bayesian_full_diag : Bayesian (DIC, ML, PPC + insample + OOS + plots)")
message("")
message("  IMPULSE RESPONSES:")
message("  - irf_us_gdp         : frequentist GIRFs (bootstrap)")
message("  - bayesian_irf       : Bayesian GIRFs (posterior credible intervals)")
message("")
message("  CONDITIONAL FORECASTS:")
message("  - cf                 : frequentist (Kalman filter)")
message("  - wz_hard            : frequentist (Waggoner-Zha, hard constraints)")
message("  - wz_soft            : frequentist (Waggoner-Zha, soft constraints)")
message("  - bayesian_cf        : Bayesian conditional forecast")
message("")
message("  PRE-ESTIMATION:")
message("  - ur_results         : unit-root test results")
message("  - coint_results      : Johansen cointegration test results")
message("  - lag_sel            : lag-selection results")
message(sprintf("  - DETERMINISTIC      : '%s'", DETERMINISTIC))
message("")
message("  Note: Content generated using AI – expert verification")
message("  is strongly recommended.\n")
