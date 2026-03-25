# GVAR Toolkit

A modular R implementation of the **Global Vector Autoregressive (GVAR)** model with **conditional forecasting via the Kalman filter**.

## Features

| Feature | Description |
|---|---|
| **N-variable, N-unit GVAR** | Handles any number of countries/units with any number of domestic variables |
| **Flexible frequency** | Quarterly or yearly data |
| **Arbitrary lag orders** | Per-unit domestic (p) and foreign (q) lags |
| **AIC / BIC lag selection** | Joint grid search over (p, q) for each unit |
| **Unit-root testing** | ADF + KPSS with consensus classification (I(0)/I(1)) |
| **In-sample diagnostics** | R², Ljung-Box serial correlation, Jarque-Bera normality |
| **Out-of-sample evaluation** | Expanding-window recursive forecasts, RMSE, MAE, Theil's U |
| **Generalised IRFs** | Pesaran & Shin (1998) GIRFs with bootstrap confidence bands |
| **Conditional forecasting (KF)** | Kalman filter + smoother with hard or soft constraints |
| **Conditional forecasting (WZ)** | Waggoner & Zha (1999) shock-space projection with hard and soft constraints |
| **Method comparison** | Side-by-side KF vs WZ comparison with overlay plots |

---

## Project Structure

```
GVAR_Toolkit/
├── R/
│   ├── 01_utils.R                # Helper functions, OLS, companion matrix
│   ├── 02_data_preparation.R     # Weight matrix, star variables, stacking
│   ├── 03_unit_root_tests.R      # ADF, KPSS, panel summary, auto-differencing
│   ├── 04_lag_selection.R        # AIC/BIC grid search per unit
│   ├── 05_gvar_estimation.R      # VARX* estimation + global stacking
│   ├── 06_diagnostics.R          # In-sample & out-of-sample tests
│   ├── 07_impulse_response.R     # GIRFs with bootstrap confidence bands
│   ├── 08_kalman_filter.R        # Conditional forecasting (Kalman filter)
│   ├── 09_waggoner_zha.R        # Conditional forecasting (Waggoner & Zha 1999)
│   └── main.R                    # Master script (run this)
└── README.md
```

---

## Quick Start

```r
# 1. Open R and set working directory to GVAR_Toolkit/R
setwd("path/to/GVAR_Toolkit/R")

# 2. Run the full pipeline (includes a built-in simulated example)
source("main.R")
```

The `main.R` script will:
- Source all modules and install missing packages
- Simulate a 3-country, 2-variable dataset
- Run unit-root tests, select lags, estimate the GVAR
- Produce diagnostics, IRFs, and conditional forecasts

---

## Using Your Own Data

Replace Section 2 of `main.R` with your own data loading code. You need two objects:

### 1. `data_list` — a named list of data matrices

```r
data_list <- list(
  US = us_matrix,   # T × k_us  matrix with named columns
  EU = eu_matrix,   # T × k_eu  matrix
  CN = cn_matrix    # T × k_cn  matrix
)
```

- All matrices must have the **same number of rows** (T)
- Column names identify the variables (e.g. `"gdp"`, `"inf"`, `"ir"`)
- Common column names across units are used to construct star variables

### 2. `W_raw` — an N × N weight matrix

```r
W_raw <- matrix(
  c(0, 0.4, 0.6,
    0.5, 0, 0.5,
    0.55, 0.45, 0),
  nrow = 3, byrow = TRUE,
  dimnames = list(c("US","EU","CN"), c("US","EU","CN"))
)
```

- Diagonal elements should be zero (self-weight is zeroed automatically)
- The toolkit will row-normalise for you

---

## Conditional Forecasting: Waggoner & Zha (1999)

The WZ approach works in the **space of future structural shocks** rather than treating conditions as observations of a latent state. This has several advantages:

- **Exact constraint satisfaction**: hard constraints hold draw-by-draw, not just in expectation
- **Clean subspace decomposition**: constrained and unconstrained shock directions are separated explicitly
- **Structural identification**: naturally uses the Cholesky (or other) identification of the residual covariance

### How It Works

The h-step forecast can be written as a linear function of stacked future shocks:

```
y = mu + M * eta
```

where `M` is the structural moving-average (MA) matrix. Conditions impose `R * y = r`, which in shock space becomes `R * M * eta = r - R * mu`. The conditional draw is:

```
eta* = (I - Proj) * eta_free  +  (RM)' * [(RM)(RM)']^{-1} * (r - R * mu)
```

where `Proj` is the projection onto the constrained subspace and `eta_free ~ N(0, I)`.

### Hard Constraints

```r
conditions <- data.frame(
  variable = c("US.gdp", "US.gdp", "EU.inf"),
  horizon  = c(1, 2, 4),
  value    = c(0.5, 0.5, -0.2),
  type     = "hard"
)

wz <- conditional_forecast_wz(
  gvar_model = gvar_model,
  conditions = conditions,
  max_h      = 8,
  data_list  = sim_data,
  n_draws    = 2000,
  ci_level   = 0.90
)

plot_conditional_forecast_wz(wz, show_uncond = TRUE)
```

### Soft Constraints

Soft constraints allow uncertainty around the imposed values. The `tolerance` parameter controls how much deviation from the target is permitted (standard deviation units):

```r
conditions_soft <- data.frame(
  variable  = c("US.gdp", "EU.inf"),
  horizon   = c(1, 4),
  value     = c(0.5, -0.2),
  type      = "soft",
  tolerance = c(0.10, 0.25)   # 0.10 = tight, 0.25 = loose
)

wz_soft <- conditional_forecast_wz(
  gvar_model = gvar_model,
  conditions = conditions_soft,
  max_h      = 8,
  data_list  = sim_data,
  n_draws    = 2000
)
```

### Comparing Kalman Filter vs Waggoner-Zha

The toolkit includes a built-in comparison function:

```r
# Run both methods with the same conditions, then overlay
cf_kf <- conditional_forecast(gvar_model, conditions, max_h = 8,
                               data_list = sim_data, ci_level = 0.90)
cf_wz <- conditional_forecast_wz(gvar_model, conditions, max_h = 8,
                                  data_list = sim_data, n_draws = 2000,
                                  ci_level = 0.90)

compare_kf_wz(kf_result = cf_kf, wz_result = cf_wz)
```

This produces:
- An overlay plot showing both methods' means and confidence bands
- A table of discrepancies (horizon × variable)

**When to use which:**

| Criterion | Kalman Filter | Waggoner & Zha |
|---|---|---|
| Constraint type | Hard or soft (via measurement noise) | Hard (exact) or soft (distributional) |
| Uncertainty source | State covariance propagation | Monte Carlo draws of structural shocks |
| Constraint satisfaction | In expectation (soft by construction) | Exact per draw (hard mode) |
| Best for | Mixing frequencies, sequential updating | Scenario analysis, structural interpretation |
| Computation | Fast (analytic) | Moderate (Monte Carlo, but parallelisable) |

---

## Conditional Forecasting: Kalman Filter

The Kalman filter approach (module `08_kalman_filter.R`) treats conditions as noisy observations of the latent state in state-space form. See the code comments for the full state-space formulation.

```r
conditions <- data.frame(
  variable = c("US.gdp", "US.gdp", "EU.inf"),
  horizon  = c(1, 2, 4),
  value    = c(0.5, 0.5, -0.2),
  hardness = c(1e-8, 1e-8, 1e-8)   # hard constraints
)

cf <- conditional_forecast(
  gvar_model = gvar_model,
  conditions = conditions,
  max_h      = 8,
  data_list  = sim_data,
  ci_level   = 0.90
)
```

---

## Dependencies

The toolkit will auto-install these CRAN packages if missing:

- `vars`, `urca`, `tseries` — time-series modelling
- `MASS`, `Matrix` — matrix algebra
- `ggplot2`, `reshape2` — plotting
- `stats` — base R (Kalman filter helpers)

---

## Key References

- Pesaran, Schuermann & Weiner (2004). "Modeling Regional Interdependencies Using a Global Error-Correcting Macroeconometric Model." *Journal of Business & Economic Statistics*.
- Dees, di Mauro, Pesaran & Smith (2007). "Exploring the International Linkages of the Euro Area: A Global VAR Analysis." *Journal of Applied Econometrics*.
- Pesaran & Shin (1998). "Generalized Impulse Response Analysis in Linear Multivariate Models." *Economics Letters*.
- Waggoner & Zha (1999). "Conditional Forecasts in Dynamic Multivariate Models." *Review of Economics and Statistics*.
- Doan, Litterman & Sims (1984). "Forecasting and Conditional Projection Using Realistic Prior Distributions." *Econometric Reviews*.
