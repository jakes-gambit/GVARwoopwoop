###############################################################################
#  00_fetch_fred_data.R  –  Fetch Quarterly Macro Data from FRED for GVAR
#
#  Variables: Real GDP (growth), CPI Inflation, Short-term Rate, 10Y Yield,
#             Real Effective Exchange Rate (REER), Oil, Equity Prices
#  Frequency: Quarterly
#  Coverage:  25+ countries (G7, Euro area, CEE, Emerging markets)
#
#  PREREQUISITES:
#  1. Get free FRED API key: https://fred.stlouisfed.org/docs/api/api_key.html
#  2. Install packages: install.packages(c("fredr", "dplyr", "tidyr"))
#
###############################################################################

# ═══════════════════════════════════════════════════════════════════════════════
# 0.  Setup
# ═══════════════════════════════════════════════════════════════════════════════

library(fredr)
library(dplyr)
library(tidyr)

# ──────────────────────────────────────────────────────────────────────────────
# SET YOUR FRED API KEY HERE
# Get free key at: https://fred.stlouisfed.org/docs/api/api_key.html
# ──────────────────────────────────────────────────────────────────────────────

FRED_API_KEY <- "ee414e1326bd3ad7ad5d4bae231901b0"   # <-- PASTE YOUR KEY

fredr_set_key(FRED_API_KEY)

# ═══════════════════════════════════════════════════════════════════════════════
# 1.  Define FRED Series IDs by Country
# ═══════════════════════════════════════════════════════════════════════════════

# These are OECD Main Economic Indicators series hosted on FRED
# GDP:     Real GDP index (compute growth)
# CPI:     Consumer price index (compute inflation)
# Rate:    3-month interbank rate
# LT_Rate: 10-year government bond yield (OECD MEI: IRLTLT01xxQ156N)
# REER:    Real effective exchange rate, CPI-based (OECD MEI: CCRETT01xxQ661N)
# EQ:      Share price index, all shares (OECD MEI: SPASTT01xxQ661N) – [NEW]
#          Returns NA for countries without OECD MEI coverage; handled gracefully
start_date = "1970-01-01"
end_date = "2019-12-31"
fred_series <- list(

  # ══════════════════════════════════════════════════════════════════════════
  # AGGREGATES
  # ══════════════════════════════════════════════════════════════════════════
  # Euro Area - using HICP index
  EA  = list(gdp = "CLVMNACSCAB1GQEA19", cpi = "CP0000EZ19M086NEST", rate = "IR3TIB01EZQ156N",
             lt_rate = "IRLTLT01EZQ156N",  reer = "CCRETT01EZQ661N",
             eq = "SPASTT01EZQ661N"),                                        # [NEW] EA equity index

  # ══════════════════════════════════════════════════════════════════════════
  # G7
  # ══════════════════════════════════════════════════════════════════════════
  USA = list(gdp = "GDPC1",              cpi = "CPIAUCSL",          rate = "TB3MS",
             lt_rate = "GS10",             reer = "CCRETT01USQ661N",
             eq = "SPASTT01USQ661N"),                                        # [NEW] US equity index
  GBR = list(gdp = "CLVMNACSCAB1GQUK",   cpi = "GBRCPIALLMINMEI",   rate = "IR3TIB01GBQ156N",
             lt_rate = "IRLTLT01GBQ156N",  reer = "CCRETT01GBQ661N",
             eq = "SPASTT01GBQ661N"),                                        # [NEW] UK equity index
  DEU = list(gdp = "CLVMNACSCAB1GQDE",   cpi = "DEUCPIALLMINMEI",   rate = "IR3TIB01DEQ156N",
             lt_rate = "IRLTLT01DEQ156N",  reer = "CCRETT01DEQ661N",
             eq = "SPASTT01DEQ661N"),                                        # [NEW] Germany equity index
  FRA = list(gdp = "CLVMNACSCAB1GQFR",   cpi = "FRACPIALLMINMEI",   rate = "IR3TIB01FRQ156N",
             lt_rate = "IRLTLT01FRQ156N",  reer = "CCRETT01FRQ661N",
             eq = "SPASTT01FRQ661N"),                                        # [NEW] France equity index
  ITA = list(gdp = "CLVMNACSCAB1GQIT",   cpi = "ITACPIALLMINMEI",   rate = "IR3TIB01ITQ156N",
             lt_rate = "IRLTLT01ITQ156N",  reer = "CCRETT01ITQ661N",
             eq = "SPASTT01ITQ661N"),                                        # [NEW] Italy equity index
  JPN = list(gdp = "JPNRGDPEXP",         cpi = "JPNCPIALLMINMEI",   rate = "IR3TIB01JPQ156N",
             lt_rate = "IRLTLT01JPQ156N",  reer = "CCRETT01JPQ661N",
             eq = "SPASTT01JPQ661N"),                                        # [NEW] Japan equity index
  CAN = list(gdp = "NGDPRSAXDCCAQ",      cpi = "CANCPIALLMINMEI",   rate = "IR3TIB01CAQ156N",
             lt_rate = "IRLTLT01CAQ156N",  reer = "CCRETT01CAQ661N",
             eq = "SPASTT01CAQ661N"),                                        # [NEW] Canada equity index

  # ══════════════════════════════════════════════════════════════════════════
  # Other Advanced
  # ══════════════════════════════════════════════════════════════════════════
  AUS = list(gdp = "CLVMNACSCAB1GQAU",   cpi = "AUSCPIALLMINMEI",   rate = "IR3TIB01AUQ156N",
             lt_rate = "IRLTLT01AUQ156N",  reer = "CCRETT01AUQ661N",
             eq = "SPASTT01AUQ661N"),                                        # [NEW] Australia equity index
  CHE = list(gdp = "CLVMNACSCAB1GQCH",   cpi = "CHECPIALLMINMEI",   rate = "IR3TIB01CHQ156N",
             lt_rate = "IRLTLT01CHQ156N",  reer = "CCRETT01CHQ661N",
             eq = "SPASTT01CHQ661N"),                                        # [NEW] Switzerland equity index
  SWE = list(gdp = "CLVMNACSCAB1GQSE",   cpi = "SWECPIALLMINMEI",   rate = "IR3TIB01SEQ156N",
             lt_rate = "IRLTLT01SEQ156N",  reer = "CCRETT01SEQ661N",
             eq = "SPASTT01SEQ661N"),                                        # [NEW] Sweden equity index
  NOR = list(gdp = "CLVMNACSCAB1GQNO",   cpi = "NORCPIALLMINMEI",   rate = "IR3TIB01NOQ156N",
             lt_rate = "IRLTLT01NOQ156N",  reer = "CCRETT01NOQ661N",
             eq = "SPASTT01NOQ661N"),                                        # [NEW] Norway equity index
  KOR = list(gdp = "CLVMNACSCAB1GQKR",   cpi = "KORCPIALLMINMEI",   rate = "IR3TIB01KRQ156N",
             lt_rate = "IRLTLT01KRQ156N",  reer = "CCRETT01KRQ661N",
             eq = "SPASTT01KRQ661N"),                                        # [NEW] Korea equity index
  DNK = list(gdp = "CLVMNACSCAB1GQDK",   cpi = "DNKCPIALLMINMEI",   rate = "IR3TIB01DKQ156N",
             lt_rate = "IRLTLT01DKQ156N",  reer = "CCRETT01DKQ661N",
             eq = "SPASTT01DKQ661N"),                                        # [NEW] Denmark equity index

  # ══════════════════════════════════════════════════════════════════════════
  # Euro Area Countries
  # ══════════════════════════════════════════════════════════════════════════
  NLD = list(gdp = "CLVMNACSCAB1GQNL",   cpi = "NLDCPIALLMINMEI",   rate = "IR3TIB01NLQ156N",
             lt_rate = "IRLTLT01NLQ156N",  reer = "CCRETT01NLQ661N",
             eq = "SPASTT01NLQ661N"),                                        # [NEW] Netherlands equity index
  BEL = list(gdp = "CLVMNACSCAB1GQBE",   cpi = "BELCPIALLMINMEI",   rate = "IR3TIB01BEQ156N",
             lt_rate = "IRLTLT01BEQ156N",  reer = "CCRETT01BEQ661N",
             eq = "SPASTT01BEQ661N"),                                        # [NEW] Belgium equity index
  AUT = list(gdp = "CLVMNACSCAB1GQAT",   cpi = "AUTCPIALLMINMEI",   rate = "IR3TIB01ATQ156N",
             lt_rate = "IRLTLT01ATQ156N",  reer = "CCRETT01ATQ661N",
             eq = "SPASTT01ATQ661N"),                                        # [NEW] Austria equity index
  ESP = list(gdp = "CLVMNACSCAB1GQES",   cpi = "ESPCPIALLMINMEI",   rate = "IR3TIB01ESQ156N",
             lt_rate = "IRLTLT01ESQ156N",  reer = "CCRETT01ESQ661N",
             eq = "SPASTT01ESQ661N"),                                        # [NEW] Spain equity index
  PRT = list(gdp = "CLVMNACSCAB1GQPT",   cpi = "PRTCPIALLMINMEI",   rate = "IR3TIB01PTQ156N",
             lt_rate = "IRLTLT01PTQ156N",  reer = "CCRETT01PTQ661N",
             eq = "SPASTT01PTQ661N"),                                        # [NEW] Portugal equity index
  GRC = list(gdp = "CLVMNACSCAB1GQGR",   cpi = "GRCCPIALLMINMEI",   rate = "IR3TIB01GRQ156N",
             lt_rate = "IRLTLT01GRQ156N",  reer = "CCRETT01GRQ661N",
             eq = "SPASTT01GRQ661N"),                                        # [NEW] Greece equity index
  IRL = list(gdp = "CLVMNACSCAB1GQIE",   cpi = "IRLCPIALLMINMEI",   rate = "IR3TIB01IEQ156N",
             lt_rate = "IRLTLT01IEQ156N",  reer = "CCRETT01IEQ661N",
             eq = "SPASTT01IEQ661N"),                                        # [NEW] Ireland equity index
  FIN = list(gdp = "CLVMNACSCAB1GQFI",   cpi = "FINCPIALLMINMEI",   rate = "IR3TIB01FIQ156N",
             lt_rate = "IRLTLT01FIQ156N",  reer = "CCRETT01FIQ661N",
             eq = "SPASTT01FIQ661N"),                                        # [NEW] Finland equity index

  # ══════════════════════════════════════════════════════════════════════════
  # CEE
  # ══════════════════════════════════════════════════════════════════════════
  POL = list(gdp = "CLVMNACSCAB1GQPL",   cpi = "POLCPIALLMINMEI",   rate = "IR3TIB01PLQ156N",
             lt_rate = "IRLTLT01PLQ156N",  reer = "CCRETT01PLQ661N",
             eq = "SPASTT01PLQ661N"),                                        # [NEW] Poland equity index
  CZE = list(gdp = "CLVMNACSCAB1GQCZ",   cpi = "CZECPIALLMINMEI",   rate = "IR3TIB01CZQ156N",
             lt_rate = "IRLTLT01CZQ156N",  reer = "CCRETT01CZQ661N",
             eq = "SPASTT01CZQ661N"),                                        # [NEW] Czech equity index
  HUN = list(gdp = "CLVMNACSCAB1GQHU",   cpi = "HUNCPIALLMINMEI",   rate = "IR3TIB01HUQ156N",
             lt_rate = "IRLTLT01HUQ156N",  reer = "CCRETT01HUQ661N",
             eq = "SPASTT01HUQ661N"),                                        # [NEW] Hungary equity index
  SVK = list(gdp = "CLVMNACSCAB1GQSK",   cpi = "SVKCPIALLMINMEI",   rate = "IR3TIB01SKQ156N",
             lt_rate = "IRLTLT01SKQ156N",  reer = "CCRETT01SKQ661N"),        # no OECD equity series
  SVN = list(gdp = "CLVMNACSCAB1GQSI",   cpi = "SVNCPIALLMINMEI",   rate = "IR3TIB01SIQ156N",
             lt_rate = "IRLTLT01SIQ156N",  reer = "CCRETT01SIQ661N"),        # no OECD equity series
  ROU = list(gdp = "CLVMNACSCAB1GQRO",   cpi = "ROUCPIALLMINMEI",   rate = "IR3TIB01ROQ156N",
             lt_rate = "IRLTLT01ROQ156N",  reer = "CCRETT01ROQ661N"),        # no OECD equity series
  BGR = list(gdp = "CLVMNACSCAB1GQBG",   cpi = "BGRCPIALLMINMEI",   rate = "IR3TIB01BGQ156N",
             lt_rate = "IRLTLT01BGQ156N",  reer = "CCRETT01BGQ661N"),        # no OECD equity series
  HRV = list(gdp = "CLVMNACSCAB1GQHR",   cpi = "HRVCPIALLMINMEI",   rate = "IR3TIB01HRQ156N",
             lt_rate = "IRLTLT01HRQ156N",  reer = "CCRETT01HRQ661N"),        # no OECD equity series

  # ══════════════════════════════════════════════════════════════════════════
  # Emerging Markets
  # ══════════════════════════════════════════════════════════════════════════
  MEX = list(gdp = "CLVMNACSCAB1GQMX",   cpi = "MEXCPIALLMINMEI",   rate = "IR3TIB01MXQ156N",
             lt_rate = "IRLTLT01MXQ156N",  reer = "CCRETT01MXQ661N",
             eq = "SPASTT01MXQ661N"),                                        # [NEW] Mexico equity index
  BRA = list(gdp = "CLVMNACSCAB1GQBR",   cpi = "BRACPIALLMINMEI",   rate = "IR3TIB01BRQ156N",
             lt_rate = "IRLTLT01BRQ156N",  reer = "CCRETT01BRQ661N",
             eq = "SPASTT01BRQ661N"),                                        # [NEW] Brazil equity index
  ZAF = list(gdp = "CLVMNACSCAB1GQZA",   cpi = "ZAFCPIALLMINMEI",   rate = "IR3TIB01ZAQ156N",
             lt_rate = "IRLTLT01ZAQ156N",  reer = "CCRETT01ZAQ661N",
             eq = "SPASTT01ZAQ661N"),                                        # [NEW] S. Africa equity index
  TUR = list(gdp = "CLVMNACSCAB1GQTR",   cpi = "TURCPIALLMINMEI",   rate = "IR3TIB01TRQ156N",
             lt_rate = "IRLTLT01TRQ156N",  reer = "CCRETT01TRQ661N",
             eq = "SPASTT01TRQ661N"),                                        # [NEW] Turkey equity index
  RUS = list(gdp = "CLVMNACSCAB1GQRU",   cpi = "RUSCPIALLMINMEI",   rate = "IR3TIB01RUQ156N",
             lt_rate = "IRLTLT01RUQ156N",  reer = "CCRETT01RUQ661N"),        # no OECD equity series
  IND = list(gdp = "CLVMNACSCAB1GQIN",   cpi = "INDCPIALLMINMEI",   rate = "IR3TIB01INQ156N",
             lt_rate = "IRLTLT01INQ156N",  reer = "CCRETT01INQ661N",
             eq = "SPASTT01INQ661N"),                                        # [NEW] India equity index
  CHN = list(gdp = "NAEXKP01CNQ657S",    cpi = "CHNCPIALLMINMEI",   rate = "IR3TIB01CNQ156N",
             lt_rate = "IRLTLT01CNQ156N",  reer = "CCRETT01CNQ661N"),        # no OECD equity series
  IDN = list(gdp = "CLVMNACSCAB1GQID",   cpi = "IDNCPIALLMINMEI",   rate = "IR3TIB01IDQ156N",
             lt_rate = "IRLTLT01IDQ156N",  reer = "CCRETT01IDQ661N",
             eq = "SPASTT01IDQ661N")                                         # [NEW] Indonesia equity index
)

# Oil price (Brent crude) - Global variable
OIL_SERIES <- "DCOILBRENTEU"

# ═══════════════════════════════════════════════════════════════════════════════
# 2.  Data Fetching Functions
# ═══════════════════════════════════════════════════════════════════════════════

#' Fetch a single FRED series
fetch_series <- function(series_id, start = start_date, end = end_date) {
  
  tryCatch({
    df <- fredr(
      series_id = series_id,
      observation_start = as.Date(start),
      observation_end = as.Date(end),
      frequency = "q"
    )
    
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    # Validate expected columns exist (guards against fredr version differences
    # and namespace conflicts where dplyr::select is masked)
    if (!all(c("date", "value") %in% names(df))) {
      cat(sprintf("[WARN] Series %s missing expected columns. Got: %s\n",
                  series_id, paste(names(df), collapse = ", ")))
      return(NULL)
    }
    
    return(df[, c("date", "value")])
    
  }, error = function(e) {
    cat(sprintf("[WARN] Failed to fetch %s: %s\n", series_id, conditionMessage(e)))
    return(NULL)
  })
  
  return(NULL)
}

#' Fetch all variables for one country
fetch_country <- function(country_code, series_list, start = start_date, end = end_date) {
  
  cat(sprintf("  Fetching %-3s ... ", country_code))
  
  results <- list()
  
  for (var in names(series_list)) {
    df <- fetch_series(series_list[[var]], start, end)
    
    if (!is.null(df) && "value" %in% names(df)) {
      names(df)[names(df) == "value"] <- var
      results[[var]] <- df
    }
    
    Sys.sleep(0.15)  # Rate limiting
  }
  
  if (length(results) == 0) {
    cat("FAILED\n")
    return(NULL)
  }
  
  # Merge all variables by date
  merged <- results[[1]]
  if (length(results) > 1) {
    for (i in 2:length(results)) {
      merged <- merge(merged, results[[i]], by = "date", all = TRUE)
    }
  }
  
  merged$country <- country_code
  
  n_complete <- sum(complete.cases(merged))
  cat(sprintf("OK (%d obs, %d complete)\n", nrow(merged), n_complete))
  
  return(merged)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 3.  Fetch All Data
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("  FRED Data Download for GVAR\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

cat("Fetching country data...\n\n")

all_data <- list()

for (ctry in names(fred_series)) {
  df <- fetch_country(ctry, fred_series[[ctry]])
  if (!is.null(df)) {
    all_data[[ctry]] <- df
  }
}

# Combine all countries
fred_raw <- dplyr::bind_rows(all_data)

cat("\n")
cat(sprintf("Raw data: %d observations, %d countries\n",
            nrow(fred_raw), length(unique(fred_raw$country))))

# ═══════════════════════════════════════════════════════════════════════════════
# 4.  Fetch Oil Price (Global Variable)
# ═══════════════════════════════════════════════════════════════════════════════

cat("\nFetching oil price (Brent crude)... ")

oil_raw <- fetch_series(OIL_SERIES)

if (!is.null(oil_raw)) {
  oil_data <- oil_raw %>%
    dplyr::arrange(date) %>%
    dplyr::mutate(
      oil          = (value / dplyr::lag(value, 4) - 1) * 100,   # YoY % growth (default)
      oil_level    = value,                                         # [NEW] price level ($/barrel)
      oil_log      = log(value),                                    # [NEW] log of price level
      oil_logdiff  = c(NA, diff(log(value)))                        # [NEW] quarterly log-change
    ) %>%
    dplyr::select(date, oil, oil_level, oil_log, oil_logdiff)

  cat(sprintf("OK (%d obs)\n", nrow(oil_data)))
} else {
  cat("FAILED\n")
  oil_data <- NULL
}

# ═══════════════════════════════════════════════════════════════════════════════
# 5.  Transform to Growth Rates / Inflation
# ═══════════════════════════════════════════════════════════════════════════════

cat("\nComputing growth rates and inflation...\n")

fred_transformed <- fred_raw %>%
  dplyr::group_by(country) %>%
  dplyr::arrange(date) %>%
  dplyr::mutate(
    # Real GDP growth (YoY %)
    gdp_growth = (gdp / dplyr::lag(gdp, 4) - 1) * 100,

    # CPI inflation (YoY %)
    inflation = (cpi / dplyr::lag(cpi, 4) - 1) * 100,

    # REER: YoY % change (index level fetched, compute growth)
    reer_growth = (reer / dplyr::lag(reer, 4) - 1) * 100,

    # [NEW] Equity price: YoY % change (log-return approximation)
    # Countries without an eq series will have NA here
    eq_growth = dplyr::if_else(
      !is.na(eq),
      (eq / dplyr::lag(eq, 4) - 1) * 100,
      NA_real_
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(country, date, gdp = gdp_growth, inf = inflation,
                rate, lt_rate, reer = reer_growth,
                eq = eq_growth)    # [NEW] equity price growth

# Add oil price
if (!is.null(oil_data)) {
  fred_transformed <- fred_transformed %>%
    dplyr::left_join(oil_data, by = "date")
}

# Remove rows with NA in key variables (first 4 quarters due to YoY computation)
fred_final <- fred_transformed %>%
  dplyr::filter(!is.na(gdp), !is.na(inf)) %>%
  dplyr::arrange(country, date)

cat(sprintf("Transformed data: %d observations\n", nrow(fred_final)))

# ═══════════════════════════════════════════════════════════════════════════════
# 6.  Check Data Coverage
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("  Data Coverage Summary\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

coverage <- fred_final %>%
  dplyr::group_by(country) %>%
  dplyr::summarise(
    n_obs = dplyr::n(),
    first_date = min(date),
    last_date = max(date),
    gdp_na = sum(is.na(gdp)),
    inf_na = sum(is.na(inf)),
    rate_na = sum(is.na(rate)),
    lt_rate_na = sum(is.na(lt_rate)),
    reer_na = sum(is.na(reer)),
    oil_na = sum(is.na(oil)),
    eq_na = sum(is.na(eq)),    # [NEW] equity coverage gap
    complete = sum(complete.cases(gdp, inf, rate)),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(complete))

print(coverage, n = 40)

# ═══════════════════════════════════════════════════════════════════════════════
# 7.  Convert to GVAR List Format
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("  Converting to GVAR Format\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

#' Convert data frame to list of matrices (GVAR format)
#'
#' @param df         Data frame with country, date, and variable columns
#' @param var_cols   Character vector of variable column names
#' @param min_obs    Minimum complete observations required (default 60)
#' @return           Named list of T x k matrices
to_gvar_list <- function(df, var_cols = c("gdp", "inf", "rate"), min_obs = 60) {
  
  countries <- unique(df$country)
  gvar_list <- list()
  
  for (ctry in countries) {
    ctry_df <- df %>%
      dplyr::filter(country == ctry) %>%
      dplyr::arrange(date)

    # Check complete cases for selected variables
    complete_mask <- complete.cases(ctry_df[, var_cols, drop = FALSE])
    n_complete <- sum(complete_mask)
    
    if (n_complete < min_obs) {
      cat(sprintf("  Skipping %s: only %d complete obs (need %d)\n",
                  ctry, n_complete, min_obs))
      next
    }
    
    # Filter to complete cases
    ctry_df <- ctry_df[complete_mask, ]
    
    # Create matrix
    mat <- ctry_df %>%
      dplyr::select(all_of(var_cols)) %>%
      as.matrix()
    
    # Format row names as quarters (e.g., "2000-Q1")
    rownames(mat) <- paste0(
      format(ctry_df$date, "%Y"),
      "-Q",
      ceiling(as.numeric(format(ctry_df$date, "%m")) / 3)
    )
    
    gvar_list[[ctry]] <- mat
    cat(sprintf("  %s: %d obs, %s to %s\n",
                ctry, nrow(mat), rownames(mat)[1], rownames(mat)[nrow(mat)]))
  }
  
  return(gvar_list)
}

# Create GVAR dataset (without oil - oil is global/exogenous)
gvar_data <- to_gvar_list(fred_final, var_cols = c("gdp", "inf", "rate", "lt_rate", "reer"), min_obs = 60)

# Also create version with oil (YoY % growth, default)
gvar_data_with_oil <- to_gvar_list(fred_final, var_cols = c("gdp", "inf", "rate", "lt_rate", "reer", "oil"), min_obs = 60)

# [NEW] Version with equity prices added
gvar_data_with_eq <- to_gvar_list(fred_final, var_cols = c("gdp", "inf", "rate", "lt_rate", "reer", "eq"), min_obs = 60)

# [NEW] Full dataset: oil + equity prices
gvar_data_full <- to_gvar_list(fred_final, var_cols = c("gdp", "inf", "rate", "lt_rate", "reer", "oil", "eq"), min_obs = 60)

cat("\n")
cat(sprintf("GVAR dataset created: %d countries\n", length(gvar_data)))
cat(sprintf("Variables: %s\n", paste(colnames(gvar_data[[1]]), collapse = ", ")))
cat(sprintf("Sample period: %s to %s\n",
            rownames(gvar_data[[1]])[1],
            rownames(gvar_data[[1]])[nrow(gvar_data[[1]])]))


# ═══════════════════════════════════════════════════════════════════════════════
# 8.  Save Data
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("  Saving Data\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

# Save as RData
save(gvar_data, gvar_data_with_oil, gvar_data_with_eq, gvar_data_full,  # [NEW] added eq variants
     fred_final, oil_data,
     file = "gvar_fred_data.RData")

cat("Data saved to: gvar_fred_data.RData\n")
cat("\nObjects saved:\n")
cat("  - gvar_data:          List of country matrices (gdp, inf, rate, lt_rate, reer)\n")
cat("  - gvar_data_with_oil: List of country matrices (gdp, inf, rate, lt_rate, reer, oil)\n")
cat("  - gvar_data_with_eq:  List of country matrices (gdp, inf, rate, lt_rate, reer, eq)\n")  # [NEW]
cat("  - gvar_data_full:     List of country matrices (gdp, inf, rate, lt_rate, reer, oil, eq)\n")  # [NEW]
cat("  - fred_final:         Full data frame\n")
cat("  - oil_data:           Oil price series (oil, oil_level, oil_log, oil_logdiff)\n")  # [NEW] updated

# ═══════════════════════════════════════════════════════════════════════════════
# 9.  Preview Data
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("  Data Preview\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

cat("First 3 countries (first 5 rows each):\n\n")

for (ctry in names(gvar_data)[1:min(3, length(gvar_data))]) {
  cat(paste0("── ", ctry, " ──\n"))
  print(head(gvar_data[[ctry]], 5))
  cat("\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# 10. Plot All Series by Country
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("  Plotting Series by Country\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

#' Plot all variables for each country
#'
#' @param df         Data frame with country, date, and variable columns
#' @param countries  Character vector of countries to plot (NULL = all)
#' @param save_pdf   If TRUE, save to PDF file
#' @param pdf_file   PDF filename
plot_country_series <- function(df, countries = NULL, save_pdf = FALSE,
                                pdf_file = "fred_data_plots.pdf") {
  
  if (is.null(countries)) {
    countries <- unique(df$country)
  }
  
  # Variables to plot
  var_cols <- c("gdp", "inf", "rate", "lt_rate", "reer", "oil", "eq")  # [NEW] eq added
  var_cols <- var_cols[var_cols %in% names(df)]
  var_labels <- c(gdp = "GDP Growth (%)", inf = "Inflation (%)",
                  rate = "Short-term Rate (%)", lt_rate = "10Y Yield (%)",
                  reer = "REER Growth (%)", oil = "Oil Price Growth (%)",
                  eq = "Equity Price Growth (%)")                              # [NEW]
  var_colors <- c(gdp = "steelblue", inf = "darkred",
                  rate = "darkgreen", lt_rate = "purple",
                  reer = "darkcyan", oil = "darkorange",
                  eq = "brown")                                                 # [NEW]
  
  if (save_pdf) {
    pdf(pdf_file, width = 10, height = 8)
  }
  
  for (ctry in countries) {
    ctry_df <- df %>%
      dplyr::filter(country == ctry) %>%
      dplyr::arrange(date)

    if (nrow(ctry_df) < 5) next
    
    # Set up plot layout (3x2 for up to 6 variables)
    n_vars <- length(var_cols)
    n_rows <- ceiling(n_vars / 2)
    par(mfrow = c(n_rows, 2), mar = c(4, 4, 2, 1), oma = c(0, 0, 3, 0))
    
    for (v in var_cols) {
      y <- ctry_df[[v]]
      x <- ctry_df$date
      
      # Skip if all NA
      if (all(is.na(y))) {
        plot.new()
        text(0.5, 0.5, paste(var_labels[v], "\n(No Data)"), cex = 1.2)
        next
      }
      
      # Calculate ylim with buffer
      ylim <- range(y, na.rm = TRUE)
      ylim <- ylim + c(-0.1, 0.1) * diff(ylim)
      
      # Plot
      plot(x, y, type = "l", col = var_colors[v], lwd = 2,
           xlab = "Date", ylab = var_labels[v],
           main = var_labels[v], ylim = ylim)
      
      # Add zero line for growth rates
      if (v %in% c("gdp", "inf", "reer", "oil", "eq")) {  # [NEW] eq added
        abline(h = 0, col = "gray", lty = 2)
      }
      
      # Add grid
      grid(col = "lightgray", lty = 3)
    }
    
    # Add country title
    mtext(ctry, outer = TRUE, cex = 1.5, font = 2)
  }
  
  if (save_pdf) {
    dev.off()
    cat(sprintf("Plots saved to: %s\n", pdf_file))
  }
  
  invisible(NULL)
}

# Generate plots
plot_country_series(fred_final, save_pdf = TRUE, pdf_file = "fred_data_plots.pdf")

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("  DONE - Data ready for GVAR estimation\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("\n")
cat("Next steps:\n")
cat("  1. Load data: load('gvar_fred_data.RData')\n")
cat("  2. Create your weight matrix (bilateral trade flows)\n")
cat("  3. Run unit root tests: ur_results <- panel_unit_root(gvar_data)\n")
cat("  4. Estimate GVAR: Follow main.R workflow\n")
cat("\n")
cat("Note: Content generated using AI - expert verification recommended.\n")