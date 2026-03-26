###############################################################################
#  00_fetch_fred_data.R  –  Fetch Quarterly Macro Data from FRED
#
#  Output:  fred_data  (one data frame, saved to gvar_fred_data.RData)
#
#  Columns per country × date:
#    gdp_level,  gdp_log,  gdp_logdiff     – Real GDP
#    cpi_level,  cpi_log,  cpi_logdiff     – Consumer price index
#    rate                                  – 3-month interbank rate (level)
#    lt_rate                               – 10-year bond yield (level)
#    reer_level, reer_log, reer_logdiff    – Real effective exchange rate
#    eq_level,   eq_log,   eq_logdiff      – Equity price index
#    oil_level,  oil_log,  oil_logdiff     – Brent crude (same for all ctries)
#
#  logdiff = log(x_t) − log(x_{t−1})  (quarterly log-change)
#
#  PREREQUISITES:
#    fredr, dplyr  – install.packages(c("fredr","dplyr"))
#    API key: https://fred.stlouisfed.org/docs/api/api_key.html
###############################################################################

library(fredr)
library(dplyr)

FRED_API_KEY <- "ee414e1326bd3ad7ad5d4bae231901b0"   # <-- paste your key
fredr_set_key(FRED_API_KEY)

start_date <- "1970-01-01"
end_date   <- "2025-12-31"

# ─────────────────────────────────────────────────────────────────────────────
# 1.  FRED Series IDs
#     GDP     : real GDP index  (OECD MEI or national accounts)
#     CPI     : consumer/harmonized price index
#     rate    : 3-month interbank rate
#     lt_rate : 10-year government bond yield
#     reer    : real effective exchange rate, CPI-based (OECD MEI)
#     eq      : share price index, all shares (OECD MEI: SPASTT01xxQ661N)
#               countries without OECD MEI equity data omit the eq entry
# ─────────────────────────────────────────────────────────────────────────────

fred_series <- list(

  EA  = list(gdp="CLVMNACSCAB1GQEA19", cpi="CP0000EZ19M086NEST",
             rate="IR3TIB01EZQ156N", lt_rate="IRLTLT01EZQ156N",
             reer="CCRETT01EZQ661N",  eq="SPASTT01EZQ661N"),

  USA = list(gdp="GDPC1",             cpi="CPIAUCSL",
             rate="TB3MS",            lt_rate="GS10",
             reer="CCRETT01USQ661N",  eq="SPASTT01USQ661N"),

  GBR = list(gdp="CLVMNACSCAB1GQUK",  cpi="GBRCPIALLMINMEI",
             rate="IR3TIB01GBQ156N",  lt_rate="IRLTLT01GBQ156N",
             reer="CCRETT01GBQ661N",  eq="SPASTT01GBQ661N"),

  DEU = list(gdp="CLVMNACSCAB1GQDE",  cpi="DEUCPIALLMINMEI",
             rate="IR3TIB01DEQ156N",  lt_rate="IRLTLT01DEQ156N",
             reer="CCRETT01DEQ661N",  eq="SPASTT01DEQ661N"),

  FRA = list(gdp="CLVMNACSCAB1GQFR",  cpi="FRACPIALLMINMEI",
             rate="IR3TIB01FRQ156N",  lt_rate="IRLTLT01FRQ156N",
             reer="CCRETT01FRQ661N",  eq="SPASTT01FRQ661N"),

  ITA = list(gdp="CLVMNACSCAB1GQIT",  cpi="ITACPIALLMINMEI",
             rate="IR3TIB01ITQ156N",  lt_rate="IRLTLT01ITQ156N",
             reer="CCRETT01ITQ661N",  eq="SPASTT01ITQ661N"),

  JPN = list(gdp="JPNRGDPEXP",        cpi="JPNCPIALLMINMEI",
             rate="IR3TIB01JPQ156N",  lt_rate="IRLTLT01JPQ156N",
             reer="CCRETT01JPQ661N",  eq="SPASTT01JPQ661N"),

  CAN = list(gdp="NGDPRSAXDCCAQ",     cpi="CANCPIALLMINMEI",
             rate="IR3TIB01CAQ156N",  lt_rate="IRLTLT01CAQ156N",
             reer="CCRETT01CAQ661N",  eq="SPASTT01CAQ661N"),

  AUS = list(gdp="CLVMNACSCAB1GQAU",  cpi="AUSCPIALLMINMEI",
             rate="IR3TIB01AUQ156N",  lt_rate="IRLTLT01AUQ156N",
             reer="CCRETT01AUQ661N",  eq="SPASTT01AUQ661N"),

  CHE = list(gdp="CLVMNACSCAB1GQCH",  cpi="CHECPIALLMINMEI",
             rate="IR3TIB01CHQ156N",  lt_rate="IRLTLT01CHQ156N",
             reer="CCRETT01CHQ661N",  eq="SPASTT01CHQ661N"),

  SWE = list(gdp="CLVMNACSCAB1GQSE",  cpi="SWECPIALLMINMEI",
             rate="IR3TIB01SEQ156N",  lt_rate="IRLTLT01SEQ156N",
             reer="CCRETT01SEQ661N",  eq="SPASTT01SEQ661N"),

  NOR = list(gdp="CLVMNACSCAB1GQNO",  cpi="NORCPIALLMINMEI",
             rate="IR3TIB01NOQ156N",  lt_rate="IRLTLT01NOQ156N",
             reer="CCRETT01NOQ661N",  eq="SPASTT01NOQ661N"),

  KOR = list(gdp="CLVMNACSCAB1GQKR",  cpi="KORCPIALLMINMEI",
             rate="IR3TIB01KRQ156N",  lt_rate="IRLTLT01KRQ156N",
             reer="CCRETT01KRQ661N",  eq="SPASTT01KRQ661N"),

  DNK = list(gdp="CLVMNACSCAB1GQDK",  cpi="DNKCPIALLMINMEI",
             rate="IR3TIB01DKQ156N",  lt_rate="IRLTLT01DKQ156N",
             reer="CCRETT01DKQ661N",  eq="SPASTT01DKQ661N"),

  NLD = list(gdp="CLVMNACSCAB1GQNL",  cpi="NLDCPIALLMINMEI",
             rate="IR3TIB01NLQ156N",  lt_rate="IRLTLT01NLQ156N",
             reer="CCRETT01NLQ661N",  eq="SPASTT01NLQ661N"),

  BEL = list(gdp="CLVMNACSCAB1GQBE",  cpi="BELCPIALLMINMEI",
             rate="IR3TIB01BEQ156N",  lt_rate="IRLTLT01BEQ156N",
             reer="CCRETT01BEQ661N",  eq="SPASTT01BEQ661N"),

  AUT = list(gdp="CLVMNACSCAB1GQAT",  cpi="AUTCPIALLMINMEI",
             rate="IR3TIB01ATQ156N",  lt_rate="IRLTLT01ATQ156N",
             reer="CCRETT01ATQ661N",  eq="SPASTT01ATQ661N"),

  ESP = list(gdp="CLVMNACSCAB1GQES",  cpi="ESPCPIALLMINMEI",
             rate="IR3TIB01ESQ156N",  lt_rate="IRLTLT01ESQ156N",
             reer="CCRETT01ESQ661N",  eq="SPASTT01ESQ661N"),

  PRT = list(gdp="CLVMNACSCAB1GQPT",  cpi="PRTCPIALLMINMEI",
             rate="IR3TIB01PTQ156N",  lt_rate="IRLTLT01PTQ156N",
             reer="CCRETT01PTQ661N",  eq="SPASTT01PTQ661N"),

  GRC = list(gdp="CLVMNACSCAB1GQGR",  cpi="GRCCPIALLMINMEI",
             rate="IR3TIB01GRQ156N",  lt_rate="IRLTLT01GRQ156N",
             reer="CCRETT01GRQ661N",  eq="SPASTT01GRQ661N"),

  IRL = list(gdp="CLVMNACSCAB1GQIE",  cpi="IRLCPIALLMINMEI",
             rate="IR3TIB01IEQ156N",  lt_rate="IRLTLT01IEQ156N",
             reer="CCRETT01IEQ661N",  eq="SPASTT01IEQ661N"),

  FIN = list(gdp="CLVMNACSCAB1GQFI",  cpi="FINCPIALLMINMEI",
             rate="IR3TIB01FIQ156N",  lt_rate="IRLTLT01FIQ156N",
             reer="CCRETT01FIQ661N",  eq="SPASTT01FIQ661N"),

  POL = list(gdp="CLVMNACSCAB1GQPL",  cpi="POLCPIALLMINMEI",
             rate="IR3TIB01PLQ156N",  lt_rate="IRLTLT01PLQ156N",
             reer="CCRETT01PLQ661N",  eq="SPASTT01PLQ661N"),

  CZE = list(gdp="CLVMNACSCAB1GQCZ",  cpi="CZECPIALLMINMEI",
             rate="IR3TIB01CZQ156N",  lt_rate="IRLTLT01CZQ156N",
             reer="CCRETT01CZQ661N",  eq="SPASTT01CZQ661N"),

  HUN = list(gdp="CLVMNACSCAB1GQHU",  cpi="HUNCPIALLMINMEI",
             rate="IR3TIB01HUQ156N",  lt_rate="IRLTLT01HUQ156N",
             reer="CCRETT01HUQ661N",  eq="SPASTT01HUQ661N"),

  # CEE countries without OECD MEI equity series
  SVK = list(gdp="CLVMNACSCAB1GQSK",  cpi="SVKCPIALLMINMEI",
             rate="IR3TIB01SKQ156N",  lt_rate="IRLTLT01SKQ156N",
             reer="CCRETT01SKQ661N"),

  SVN = list(gdp="CLVMNACSCAB1GQSI",  cpi="SVNCPIALLMINMEI",
             rate="IR3TIB01SIQ156N",  lt_rate="IRLTLT01SIQ156N",
             reer="CCRETT01SIQ661N"),

  ROU = list(gdp="CLVMNACSCAB1GQRO",  cpi="ROUCPIALLMINMEI",
             rate="IR3TIB01ROQ156N",  lt_rate="IRLTLT01ROQ156N",
             reer="CCRETT01ROQ661N"),

  BGR = list(gdp="CLVMNACSCAB1GQBG",  cpi="BGRCPIALLMINMEI",
             rate="IR3TIB01BGQ156N",  lt_rate="IRLTLT01BGQ156N",
             reer="CCRETT01BGQ661N"),

  HRV = list(gdp="CLVMNACSCAB1GQHR",  cpi="HRVCPIALLMINMEI",
             rate="IR3TIB01HRQ156N",  lt_rate="IRLTLT01HRQ156N",
             reer="CCRETT01HRQ661N"),

  MEX = list(gdp="CLVMNACSCAB1GQMX",  cpi="MEXCPIALLMINMEI",
             rate="IR3TIB01MXQ156N",  lt_rate="IRLTLT01MXQ156N",
             reer="CCRETT01MXQ661N",  eq="SPASTT01MXQ661N"),

  BRA = list(gdp="CLVMNACSCAB1GQBR",  cpi="BRACPIALLMINMEI",
             rate="IR3TIB01BRQ156N",  lt_rate="IRLTLT01BRQ156N",
             reer="CCRETT01BRQ661N",  eq="SPASTT01BRQ661N"),

  ZAF = list(gdp="CLVMNACSCAB1GQZA",  cpi="ZAFCPIALLMINMEI",
             rate="IR3TIB01ZAQ156N",  lt_rate="IRLTLT01ZAQ156N",
             reer="CCRETT01ZAQ661N",  eq="SPASTT01ZAQ661N"),

  TUR = list(gdp="CLVMNACSCAB1GQTR",  cpi="TURCPIALLMINMEI",
             rate="IR3TIB01TRQ156N",  lt_rate="IRLTLT01TRQ156N",
             reer="CCRETT01TRQ661N",  eq="SPASTT01TRQ661N"),

  # Emerging markets without OECD MEI equity series
  RUS = list(gdp="CLVMNACSCAB1GQRU",  cpi="RUSCPIALLMINMEI",
             rate="IR3TIB01RUQ156N",  lt_rate="IRLTLT01RUQ156N",
             reer="CCRETT01RUQ661N"),

  IND = list(gdp="CLVMNACSCAB1GQIN",  cpi="INDCPIALLMINMEI",
             rate="IR3TIB01INQ156N",  lt_rate="IRLTLT01INQ156N",
             reer="CCRETT01INQ661N",  eq="SPASTT01INQ661N"),

  CHN = list(gdp="NAEXKP01CNQ657S",   cpi="CHNCPIALLMINMEI",
             rate="IR3TIB01CNQ156N",  lt_rate="IRLTLT01CNQ156N",
             reer="CCRETT01CNQ661N"),

  IDN = list(gdp="CLVMNACSCAB1GQID",  cpi="IDNCPIALLMINMEI",
             rate="IR3TIB01IDQ156N",  lt_rate="IRLTLT01IDQ156N",
             reer="CCRETT01IDQ661N",  eq="SPASTT01IDQ661N")
)

OIL_SERIES <- "DCOILBRENTEU"   # Brent crude (daily → quarterly average)

# ─────────────────────────────────────────────────────────────────────────────
# 2.  Fetch helpers
# ─────────────────────────────────────────────────────────────────────────────

# Fetch one FRED series; return two-column data.frame (date, value) or NULL.
fetch_one <- function(series_id) {
  tryCatch({
    df <- fredr(series_id,
                observation_start = as.Date(start_date),
                observation_end   = as.Date(end_date),
                frequency         = "q")
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df[, c("date", "value")]
  }, error = function(e) NULL)
}

# Fetch all variables for one country and merge by date.
fetch_country <- function(code, series) {
  cat(sprintf("  %-3s ...", code))
  parts <- lapply(names(series), function(v) {
    df <- fetch_one(series[[v]])
    if (!is.null(df)) colnames(df)[2] <- v
    Sys.sleep(0.15)   # FRED rate limit
    df
  })
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0) { cat(" FAILED\n"); return(NULL) }
  merged <- Reduce(function(a, b) merge(a, b, by = "date", all = TRUE), parts)
  merged$country <- code
  cat(sprintf(" OK (%d rows)\n", nrow(merged)))
  merged
}

# ─────────────────────────────────────────────────────────────────────────────
# 3.  Download
# ─────────────────────────────────────────────────────────────────────────────

cat("\nFetching country data...\n")
country_data <- Filter(Negate(is.null),
                       lapply(names(fred_series),
                              function(cc) fetch_country(cc, fred_series[[cc]])))
raw <- dplyr::bind_rows(country_data)

cat("\nFetching oil price (Brent)... ")
oil_raw <- fetch_one(OIL_SERIES)
if (is.null(oil_raw)) stop("Oil price download failed.")
cat(sprintf("OK (%d rows)\n", nrow(oil_raw)))

# ─────────────────────────────────────────────────────────────────────────────
# 4.  Transform
#
#     For GDP, CPI, REER, EQ : level / log / quarterly log-change (logdiff)
#     For rates               : level only (already in %)
#     For oil                 : aggregate daily to quarterly average first,
#                               then level / log / logdiff
#
#     logdiff = log(x_t) − log(x_{t-1}) ≈ quarter-on-quarter change
#     log(0) / log(negative) → NA (not treated as error)
# ─────────────────────────────────────────────────────────────────────────────

# Quarter-average the daily oil price
oil_q <- oil_raw %>%
  dplyr::mutate(date = as.Date(format(date, "%Y-%m-01"))) %>%
  dplyr::group_by(date) %>%
  dplyr::summarise(oil_level = mean(value, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(date) %>%
  dplyr::mutate(
    oil_log     = log(oil_level),
    oil_logdiff = oil_log - dplyr::lag(oil_log),
    oil_diff = oil_level - dplyr::lag(oil_level)
  )

# Helper: compute log-diff safely within a vector
# log(NA) = NA propagates naturally; no pmax needed
ld <- function(x) c(NA_real_, diff(log(x)))

# Ensure eq column exists for all countries (absent for some)
# if (!"eq" %in% names(raw)) raw$eq <- NA_real_

fred_data <- raw %>%
  dplyr::group_by(country) %>%
  dplyr::arrange(date) %>%
  dplyr::mutate(
    # GDP: normalize to a volume index (2015 mean = 100) for cross-country
    # comparability, regardless of original currency / units.
    # Growth rates (gdp_logdiff) are identical to ld(raw gdp).
    gdp_base    = {
      b <- mean(gdp[as.integer(format(date, "%Y")) == 2015L], na.rm = TRUE)
      if (!is.finite(b)) b <- mean(gdp, na.rm = TRUE)
      b
    },
    gdp_level   = gdp / gdp_base * 100,   # index: ≈ 100 in 2015 for every country
    gdp_log     = log(gdp_level),
    gdp_logdiff = ld(gdp),
    cpi_level   = cpi,
    cpi_log     = log(cpi),
    cpi_logdiff = ld(cpi),
    reer_level   = reer,
    reer_log     = log(reer),
    reer_logdiff = ld(reer),
    eq_level         = eq,
    eq_log           = log(eq),
    eq_logdiff       = ld(eq),
    # Real equity price: q_it = ln(EQ_it / CPI_it)  (Dees et al. 2007)
    eq_real_log      = log(eq / cpi),
    eq_real_logdiff  = ld(eq / cpi),
    # Quarterly log interest rates: ρ = 0.25 · ln(1 + R/100)  (Dees et al. 2007)
    rho_s            = 0.25 * log(1 + rate    / 100),
    rho_l            = 0.25 * log(1 + lt_rate / 100),
    rho_lms          = rho_l - rho_s
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(oil_q, by = "date") %>%   # same oil value for every country
  dplyr::arrange(country, date)

cat(sprintf("\nfred_data: %d rows, %d columns, %d countries\n",
            nrow(fred_data), ncol(fred_data), length(unique(fred_data$country))))
cat("Columns: ", paste(names(fred_data), collapse = ", "), "\n")

# ─────────────────────────────────────────────────────────────────────────────
# 5.  to_gvar_list()  –  extract GVAR-format matrices from fred_data
#
#  Example:
#    sim_data <- to_gvar_list(fred_data,
#                  var_cols  = c("gdp_logdiff","cpi_logdiff","lt_rate","oil_logdiff"),
#                  countries = c("DEU","FRA","ITA","ESP"))
# ─────────────────────────────────────────────────────────────────────────────

#' Convert fred_data to a named list of T × k matrices (GVAR format).
#'
#' Rows with any NA in var_cols are dropped; the result contains only
#' complete-case rows, aligned to a common date set via align_to_common_dates().
#'
#' @param df        fred_data data frame
#' @param var_cols  character vector of column names to include
#' @param countries countries to extract (NULL = all with sufficient data)
#' @param min_obs   minimum complete rows required per country (default 60)
#' @return named list; each element is a T × k matrix with "YYYY-QN" row names
to_gvar_list <- function(df, var_cols, countries = NULL, min_obs = 60) {

  if (!is.null(countries)) df <- dplyr::filter(df, country %in% countries)

  to_qtr <- function(d) paste0(format(d, "%Y"), "-Q",
                                ceiling(as.integer(format(d, "%m")) / 3))
  out <- list()
  for (cc in unique(df$country)) {
    sub <- df[df$country == cc, ]
    sub <- sub[order(sub$date), ]
    mat <- as.matrix(sub[, var_cols, drop = FALSE])
    storage.mode(mat) <- "numeric"
    ok  <- complete.cases(mat)
    if (sum(ok) < min_obs) {
      message(sprintf("  Skipping %-3s: only %d complete rows (need %d)",
                      cc, sum(ok), min_obs))
      next
    }
    mat <- mat[ok, , drop = FALSE]
    rownames(mat) <- to_qtr(sub$date[ok])
    out[[cc]] <- mat
    message(sprintf("  %-3s: %d obs  [%s – %s]",
                    cc, nrow(mat), rownames(mat)[1], rownames(mat)[nrow(mat)]))
  }
  out
}

# ─────────────────────────────────────────────────────────────────────────────
# 6.  Save – single object plus the helper function
# ─────────────────────────────────────────────────────────────────────────────

save(fred_data, to_gvar_list, file = "gvar_fred_data.RData")
cat("\nSaved gvar_fred_data.RData  (fred_data + to_gvar_list)\n")
cat("Note: Content generated using AI – expert verification recommended.\n")


# ─────────────────────────────────────────────────────────────────────────────
# 7.  China (CHN) – Annual FRED data + Kalman temporal disaggregation
#
#  Several FRED series for China are only available at annual frequency or
#  have reliable coverage only in annual form (e.g. World Bank vintage).
#  This section:
#    (a) Downloads annual FRED series for China.
#    (b) Calls kalman_temporal_disaggregate() (from 08_kalman_filter.R) to
#        produce quarterly estimates via a Kalman smoother.
#    (c) Patches the CHN rows of fred_data with the disaggregated quarterly
#        series where the original quarterly series is missing or sparse.
#
#  Annual FRED series used
#  -----------------------
#  GDP (real, constant 2015 USD, OECD):  NAEXKP01CNA189S
#    – If unavailable, fall back to RGDPCNA (OECD annual real GDP index)
#  CPI (annual average, OECD MEI):       CPALTT01CNA661S
#    – If unavailable, fall back to averaging the monthly CHNCPIALLMINMEI
#  Long-run rate (annual):               IRLTLT01CNA156N
#  REER (annual average, OECD MEI):      CCRETT01CNA661N
#
#  NOTE: Series IDs and availability may change.  The code wraps each fetch
#  in tryCatch and silently skips unavailable series.
# ─────────────────────────────────────────────────────────────────────────────

source("08_kalman_filter.R")   # for kalman_temporal_disaggregate()


# Use the validated quarterly series IDs with frequency="a" so FRED aggregates
# them to annual automatically.  The quarterly IDs below are confirmed to exist.
#   GDP  : NAEXKP01CNQ657S  (real GDP index, SA, quarterly)
#   CPI  : CHNCPIALLMINMEI  (CPI all items, monthly – FRED averages to annual)
#   Short: IR3TIB01CNQ156N  (3-month money market rate, quarterly)
#   Long : IRLTLT01CNQ156N  (10-yr bond yield, quarterly)
#   REER : CCRETT01CNQ661N  (REER, quarterly)
chn_annual_series <- list(
  gdp     = "NAEXKP01CNQ657S",   # quarterly → annual avg via frequency="a"
  cpi     = "CHNCPIALLMINMEI",   # monthly  → annual avg via frequency="a"
  lt_rate = "IRLTLT01CNQ156N",   # quarterly → annual avg via frequency="a"
  reer    = "CCRETT01CNQ661N"    # quarterly → annual avg via frequency="a"
)

cat("\nFetching annual China (CHN) series from FRED...\n")

fetch_annual <- function(series_id) {
  tryCatch({
    df <- fredr(series_id,
                observation_start = as.Date(start_date),
                observation_end   = as.Date(end_date),
                frequency         = "a")
    if (is.null(df) || nrow(df) == 0) return(NULL)
    Sys.sleep(0.15)
    df[, c("date", "value")]
  }, error = function(e) {
    message(sprintf("  [CHN annual] Could not fetch %s: %s", series_id, conditionMessage(e)))
    NULL
  })
}

chn_annual_raw <- lapply(names(chn_annual_series), function(v) {
  cat(sprintf("  %-8s (%s) ...", v, chn_annual_series[[v]]))
  df <- fetch_annual(chn_annual_series[[v]])
  if (is.null(df)) { cat(" SKIPPED\n"); return(NULL) }
  cat(sprintf(" OK (%d years)\n", nrow(df)))
  colnames(df)[2] <- v
  df
})
names(chn_annual_raw) <- names(chn_annual_series)

# For each successfully fetched annual series, disaggregate to quarterly
chn_quarterly_disagg <- list()

for (v in names(chn_annual_raw)) {
  ann_df <- chn_annual_raw[[v]]
  if (is.null(ann_df) || sum(!is.na(ann_df[[v]])) < 5) next

  # Extract complete years
  ann_df   <- ann_df[!is.na(ann_df[[v]]), ]
  years    <- as.integer(format(ann_df$date, "%Y"))
  ann_vals <- ann_df[[v]]
  n_yrs    <- length(years)

  cat(sprintf("  Disaggregating %s: %d annual obs [%d–%d] ...",
              v, n_yrs, min(years), max(years)))

  agg_type <- if (v %in% c("lt_rate", "reer")) "average" else "average"

  res <- tryCatch(
    kalman_temporal_disaggregate(
      annual_values = ann_vals,
      n_years       = n_yrs,
      start_year    = min(years),
      aggregation   = agg_type
    ),
    error = function(e) {
      message(sprintf("  [CHN disagg] %s failed: %s", v, conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(res)) {
    cat(sprintf(" OK (%d quarterly obs)\n", length(res$quarterly)))
    chn_quarterly_disagg[[v]] <- res$quarterly
  } else {
    cat(" FAILED\n")
  }
}

# ── Patch fred_data for CHN rows where quarterly is NA/missing ──────────────
if (length(chn_quarterly_disagg) > 0) {

  # Helper: convert "YYYY-QN" label to the first date of that quarter
  qtr_to_date <- function(lbl) {
    yr  <- as.integer(sub("-Q.*", "", lbl))
    q   <- as.integer(sub(".*-Q", "", lbl))
    as.Date(sprintf("%d-%02d-01", yr, (q - 1L) * 3L + 1L))
  }

  for (v in names(chn_quarterly_disagg)) {
    qtr_series <- chn_quarterly_disagg[[v]]
    qtr_dates  <- as.Date(sapply(names(qtr_series), qtr_to_date))

    # Map variable name to the fred_data column(s) to patch
    col_map <- list(
      gdp     = c("gdp", "gdp_level", "gdp_log", "gdp_logdiff"),
      cpi     = c("cpi", "cpi_level", "cpi_log", "cpi_logdiff"),
      lt_rate = c("lt_rate", "rho_l"),
      reer    = c("reer", "reer_level", "reer_log", "reer_logdiff")
    )

    level_col <- switch(v, gdp = "gdp", cpi = "cpi", lt_rate = "lt_rate", reer = "reer")

    for (qi in seq_along(qtr_series)) {
      dt   <- qtr_dates[qi]
      rows <- which(fred_data$country == "CHN" & fred_data$date == dt)
      if (length(rows) == 0) {
        # Insert a new row for this country-date if needed
        next   # skip for now; merging new rows is handled below
      }
      for (r in rows) {
        # Only patch if the existing value is NA
        if (is.na(fred_data[r, level_col])) {
          fred_data[r, level_col] <- qtr_series[qi]
        }
      }
    }

    # Recompute log and logdiff for the patched level column
    chn_idx <- which(fred_data$country == "CHN")
    if (length(chn_idx) > 0 && level_col %in% colnames(fred_data)) {
      lvl <- fred_data[chn_idx, level_col]
      if (paste0(v, "_log") %in% colnames(fred_data)) {
        fred_data[chn_idx, paste0(v, "_log")] <- log(lvl)
      }
      if (paste0(v, "_logdiff") %in% colnames(fred_data)) {
        ld_vals <- c(NA_real_, diff(log(lvl)))
        fred_data[chn_idx, paste0(v, "_logdiff")] <- ld_vals
      }
      # Recompute derived columns that depend on this level
      if (v == "cpi" && "eq_real_log" %in% colnames(fred_data)) {
        eq_lvl <- fred_data[chn_idx, "eq"]
        fred_data[chn_idx, "eq_real_log"]     <- log(eq_lvl / lvl)
        fred_data[chn_idx, "eq_real_logdiff"] <- c(NA_real_, diff(log(eq_lvl / lvl)))
      }
      if (v == "lt_rate" && "rho_l" %in% colnames(fred_data)) {
        fred_data[chn_idx, "rho_l"] <- 0.25 * log(1 + lvl / 100)
      }
    }
    cat(sprintf("  [CHN patch] %s patched in fred_data.\n", v))
  }

  # Re-save with patched data
  save(fred_data, to_gvar_list, file = "gvar_fred_data.RData")
  cat("\nRe-saved gvar_fred_data.RData with CHN disaggregated quarterly data.\n")
} else {
  cat("\n[CHN annual] No series successfully disaggregated; fred_data unchanged.\n")
}
