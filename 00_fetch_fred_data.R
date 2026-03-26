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
end_date   <- "2023-12-31"

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
    oil_logdiff = oil_log - dplyr::lag(oil_log)
  )

# Helper: compute log-diff safely within a vector
# log(NA) = NA propagates naturally; no pmax needed
ld <- function(x) c(NA_real_, diff(log(x)))

# Ensure eq column exists for all countries (absent for some)
if (!"eq" %in% names(raw)) raw$eq <- NA_real_

fred_data <- raw %>%
  dplyr::group_by(country) %>%
  dplyr::arrange(date) %>%
  dplyr::transmute(
    date,
    gdp_level   = gdp,
    gdp_log     = log(gdp),
    gdp_logdiff = ld(gdp),
    cpi_level   = cpi,
    cpi_log     = log(cpi),
    cpi_logdiff = ld(cpi),
    rate,
    lt_rate,
    reer_level   = reer,
    reer_log     = log(reer),
    reer_logdiff = ld(reer),
    eq_level     = eq,
    eq_log       = log(eq),
    eq_logdiff   = ld(eq)
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
