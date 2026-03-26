###############################################################################
#  00b_fetch_weight_data.R  –  Bilateral Trade Weight Matrix from Empirical Data
#
#  Builds a row-normalised N×N weight matrix W for any set of countries.
#  Data sources, tried in order:
#    1. OECD Bilateral Trade in Goods (OECD package, OECD.Stat API)
#    2. IMF Direction of Trade Statistics (DOTS, JSON REST API – no key needed)
#
#  Main entry point (called from main.R):
#    W_raw <- fetch_trade_weights(countries, avg_years = 2014:2016)
#
#  The matrix can optionally be blended with BIS banking-claims weights to
#  capture financial linkages alongside trade linkages.
#
#  PREREQUISITES:
#    dplyr, httr, jsonlite  – always needed
#    OECD                   – install.packages("OECD")  [optional, preferred]
#    BIS                    – install.packages("BIS")   [optional, blending]
###############################################################################

library(dplyr)

# ─────────────────────────────────────────────────────────────────────────────
# ISO-3 ↔ IMF-area-code mapping
# The IMF DOTS API uses its own 2-letter area codes (mostly ISO-2, a few differ).
# ─────────────────────────────────────────────────────────────────────────────

ISO3_TO_IMF2 <- c(
  USA="US", GBR="GB", DEU="DE", FRA="FR", ITA="IT", JPN="JP", CAN="CA",
  AUS="AU", CHE="CH", SWE="SE", NOR="NO", KOR="KR", DNK="DK",
  NLD="NL", BEL="BE", AUT="AT", ESP="ES", PRT="PT", GRC="GR",
  IRL="IE", FIN="FI", POL="PL", CZE="CZ", HUN="HU", SVK="SK",
  SVN="SI", ROU="RO", BGR="BG", HRV="HR", MEX="MX", BRA="BR",
  ZAF="ZA", TUR="TR", RUS="RU", IND="IN", CHN="CN", IDN="ID",
  EA ="U2"
)

# Row-normalise a weight matrix (zero diagonal, each row sums to 1).
.row_norm <- function(W) {
  diag(W) <- 0
  rs <- rowSums(W, na.rm = TRUE)
  rs[rs == 0] <- 1   # isolated country → equal weights after norm
  W / rs
}

# Expand a smaller matrix to the full country set, filling missing pairs with 0.
.pad <- function(W, countries) {
  N   <- length(countries)
  out <- matrix(0, N, N, dimnames = list(countries, countries))
  hit <- intersect(rownames(W), countries)
  hit_c <- intersect(colnames(W), countries)
  if (length(hit) > 0 && length(hit_c) > 0)
    out[hit, hit_c] <- W[hit, hit_c, drop = FALSE]
  out
}

# ─────────────────────────────────────────────────────────────────────────────
# Source 1 : OECD Bilateral Trade in Goods
#   Dataset : TRADE_IN_GOODS  (annual, USD millions)
#   Requires: install.packages("OECD")
# ─────────────────────────────────────────────────────────────────────────────

.oecd_trade <- function(countries, avg_years) {

  if (!requireNamespace("OECD", quietly = TRUE)) return(NULL)

  cat("[Weights] Fetching bilateral trade from OECD...\n")

  raw <- tryCatch(
    OECD::get_dataset(
      dataset    = "TRADE_IN_GOODS",
      filter     = list(REPORTER = countries,
                        PARTNER  = countries,
                        FLOW     = "EXP",
                        MEASURE  = "USD"),
      start_time = min(avg_years),
      end_time   = max(avg_years)
    ),
    error = function(e) { message("  OECD error: ", e$message); NULL }
  )

  if (is.null(raw) || nrow(raw) == 0) return(NULL)

  # OECD package may return "obsValue" or "Value" depending on version
  val_col <- intersect(c("obsValue", "Value", "value"), names(raw))[1]
  if (is.na(val_col)) return(NULL)

  bilateral <- raw %>%
    dplyr::filter(REPORTER %in% countries, PARTNER %in% countries,
                  REPORTER != PARTNER) %>%
    dplyr::group_by(REPORTER, PARTNER) %>%
    dplyr::summarise(trade = mean(as.numeric(.data[[val_col]]), na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::filter(is.finite(trade), trade > 0)

  if (nrow(bilateral) == 0) return(NULL)

  W <- tidyr::pivot_wider(bilateral, names_from = PARTNER, values_from = trade,
                           values_fill = 0) %>%
    tibble::column_to_rownames("REPORTER") %>%
    as.matrix()

  W <- .pad(W, countries)
  covered <- mean(rowSums(W) > 0)
  cat(sprintf("  OECD: %.0f%% of countries covered.\n", covered * 100))
  .row_norm(W)
}

# ─────────────────────────────────────────────────────────────────────────────
# Source 2 : IMF Direction of Trade Statistics (DOTS)
#   Endpoint: https://data.imf.org/api/SDMX_JSON/CompactData/DOT/
#   No API key required.  Uses annual exports (TXG_FOB_USD).
# ─────────────────────────────────────────────────────────────────────────────

.imf_dots <- function(countries, avg_years) {

  if (!requireNamespace("httr",     quietly = TRUE) ||
      !requireNamespace("jsonlite", quietly = TRUE)) {
    message("  [Weights] httr/jsonlite not installed – skipping IMF DOTS.")
    return(NULL)
  }

  cat("[Weights] Fetching bilateral trade from IMF DOTS...\n")

  # Map ISO-3 → IMF-2; warn about unmapped codes
  imf2 <- ISO3_TO_IMF2[countries]
  missing <- countries[is.na(imf2)]
  if (length(missing) > 0)
    message("  [Weights] No IMF-2 code for: ", paste(missing, collapse=", "),
            " – these will get zero bilateral weight.")
  imf2  <- imf2[!is.na(imf2)]
  if (length(imf2) == 0) return(NULL)

  partners_str <- paste(imf2, collapse = "+")
  N     <- length(countries)
  W_imf <- matrix(0, length(imf2), length(imf2),
                   dimnames = list(imf2, imf2))

  for (rep in imf2) {
    url <- sprintf(
      "https://data.imf.org/api/SDMX_JSON/CompactData/DOT/A.%s.TXG_FOB_USD.%s?startPeriod=%d&endPeriod=%d",
      rep, partners_str, min(avg_years), max(avg_years)
    )

    resp <- tryCatch(
      httr::GET(url, httr::timeout(30)),
      error = function(e) NULL
    )
    if (is.null(resp) || httr::status_code(resp) != 200) {
      message(sprintf("  [Weights] IMF DOTS failed for reporter %s", rep))
      Sys.sleep(1)
      next
    }

    parsed <- tryCatch(
      jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"),
                         simplifyVector = FALSE),
      error = function(e) NULL
    )
    if (is.null(parsed)) next

    # Navigate SDMX-JSON structure: dataSets[[1]]$series
    series <- tryCatch(parsed$dataSets[[1]]$series, error = function(e) NULL)
    if (is.null(series)) next

    # Each key is "freq:reporter:indicator:partner" (0-indexed positions)
    # We need the partner dimension positions from structure
    dims <- tryCatch(
      parsed$structure$dimensions$series,
      error = function(e) NULL
    )
    if (is.null(dims)) next

    # Find partner dimension index (last dimension in our query)
    partner_dim <- dims[[length(dims)]]$values   # list of {id, name}
    partner_ids <- sapply(partner_dim, `[[`, "id")

    obs_dim <- parsed$structure$dimensions$observation
    time_ids <- sapply(obs_dim[[1]]$values, `[[`, "id")

    for (key in names(series)) {
      key_parts <- as.integer(strsplit(key, ":")[[1]]) + 1L  # 0-indexed → 1-indexed
      par_idx   <- key_parts[length(key_parts)]
      partner   <- partner_ids[par_idx]

      if (is.na(partner) || !(partner %in% imf2) || partner == rep) next

      obs <- series[[key]]$observations
      vals <- sapply(names(obs), function(t_idx) {
        yr <- time_ids[as.integer(t_idx) + 1L]
        if (is.null(yr)) return(NA_real_)
        if (as.integer(yr) %in% avg_years) as.numeric(obs[[t_idx]][[1]]) else NA_real_
      })
      vals <- vals[!is.na(vals) & is.finite(vals)]
      if (length(vals) > 0)
        W_imf[rep, partner] <- mean(vals)
    }

    Sys.sleep(0.5)   # polite rate limiting
  }

  # Rename back from IMF-2 to ISO-3
  iso3_present <- names(ISO3_TO_IMF2)[ISO3_TO_IMF2 %in% imf2]
  rownames(W_imf) <- names(ISO3_TO_IMF2)[match(rownames(W_imf), ISO3_TO_IMF2)]
  colnames(W_imf) <- names(ISO3_TO_IMF2)[match(colnames(W_imf), ISO3_TO_IMF2)]

  W_full <- .pad(W_imf, countries)
  covered <- mean(rowSums(W_full) > 0)
  cat(sprintf("  IMF DOTS: %.0f%% of countries covered.\n", covered * 100))
  .row_norm(W_full)
}

# ─────────────────────────────────────────────────────────────────────────────
# Source 3 : BIS Consolidated Banking Statistics (optional blending)
#   Dataset : CBS immediate-counterparty basis (Table B1)
#   Requires: install.packages("BIS")
# ─────────────────────────────────────────────────────────────────────────────

.bis_finance <- function(countries, avg_years) {

  if (!requireNamespace("BIS", quietly = TRUE)) return(NULL)

  cat("[Weights] Fetching financial linkages from BIS CBS...\n")

  raw <- tryCatch(
    BIS::get_bis("CBS"),
    error = function(e) { message("  BIS error: ", e$message); NULL }
  )
  if (is.null(raw) || nrow(raw) == 0) return(NULL)

  # Identify column names (BIS package column names can vary across versions)
  rep_col <- grep("reporting", names(raw), ignore.case = TRUE, value = TRUE)[1]
  ctp_col <- grep("counterpart", names(raw), ignore.case = TRUE, value = TRUE)[1]
  val_col <- grep("obs_value|value", names(raw), ignore.case = TRUE, value = TRUE)[1]
  dat_col <- grep("^date$|time", names(raw), ignore.case = TRUE, value = TRUE)[1]

  if (any(is.na(c(rep_col, ctp_col, val_col, dat_col)))) {
    message("  BIS: unexpected column layout – skipping.")
    return(NULL)
  }

  bilateral <- raw %>%
    dplyr::filter(
      .data[[rep_col]] %in% countries,
      .data[[ctp_col]] %in% countries,
      .data[[rep_col]] != .data[[ctp_col]],
      as.integer(format(as.Date(.data[[dat_col]]), "%Y")) %in% avg_years
    ) %>%
    dplyr::group_by(reporter = .data[[rep_col]], counterparty = .data[[ctp_col]]) %>%
    dplyr::summarise(claims = mean(as.numeric(.data[[val_col]]), na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::filter(is.finite(claims), claims > 0)

  if (nrow(bilateral) == 0) return(NULL)

  W <- tidyr::pivot_wider(bilateral, names_from = counterparty, values_from = claims,
                           values_fill = 0) %>%
    tibble::column_to_rownames("reporter") %>%
    as.matrix()

  W_full <- .pad(W, countries)
  cat(sprintf("  BIS: %.0f%% of countries covered.\n",
              mean(rowSums(W_full) > 0) * 100))
  .row_norm(W_full)
}

# ─────────────────────────────────────────────────────────────────────────────
# Source 4 : World Bank GDP (fallback when bilateral trade data unavailable)
#   API: https://api.worldbank.org/v2/country/{iso2}/indicator/NY.GDP.MKTP.CD
#   No API key required.  Constructs GDP-proportional weights: each country
#   gives weight proportional to partner GDP (rough proxy for trade links).
# ─────────────────────────────────────────────────────────────────────────────

# ISO-3 → ISO-2 mapping for World Bank API
ISO3_TO_ISO2 <- c(
  USA="US", GBR="GB", DEU="DE", FRA="FR", ITA="IT", JPN="JP", CAN="CA",
  AUS="AU", CHE="CH", SWE="SE", NOR="NO", KOR="KR", DNK="DK",
  NLD="NL", BEL="BE", AUT="AT", ESP="ES", PRT="PT", GRC="GR",
  IRL="IE", FIN="FI", POL="PL", CZE="CZ", HUN="HU", SVK="SK",
  SVN="SI", ROU="RO", BGR="BG", HRV="HR", MEX="MX", BRA="BR",
  ZAF="ZA", TUR="TR", RUS="RU", IND="IN", CHN="CN", IDN="ID",
  EA ="XC"
)

.wb_gdp_weights <- function(countries, avg_years) {

  if (!requireNamespace("httr",     quietly = TRUE) ||
      !requireNamespace("jsonlite", quietly = TRUE)) {
    message("  [Weights] httr/jsonlite not installed – skipping World Bank GDP.")
    return(NULL)
  }

  cat("[Weights] Fetching GDP from World Bank API (GDP-proportional fallback)...\n")

  iso2 <- ISO3_TO_ISO2[countries]
  missing <- countries[is.na(iso2)]
  if (length(missing) > 0)
    message("  [Weights] No ISO-2 code for: ", paste(missing, collapse=", "),
            " – these will get zero weight.")
  iso2_valid <- iso2[!is.na(iso2)]
  if (length(iso2_valid) == 0) return(NULL)

  gdp_vals <- setNames(rep(NA_real_, length(countries)), countries)

  for (i in seq_along(iso2_valid)) {
    cc3   <- names(iso2_valid)[i]
    cc2   <- iso2_valid[i]
    url   <- sprintf(
      "https://api.worldbank.org/v2/country/%s/indicator/NY.GDP.MKTP.CD?format=json&date=%d:%d&per_page=100",
      cc2, min(avg_years), max(avg_years)
    )
    resp <- tryCatch(httr::GET(url, httr::timeout(20)), error = function(e) NULL)
    if (is.null(resp) || httr::status_code(resp) != 200) {
      message(sprintf("  [Weights] WB GDP failed for %s", cc3))
      Sys.sleep(0.5)
      next
    }
    parsed <- tryCatch(
      jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"),
                         simplifyVector = TRUE),
      error = function(e) NULL
    )
    if (is.null(parsed) || length(parsed) < 2) next
    data_df <- parsed[[2]]
    if (is.null(data_df) || nrow(data_df) == 0) next
    vals <- suppressWarnings(as.numeric(data_df$value))
    vals <- vals[!is.na(vals) & is.finite(vals)]
    if (length(vals) > 0)
      gdp_vals[cc3] <- mean(vals)
    Sys.sleep(0.3)
  }

  gdp_known <- gdp_vals[!is.na(gdp_vals) & gdp_vals > 0]
  if (length(gdp_known) == 0) {
    message("  [Weights] No GDP data retrieved from World Bank.")
    return(NULL)
  }

  # GDP-proportional weight: w_ij = gdp_j / sum_{k≠i} gdp_k
  N   <- length(countries)
  W   <- matrix(0, N, N, dimnames = list(countries, countries))
  for (i in seq_len(N)) {
    for (j in seq_len(N)) {
      if (i != j && !is.na(gdp_vals[countries[j]]) && gdp_vals[countries[j]] > 0)
        W[i, j] <- gdp_vals[countries[j]]
    }
  }

  covered <- mean(rowSums(W) > 0)
  cat(sprintf("  WB GDP: %.0f%% of countries covered.\n", covered * 100))
  .row_norm(W)
}

# ─────────────────────────────────────────────────────────────────────────────
# Source 5 : Hardcoded approximate GDP shares (last-resort fallback)
#   Uses approximate 2024 nominal GDP (USD trillion) from IMF WEO Oct-2024.
#   Countries not in the table get the median share among knowns.
# ─────────────────────────────────────────────────────────────────────────────

.APPROX_GDP_2024 <- c(
  USA=28.8, CHN=18.5, DEU=4.5,  JPN=4.1,  IND=3.9,  GBR=3.1,  FRA=3.0,
  ITA=2.3,  BRA=2.2,  CAN=2.2,  RUS=2.2,  KOR=1.8,  MEX=1.8,  AUS=1.7,
  ESP=1.6,  IDN=1.4,  NLD=1.1,  TUR=1.1,  CHE=0.9,  IRL=0.6,  NOR=0.6,
  POL=0.8,  SWE=0.6,  BEL=0.6,  AUT=0.5,  ZAF=0.4,  DNK=0.4,  ROU=0.4,
  FIN=0.3,  CZE=0.3,  GRC=0.3,  PRT=0.3,  HUN=0.2,  BGR=0.1,  HRV=0.08,
  SVK=0.12, SVN=0.06, EA=16.0
)

.hardcoded_gdp_weights <- function(countries) {
  cat("[Weights] Using hardcoded approximate GDP shares (last resort).\n")
  N   <- length(countries)
  gdp <- .APPROX_GDP_2024[countries]
  # Countries missing from table get median share
  gdp[is.na(gdp)] <- median(.APPROX_GDP_2024, na.rm = TRUE)
  names(gdp) <- countries

  W <- matrix(0, N, N, dimnames = list(countries, countries))
  for (i in seq_len(N)) {
    for (j in seq_len(N)) {
      if (i != j) W[i, j] <- gdp[countries[j]]
    }
  }
  .row_norm(W)
}

# ─────────────────────────────────────────────────────────────────────────────
# Main entry point
# ─────────────────────────────────────────────────────────────────────────────

#' Fetch empirical bilateral trade data and build a GVAR weight matrix.
#'
#' Tries data sources in order: OECD → IMF DOTS → World Bank GDP →
#' hardcoded approximate GDP shares.  Optionally blends with BIS
#' financial-claims weights (set finance_alpha > 0 to activate).
#'
#' @param countries     ISO-3 character vector matching names(sim_data)
#' @param avg_years     years to average (default 2014–2016 pre-COVID window)
#' @param finance_alpha share of BIS financial weights in the blend (0 = trade
#'                      only, 0.5 = equal trade/finance blend).  Requires the
#'                      BIS package.
#' @return Named N×N row-normalised weight matrix
fetch_trade_weights <- function(countries,
                                avg_years     = 2014:2016,
                                finance_alpha = 0) {

  cat("══════════════════════════════════════════════════════════════\n")
  cat(sprintf("  Building weight matrix for: %s\n",
              paste(countries, collapse=", ")))
  cat(sprintf("  Average window: %d–%d\n", min(avg_years), max(avg_years)))
  cat("══════════════════════════════════════════════════════════════\n")

  # --- Trade weights: OECD → IMF DOTS → World Bank GDP → hardcoded ---
  W_trade <- .oecd_trade(countries, avg_years)

  if (is.null(W_trade)) {
    W_trade <- .imf_dots(countries, avg_years)
  }

  if (is.null(W_trade)) {
    W_trade <- .wb_gdp_weights(countries, avg_years)
  }

  if (is.null(W_trade)) {
    W_trade <- .hardcoded_gdp_weights(countries)
  }

  # --- Optional financial-linkage blend ---
  if (finance_alpha > 0) {
    W_fin <- .bis_finance(countries, avg_years)
    if (!is.null(W_fin)) {
      W_trade <- (1 - finance_alpha) * W_trade + finance_alpha * W_fin
      W_trade <- .row_norm(W_trade)
      cat(sprintf("  Blended: %.0f%% trade + %.0f%% BIS financial.\n",
                  (1 - finance_alpha) * 100, finance_alpha * 100))
    } else {
      message("  BIS data unavailable – using trade weights only.")
    }
  }

  cat("  Weight matrix ready.\n")
  W_trade
}
