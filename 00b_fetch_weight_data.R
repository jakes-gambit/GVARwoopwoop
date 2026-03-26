###############################################################################
#  00b_fetch_weight_data.R  –  Fetch Real Bilateral Trade & Financial Flow
#                               Data for GVAR Weight Matrix Construction
#
#  Sources supported:
#    1. OECD Bilateral Trade in Goods  (OECD.Stat via 'OECD' R package)
#    2. BIS Consolidated Banking Statistics  (cross-border claims, 'BIS' package)
#    3. IMF Direction of Trade Statistics  (DOTS, fallback via IMF Data API)
#    4. FRED bilateral import proxies  (for a quick FREd-only fallback)
#
#  Output: a named N×N weight matrix ready for build_weight_matrix() in
#          02_data_preparation.R.
#
#  PREREQUISITES (install once):
#    install.packages(c("OECD", "BIS", "dplyr", "tidyr", "fredr"))
#
#  USAGE:
#    source("00b_fetch_weight_data.R")
#    W_trade     <- build_oecd_trade_weights(countries, avg_years = 2014:2016)
#    W_finance   <- build_bis_financial_weights(countries, avg_years = 2014:2016)
#    W_combined  <- combine_weight_matrices(W_trade, W_finance, alpha = 0.5)
###############################################################################

library(dplyr)
library(tidyr)

# ─────────────────────────────────────────────────────────────────────────────
# Helper: row-normalise a matrix (zero diagonal, rows sum to 1)
# ─────────────────────────────────────────────────────────────────────────────
row_normalise <- function(W) {
  diag(W) <- 0
  rs <- rowSums(W, na.rm = TRUE)
  rs[rs == 0] <- 1   # avoid division by zero (isolated country → equal weights)
  W / rs
}


# ═════════════════════════════════════════════════════════════════════════════
# 1.  OECD Bilateral Trade in Goods
#     Dataset: "OECD ITCS" – SITC or HS based bilateral trade flows
#     Package: OECD  (install.packages("OECD"))
# ═════════════════════════════════════════════════════════════════════════════

#' Build a trade-flow weight matrix from OECD bilateral trade statistics.
#'
#' Uses the OECD "TRADE_IN_GOODS" dataset (Total Trade, USD, annual).
#' Weights are constructed as trade_ij / (sum_j trade_ij) — the share of
#' country j in country i's total bilateral trade.
#'
#' @param countries  Character vector of ISO-3 country codes matching your
#'                   GVAR unit names  (e.g. c("DEU","FRA","ITA","ESP"))
#' @param avg_years  Integer vector of years to average (default 3-year avg)
#' @return           Named N×N row-normalised weight matrix
build_oecd_trade_weights <- function(countries,
                                     avg_years = 2014:2016) {

  if (!requireNamespace("OECD", quietly = TRUE))
    stop("Package 'OECD' required. Run: install.packages('OECD')")

  cat("[OECD Trade] Fetching bilateral trade data ...\n")

  # OECD ITCS dataset ID
  dataset_id <- "TRADE_IN_GOODS"

  # Filter for selected countries as reporters and partners
  filter_list <- list(
    REPORTER = countries,
    PARTNER  = countries,
    MEASURE  = "USD",          # USD millions
    TIME     = as.character(avg_years)
  )

  raw <- tryCatch(
    OECD::get_dataset(dataset_id, filter = filter_list, start_time = min(avg_years),
                      end_time = max(avg_years)),
    error = function(e) {
      message("[OECD Trade] Download failed: ", e$message)
      return(NULL)
    }
  )

  if (is.null(raw) || nrow(raw) == 0) {
    message("[OECD Trade] No data returned. Returning NULL.")
    return(NULL)
  }

  # Sum exports + imports as measure of bilateral linkage, average over years
  bilateral <- raw %>%
    dplyr::filter(REPORTER %in% countries, PARTNER %in% countries,
                  REPORTER != PARTNER) %>%
    dplyr::group_by(REPORTER, PARTNER) %>%
    dplyr::summarise(trade = mean(as.numeric(obsValue), na.rm = TRUE),
                     .groups = "drop")

  # Pivot to N×N matrix
  W_mat <- bilateral %>%
    tidyr::pivot_wider(names_from = PARTNER, values_from = trade,
                       values_fill = 0) %>%
    tibble::column_to_rownames("REPORTER") %>%
    as.matrix()

  # Ensure all requested countries appear (add zeros for missing pairs)
  W_mat <- pad_weight_matrix(W_mat, countries)

  cat(sprintf("[OECD Trade] Weight matrix built (%d x %d), years %s.\n",
              nrow(W_mat), ncol(W_mat), paste(range(avg_years), collapse = "–")))

  return(row_normalise(W_mat))
}


# ═════════════════════════════════════════════════════════════════════════════
# 2.  BIS Consolidated Banking Statistics (Cross-border Financial Flows)
#     Dataset: CBS "Immediate counterparty" basis (Table B1)
#     Package: BIS  (install.packages("BIS"))
# ═════════════════════════════════════════════════════════════════════════════

#' Build a financial-flow weight matrix from BIS consolidated banking stats.
#'
#' Uses BIS CBS total claims (USD millions) as a measure of cross-border
#' financial linkages between countries.  The resulting matrix captures
#' banking-sector exposures, which complement trade-based weights.
#'
#' @param countries  Character vector of ISO-3 country codes
#' @param avg_years  Integer vector of years to average (default 3-year avg)
#' @return           Named N×N row-normalised weight matrix
build_bis_financial_weights <- function(countries,
                                        avg_years = 2014:2016) {

  if (!requireNamespace("BIS", quietly = TRUE))
    stop("Package 'BIS' required. Run: install.packages('BIS')")

  cat("[BIS] Fetching consolidated banking statistics ...\n")

  # BIS provides a pre-built tidy dataset; download and filter
  bis_raw <- tryCatch(
    BIS::get_bis("CBS"),   # Consolidated Banking Statistics
    error = function(e) {
      message("[BIS] Download failed: ", e$message)
      return(NULL)
    }
  )

  if (is.null(bis_raw) || nrow(bis_raw) == 0) {
    message("[BIS] No data returned. Returning NULL.")
    return(NULL)
  }

  # Filter: immediate counterparty basis, total claims (all sectors, all maturities)
  bilateral <- bis_raw %>%
    dplyr::filter(
      measure       == "B",         # Immediate counterparty basis
      reporting_country %in% countries,
      counterparty_country %in% countries,
      reporting_country != counterparty_country,
      currency  == "TO1",           # All currencies
      date >= as.Date(paste0(min(avg_years), "-01-01")),
      date <= as.Date(paste0(max(avg_years), "-12-31"))
    ) %>%
    dplyr::group_by(reporting_country, counterparty_country) %>%
    dplyr::summarise(claims = mean(as.numeric(obs_value), na.rm = TRUE),
                     .groups = "drop")

  # Pivot to N×N matrix
  W_mat <- bilateral %>%
    tidyr::pivot_wider(names_from = counterparty_country, values_from = claims,
                       values_fill = 0) %>%
    tibble::column_to_rownames("reporting_country") %>%
    as.matrix()

  W_mat <- pad_weight_matrix(W_mat, countries)

  cat(sprintf("[BIS] Weight matrix built (%d x %d), years %s.\n",
              nrow(W_mat), ncol(W_mat), paste(range(avg_years), collapse = "–")))

  return(row_normalise(W_mat))
}


# ═════════════════════════════════════════════════════════════════════════════
# 3.  IMF Direction of Trade Statistics (DOTS) via IMF Data API
#     No R package needed – plain HTTP JSON request via httr / jsonlite
# ═════════════════════════════════════════════════════════════════════════════

#' Build a trade weight matrix using the IMF DOTS API.
#'
#' Falls back to the IMF REST API when the OECD package is unavailable.
#' Fetches annual exports (TXG_FOB_USD) for all reporter–partner pairs.
#'
#' @param countries  Character vector of ISO-2 country codes (IMF format)
#'                   Note: IMF DOTS uses 2-digit ISO codes (e.g. "DE","FR")
#' @param avg_years  Integer vector of years to average
#' @return           Named N×N row-normalised weight matrix
build_imf_dots_weights <- function(countries_iso2,
                                   avg_years = 2014:2016) {

  if (!requireNamespace("httr", quietly = TRUE) ||
      !requireNamespace("jsonlite", quietly = TRUE))
    stop("Packages 'httr' and 'jsonlite' required. Run: install.packages(c('httr','jsonlite'))")

  cat("[IMF DOTS] Fetching bilateral exports via IMF Data API ...\n")

  N     <- length(countries_iso2)
  W_mat <- matrix(0, N, N,
                  dimnames = list(countries_iso2, countries_iso2))

  base_url <- "https://www.imf.org/external/datamapper/api/v1/TXG_FOB_USD"

  for (reporter in countries_iso2) {
    url <- sprintf("%s/%s", base_url, reporter)

    resp <- tryCatch(
      httr::GET(url, httr::timeout(30)),
      error = function(e) NULL
    )

    if (is.null(resp) || httr::status_code(resp) != 200) {
      message(sprintf("[IMF DOTS] Skipping %s: API error", reporter))
      next
    }

    parsed <- jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"),
                                  simplifyVector = FALSE)

    partners_data <- parsed$values$TXG_FOB_USD[[reporter]]

    if (is.null(partners_data)) next

    for (partner in countries_iso2) {
      if (reporter == partner) next
      vals <- partners_data[[partner]]
      if (is.null(vals)) next

      year_vals <- unlist(vals[as.character(avg_years)])
      if (length(year_vals) > 0) {
        W_mat[reporter, partner] <- mean(as.numeric(year_vals), na.rm = TRUE)
      }
    }

    Sys.sleep(0.3)   # polite rate limiting
  }

  cat(sprintf("[IMF DOTS] Weight matrix built (%d x %d).\n", N, N))
  return(row_normalise(W_mat))
}


# ═════════════════════════════════════════════════════════════════════════════
# 4.  FRED-Based Fallback – Approximate Bilateral Weights
#     Uses country-level total imports/exports from FRED as a
#     simple proxy when bilateral sources are unavailable.
#     NOTE: This produces ONLY approximate weights; prefer OECD/BIS/IMF above.
# ═════════════════════════════════════════════════════════════════════════════

#' Build approximate trade weights from FRED country-level trade series.
#'
#' When no bilateral data source is available, this helper constructs a
#' symmetric weight matrix using each country's total trade (imports + exports)
#' from FRED.  Weights are proportional to trade size — larger economies receive
#' higher weight.  This is a last-resort approximation only.
#'
#' @param countries  Character vector of country codes matching fred_series names
#' @param avg_years  Integer vector of years to average
#' @param api_key    FRED API key (uses fredr_set_key() if already set)
#' @return           Named N×N symmetric row-normalised weight matrix
build_fred_approx_weights <- function(countries,
                                      avg_years = 2014:2016,
                                      api_key   = NULL) {

  if (!requireNamespace("fredr", quietly = TRUE))
    stop("Package 'fredr' required.")

  if (!is.null(api_key)) fredr::fredr_set_key(api_key)

  # FRED total trade series (IMP + EXP in billions, SAAR) for major countries
  # Extend with additional series as needed
  fred_trade_series <- c(
    USA = "BOPGSTB",     # US trade balance → use separate IMP/EXP
    DEU = "XTEXVA01DEQ667S",
    FRA = "XTEXVA01FRQ667S",
    ITA = "XTEXVA01ITQ667S",
    GBR = "XTEXVA01GBQ667S",
    JPN = "XTEXVA01JPQ667S",
    CAN = "XTEXVA01CAQ667S",
    AUS = "XTEXVA01AUQ667S",
    ESP = "XTEXVA01ESQ667S"
  )

  # Filter to requested countries
  avail <- intersect(countries, names(fred_trade_series))

  trade_size <- setNames(rep(NA_real_, length(countries)), countries)

  for (ctry in avail) {
    df <- tryCatch(
      fredr::fredr(series_id   = fred_trade_series[[ctry]],
                   observation_start = as.Date(paste0(min(avg_years), "-01-01")),
                   observation_end   = as.Date(paste0(max(avg_years), "-12-31")),
                   frequency   = "a"),
      error = function(e) NULL
    )

    if (!is.null(df) && nrow(df) > 0) {
      trade_size[[ctry]] <- mean(abs(df$value), na.rm = TRUE)
    }
    Sys.sleep(0.15)
  }

  # For missing countries: impute with the median of available values
  med_val <- median(trade_size, na.rm = TRUE)
  if (is.na(med_val)) med_val <- 1
  trade_size[is.na(trade_size)] <- med_val

  # Build outer-product symmetric approximation
  N <- length(countries)
  W_mat <- outer(trade_size, trade_size) / sum(trade_size)
  dimnames(W_mat) <- list(countries, countries)

  cat(sprintf("[FRED Approx] Approximate weight matrix built (%d x %d).\n", N, N))
  message("[FRED Approx] WARNING: These are approximate weights based on trade size,",
          " NOT bilateral flows. Use OECD/BIS/IMF sources when possible.")

  return(row_normalise(W_mat))
}


# ═════════════════════════════════════════════════════════════════════════════
# 5.  Combine Multiple Weight Matrices (Trade + Financial)
# ═════════════════════════════════════════════════════════════════════════════

#' Combine trade and financial weight matrices with a user-specified blend.
#'
#' Following Dees et al. (2007), weights can blend trade and financial linkages:
#'   W_combined = alpha * W_trade + (1 - alpha) * W_finance
#'
#' @param W_trade    N×N trade weight matrix (from build_oecd_trade_weights)
#' @param W_finance  N×N financial weight matrix (from build_bis_financial_weights)
#' @param alpha      Blending weight for trade (0 = finance only, 1 = trade only)
#' @return           N×N row-normalised combined weight matrix
combine_weight_matrices <- function(W_trade, W_finance, alpha = 0.5) {

  stopifnot(alpha >= 0, alpha <= 1)

  # Align country order
  countries <- rownames(W_trade)
  if (!all(countries == rownames(W_finance)))
    W_finance <- W_finance[countries, countries]

  W_combined <- alpha * W_trade + (1 - alpha) * W_finance

  cat(sprintf("[Weights] Combined matrix: %.0f%% trade + %.0f%% financial.\n",
              alpha * 100, (1 - alpha) * 100))

  return(row_normalise(W_combined))
}


# ═════════════════════════════════════════════════════════════════════════════
# 6.  Utility: Pad a weight matrix to include all requested countries
# ═════════════════════════════════════════════════════════════════════════════

#' Ensure a weight matrix contains all requested countries.
#'
#' Countries missing from the data receive zero bilateral weights; after
#' row-normalisation they will have equal weight assigned to all partners.
#'
#' @param W_mat     Existing N×N matrix (possibly smaller than requested)
#' @param countries Full set of requested country codes
#' @return          Square matrix with dimnames = countries
pad_weight_matrix <- function(W_mat, countries) {

  N       <- length(countries)
  W_full  <- matrix(0, N, N, dimnames = list(countries, countries))

  present <- intersect(rownames(W_mat), countries)

  if (length(present) > 0) {
    W_full[present, intersect(colnames(W_mat), countries)] <-
      W_mat[present, intersect(colnames(W_mat), countries), drop = FALSE]
  }

  missing_ctry <- setdiff(countries, rownames(W_mat))
  if (length(missing_ctry) > 0) {
    message("[Weights] Countries not found in raw data (will get equal weights): ",
            paste(missing_ctry, collapse = ", "))
  }

  return(W_full)
}


# ═════════════════════════════════════════════════════════════════════════════
# 7.  Master Function: Build Weight Matrix from Best Available Source
# ═════════════════════════════════════════════════════════════════════════════

#' Build a GVAR weight matrix from the best available bilateral data source.
#'
#' Tries sources in order: OECD → IMF DOTS → FRED approx.
#' Optionally blends with BIS financial-flow weights.
#'
#' @param countries       Character vector of ISO-3 country codes (GVAR units)
#' @param avg_years       Integer vector of years to average for weights
#' @param include_finance Logical; if TRUE, blend trade with BIS financial weights
#' @param alpha           Trade weight blend (1 = trade only, 0 = finance only)
#' @param countries_iso2  ISO-2 codes for IMF DOTS (same order as countries)
#' @return                Named N×N row-normalised weight matrix
build_gvar_weights <- function(countries,
                                avg_years       = 2014:2016,
                                include_finance = FALSE,
                                alpha           = 0.5,
                                countries_iso2  = NULL) {

  cat("═══════════════════════════════════════════════════════════════════════\n")
  cat("  Building GVAR weight matrix from bilateral data\n")
  cat("═══════════════════════════════════════════════════════════════════════\n\n")

  W_trade <- NULL

  # --- Try OECD first ---
  if (requireNamespace("OECD", quietly = TRUE)) {
    W_trade <- tryCatch(
      build_oecd_trade_weights(countries, avg_years),
      error = function(e) {
        message("[Weights] OECD failed: ", e$message)
        NULL
      }
    )
  }

  # --- Fallback: IMF DOTS ---
  if (is.null(W_trade) && !is.null(countries_iso2)) {
    W_trade <- tryCatch(
      build_imf_dots_weights(countries_iso2, avg_years),
      error = function(e) {
        message("[Weights] IMF DOTS failed: ", e$message)
        NULL
      }
    )
  }

  # --- Last resort: FRED approximation ---
  if (is.null(W_trade)) {
    message("[Weights] Falling back to FRED approximate weights.")
    W_trade <- build_fred_approx_weights(countries, avg_years)
  }

  # --- Optionally blend with BIS financial weights ---
  if (include_finance && requireNamespace("BIS", quietly = TRUE)) {
    W_finance <- tryCatch(
      build_bis_financial_weights(countries, avg_years),
      error = function(e) {
        message("[Weights] BIS failed, using trade weights only: ", e$message)
        NULL
      }
    )

    if (!is.null(W_finance)) {
      W_trade <- combine_weight_matrices(W_trade, W_finance, alpha = alpha)
    }
  }

  cat("\n[Weights] Done. Pass this matrix to prepare_gvar_dataset() as 'weights'.\n")
  return(W_trade)
}


# ═════════════════════════════════════════════════════════════════════════════
# 8.  Example Usage (commented out – run manually)
# ═════════════════════════════════════════════════════════════════════════════

# selected_countries <- c("DEU", "FRA", "ITA", "ESP")
#
# # Option A: Trade weights only (OECD → IMF fallback → FRED approx)
# W_raw <- build_gvar_weights(
#   countries = selected_countries,
#   avg_years = 2014:2016
# )
#
# # Option B: Trade + financial blend (50/50)
# W_raw <- build_gvar_weights(
#   countries       = selected_countries,
#   avg_years       = 2014:2016,
#   include_finance = TRUE,
#   alpha           = 0.5
# )
#
# # Then pass to GVAR preparation (02_data_preparation.R):
# gvar_data_obj <- prepare_gvar_dataset(
#   data_list   = sim_data,
#   weights     = W_raw,
#   p_lag       = p_vec,
#   q_lag       = q_vec,
#   global_vars = GLOBAL_VARS
# )

cat("Note: Content generated using AI - expert verification recommended.\n")
