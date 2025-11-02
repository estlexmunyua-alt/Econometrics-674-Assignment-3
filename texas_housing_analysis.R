# =============================================================================
# TEXAS HOUSING MARKET ANALYSIS - VECTOR AUTOREGRESSIVE MODEL
# Author: Alex Munyua 
# Purpose: Econometrics 674 Assignment 3 - VAR Forecasting
# Date: October 2024
# =============================================================================

# ----------------------------
# REQUIRED LIBRARIES
# ----------------------------
libs <- c("vars", "zoo", "ggplot2", "dplyr", "tidyr", "readxl", "writexl")
invisible(lapply(libs, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}))

# ----------------------------
# SAFE P-VALUE EXTRACTOR
# ----------------------------
safe_p_value <- function(obj) {
  # Recursively search for a named element "p.value" or "p.value" within htest-like objects
  finder <- function(x) {
    if (is.null(x)) return(NA_real_)
    if (is.list(x)) {
      if (!is.null(x$p.value) && is.numeric(x$p.value)) return(as.numeric(x$p.value))
      # search children
      for (el in x) {
        pv <- finder(el)
        if (!is.na(pv)) return(pv)
      }
      return(NA_real_)
    }
    if (inherits(x, "htest") && !is.null(x$p.value)) return(as.numeric(x$p.value))
    # fallback: look for numeric entries named like "p.value" in attributes
    nm <- names(x)
    if (!is.null(nm) && any(grepl("p.value|p_value|pval", nm, ignore.case = TRUE))) {
      return(as.numeric(x[grepl("p.value|p_value|pval", nm, ignore.case = TRUE)][[1]]))
    }
    return(NA_real_)
  }
  tryCatch(finder(obj), error = function(e) NA_real_)
}

# ----------------------------
# RUN DIAGNOSTICS (robust)
# ----------------------------
message("Running model diagnostics...")
diagnostic_results <- list()

# Serial correlation (adaptive lags)
lags_choice <- min(16, floor(effective_obs / 3))
diagnostic_results$serial_correlation <- tryCatch(
  serial.test(var_model, lags.pt = lags_choice, type = "PT.asymptotic"),
  error = function(e) { warning("serial.test failed: ", e$message); NULL }
)

diagnostic_results$arch_effects <- tryCatch(
  arch.test(var_model, lags.multi = 5),
  error = function(e) { warning("arch.test failed: ", e$message); NULL }
)

diagnostic_results$normality <- tryCatch(
  normality.test(var_model),
  error = function(e) { warning("normality.test failed: ", e$message); NULL }
)

serial_p <- safe_p_value(diagnostic_results$serial_correlation)
arch_p <- safe_p_value(diagnostic_results$arch_effects)

message("  Serial correlation p-value: ", ifelse(is.na(serial_p), "NA", formatC(serial_p, digits = 4)))
message("  ARCH effects p-value: ", ifelse(is.na(arch_p), "NA", formatC(arch_p, digits = 4)))

# ----------------------------
# STRUCTURAL BREAKS (defensive)
# ----------------------------
message("Analyzing structural breaks in the time series...")
break_analysis <- list()

for (series_name in c("Sales", "Permits")) {
  series_data <- tryCatch({
    # ensure numeric vector / ts accepted by breakpoints()
    as.numeric(var_data[, series_name])
  }, error = function(e) {
    warning("Could not extract series ", series_name, ": ", e$message)
    NULL
  })
  if (is.null(series_data)) next

  breaks <- tryCatch(breakpoints(series_data ~ 1), error = function(e) {
    warning("Break detection failed for ", series_name, ": ", e$message)
    NULL
  })

  if (!is.null(breaks) && !all(is.na(breaks$breakpoints))) {
    # convert break indexes to dates if var_data has rownames / dates
    b_idx <- breaks$breakpoints[!is.na(breaks$breakpoints)]
    if (!is.null(time(var_data))) {
      # time() for ts; otherwise attempt zoo index
      idx_time <- tryCatch(time(var_data)[b_idx], error = function(e) NA)
      break_analysis[[series_name]] <- idx_time
      message("  ", series_name, " breaks at indexes: ", paste(b_idx, collapse = ", "))
    } else {
      break_analysis[[series_name]] <- b_idx
      message("  ", series_name, " breaks at indexes: ", paste(b_idx, collapse = ", "))
    }
  }
}

# Example system-wide detection function might be custom; keep existing call but safe
system_breaks <- tryCatch(detect_system_breaks(var_model), error = function(e) {
  warning("detect_system_breaks failed: ", e$message); NULL
})
if (!is.null(system_breaks)) {
  break_analysis$system_wide <- system_breaks
  message("  System-wide breaks at: ", paste(system_breaks, collapse = ", "))
}

# ----------------------------
# FORECAST GENERATION & BACK-TRANSFORMATION
# ----------------------------
message("Generating forecasts...")

# 1. predict from VAR
var_forecasts <- tryCatch(predict(var_model, n.ahead = analysis_config$forecast_months), error = function(e) {
  stop("VAR predict failed: ", e$message)
})

# 2. robust converter from log-diff forecasts to levels
#    fcst_vec: numeric vector of forecasted log-differences (or growth rates); last_actual_value: numeric scalar
convert_log_diff_forecast <- function(fcst_vec, last_actual_value, type = c("logdiff", "pct")) {
  type <- match.arg(type)
  if (length(fcst_vec) == 0 || is.null(last_actual_value) || !is.finite(last_actual_value)) return(rep(NA_real_, length(fcst_vec)))
  # ensure numeric
  fcst <- as.numeric(fcst_vec)
  if (type == "logdiff") {
    # if fcst are forecasts of log(d_t) - log(d_{t-1}), cumulative levels:
    cumulative_growth <- exp(cumsum(fcst))
  } else {
    # pct: fcst are decimal percentage changes, e.g. 0.01 for 1%
    cumulative_growth <- cumprod(1 + fcst)
  }
  as.numeric(last_actual_value * cumulative_growth)
}

# 3. extract fcst matrices carefully (vars::predict returns list with $fcst entries)
if (!is.null(var_forecasts$fcst)) {
  # Check names
  fcst_sales <- tryCatch(var_forecasts$fcst$Sales[, "fcst"], error = function(e) {
    stop("Cannot extract Sales forecasts: ", e$message)
  })
  fcst_sales_lower <- tryCatch(var_forecasts$fcst$Sales[, "lower"], error = function(e) NA_real_)
  fcst_sales_upper <- tryCatch(var_forecasts$fcst$Sales[, "upper"], error = function(e) NA_real_)

  fcst_permits <- tryCatch(var_forecasts$fcst$Permits[, "fcst"], error = function(e) {
    stop("Cannot extract Permits forecasts: ", e$message)
  })
  fcst_permits_lower <- tryCatch(var_forecasts$fcst$Permits[, "lower"], error = function(e) NA_real_)
  fcst_permits_upper <- tryCatch(var_forecasts$fcst$Permits[, "upper"], error = function(e) NA_real_)

  # 4. last actual values (ensure numeric scalar)
  last_sales <- as.numeric(tail(TexasSales, 1))
  last_permits <- as.numeric(tail(TexasPermits, 1))

  # 5. back-transform (assume forecasts are log-differences; change `type` if not)
  sales_forecast_levels <- convert_log_diff_forecast(fcst_sales, last_sales, type = "logdiff")
  permits_forecast_levels <- convert_log_diff_forecast(fcst_permits, last_permits, type = "logdiff")

  sales_lower <- convert_log_diff_forecast(fcst_sales_lower, last_sales, type = "logdiff")
  sales_upper <- convert_log_diff_forecast(fcst_sales_upper, last_sales, type = "logdiff")
  permits_lower <- convert_log_diff_forecast(fcst_permits_lower, last_permits, type = "logdiff")
  permits_upper <- convert_log_diff_forecast(fcst_permits_upper, last_permits, type = "logdiff")

  # 6. build forecast_dates (requires zoo::as.yearmon on time(TexasSales))
  last_historical_date <- time(TexasSales)[length(TexasSales)]
  forecast_dates <- as.Date(as.yearmon(last_historical_date) + (1:analysis_config$forecast_months)/12)

  # 7. assemble forecast_results dataframe
  forecast_results <- data.frame(
    Date = forecast_dates,
    Sales_Forecast = sales_forecast_levels,
    Sales_Lower = sales_lower,
    Sales_Upper = sales_upper,
    Permits_Forecast = permits_forecast_levels,
    Permits_Lower = permits_lower,
    Permits_Upper = permits_upper
  )
} else {
  stop("var_forecasts$fcst is NULL; cannot build forecasts.")
}

message("Forecasts generated for period: ", format(forecast_results$Date[1]), " to ", format(tail(forecast_results$Date, 1)))

# ----------------------------
# PLOTTING (safe)
# ----------------------------
message("Creating visualization...")

create_plot_data <- function(ts_object, series_name) {
  data.frame(
    Date = as.Date(as.yearmon(time(ts_object))),
    Value = as.numeric(ts_object),
    Series = series_name
  )
}

historical_data <- bind_rows(
  create_plot_data(TexasSales, "Historical Sales"),
  create_plot_data(TexasPermits, "Historical Permits")
)

# ensure pivot_longer available (we loaded tidyr above)
forecast_data <- forecast_results %>%
  select(Date, Sales_Forecast, Permits_Forecast) %>%
  tidyr::pivot_longer(cols = -Date, names_to = "Series", values_to = "Value") %>%
  dplyr::mutate(Series = dplyr::case_when(
    Series == "Sales_Forecast" ~ "Forecast Sales",
    Series == "Permits_Forecast" ~ "Forecast Permits",
    TRUE ~ Series
  ))

forecast_start <- as.Date(as.yearmon(last_historical_date))

final_plot <- ggplot() +
  geom_line(data = historical_data, aes(x = Date, y = Value, color = Series), linewidth = 0.8) +
  geom_line(data = forecast_data, aes(x = Date, y = Value, color = Series), linewidth = 1.2, linetype = "dashed") +
  geom_ribbon(data = forecast_results, aes(x = Date, ymin = Sales_Lower, ymax = Sales_Upper), inherit.aes = FALSE, alpha = 0.2) +
  geom_ribbon(data = forecast_results, aes(x = Date, ymin = Permits_Lower, ymax = Permits_Upper), inherit.aes = FALSE, alpha = 0.2) +
  geom_vline(xintercept = forecast_start, linetype = "dotdash", color = "black", alpha = 0.6, linewidth = 0.8) +
  annotate("text", x = forecast_start, y = Inf, label = "Forecast Start", vjust = 2, hjust = 1.1, size = 3, fontface = "italic") +
  labs(
    title = paste0("Texas Housing Market Analysis: VAR(", optimal_lags, ") Forecast"),
    subtitle = paste0("Serial p=", ifelse(is.na(serial_p), "NA", round(serial_p, 4)),
                      " | ARCH p=", ifelse(is.na(arch_p), "NA", round(arch_p, 4))),
    y = "Units",
    x = "Date",
    caption = "Dashed = forecast; ribbons = 95% intervals"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Historical Sales" = "red", "Historical Permits" = "blue", "Forecast Sales" = "darkred", "Forecast Permits" = "darkblue"))

print(final_plot)

# ----------------------------
# RESULTS EXPORT (safe)
# ----------------------------
message("Saving analysis results...")

analysis_results <- list(
  metadata = list(
    analysis_date = Sys.Date(),
    data_period = paste0(as.Date(as.yearmon(time(TexasSales))[1]), " to ", as.Date(as.yearmon(tail(time(TexasSales), 1)))),
    total_observations = length(TexasSales),
    var_lag_order = optimal_lags
  ),
  var_model = var_model,
  forecasts = forecast_results,
  diagnostics = list(serial_correlation_p = serial_p, arch_effects_p = arch_p, full = diagnostic_results),
  structural_breaks = break_analysis,
  configuration = analysis_config
)

saveRDS(analysis_results, "texas_housing_analysis_results.rds")
write.csv(forecast_results, "texas_housing_forecasts.csv", row.names = FALSE)

# reproducibility report
repro_info <- c(
  "ANALYSIS REPRODUCIBILITY INFORMATION",
  "====================================",
  paste("Analysis performed:", as.character(Sys.time())),
  paste("Data file:", analysis_config$data_file),
  paste("Analysis period:", analysis_results$metadata$data_period),
  paste("VAR lag order:", optimal_lags),
  "",
  capture.output(sessionInfo())
)
writeLines(repro_info, "reproducibility_report.txt")

# nicer separator printing
cat(paste0("\n", strrep("=", 60), "\n"))
cat("TEXAS HOUSING MARKET ANALYSIS COMPLETE\n")
cat(paste0(strrep("=", 60), "\n"))
