# Time Series Analysis Utilities
# Helper functions for time series network analysis

#' Check Time Series Data
#'
#' Validates a time series vector and provides a summary of its characteristics.
#' This function is useful for initial data inspection before further analysis.
#'
#' @param ts_data A numeric vector representing the time series.
#' @param min_length Integer. The minimum required length of the time series.
#'   If the time series is shorter than this, an error will be thrown. Defaults to 10.
#'
#' @return A list containing:
#'   \item{valid}{Logical. `TRUE` if the input `ts_data` passes basic validation checks.}
#'   \item{summary}{A list of summary statistics for the time series, including `length`,
#'     `mean`, `sd` (standard deviation), `min`, `max`, `na_count` (number of missing values),
#'     and `range` (max - min).}
#'   \item{has_missing}{Logical. `TRUE` if the time series contains any `NA` values, `FALSE` otherwise.}
#' @export
ts_check <- function(ts_data, min_length = 10) {
  # Basic validation
  if (is.null(ts_data)) {
    stop("ts_data cannot be NULL")
  }
  if (!is.numeric(ts_data)) {
    stop("ts_data must be numeric")
  }
  if (length(ts_data) < min_length) {
    stop("ts_data must have at least ", min_length, " observations")
  }
  # Summary statistics
  summary_stats <- list(
    length = length(ts_data),
    mean = mean(ts_data, na.rm = TRUE),
    sd = sd(ts_data, na.rm = TRUE),
    min = min(ts_data, na.rm = TRUE),
    max = max(ts_data, na.rm = TRUE),
    na_count = sum(is.na(ts_data)),
    range = max(ts_data, na.rm = TRUE) - min(ts_data, na.rm = TRUE)
  )
  return(list(
    valid = TRUE,
    summary = summary_stats,
    has_missing = summary_stats$na_count > 0
  ))
}

#' Calculate Optimal Window Size
#'
#' Suggests an optimal window size for rolling window analysis based on the
#' length of the time series and a desired number of resulting windows.
#'
#' @param ts_data A numeric vector representing the time series.
#' @param max_windows Integer. The maximum desired number of windows to be created.
#'   Defaults to 20.
#' @param min_window_size Integer. The minimum allowed size for a window.
#'   Defaults to 5.
#'
#' @return A list containing:
#'   \item{window_size}{The suggested optimal window size.}
#'   \item{resulting_windows}{The actual number of windows that would be created
#'     with the suggested `window_size` and a step size equal to `window_size`.}
#'   \item{data_length}{The length of the input time series.}
#'   \item{recommendation}{A character string providing a human-readable recommendation.}
#' @export
ts_window_size <- function(ts_data, max_windows = 20, min_window_size = 5) {
  n <- length(ts_data)
  # Calculate window size that gives approximately max_windows
  suggested_size <- floor(n / max_windows)
  # Ensure minimum size
  suggested_size <- max(suggested_size, min_window_size)
  # Ensure not too large
  suggested_size <- min(suggested_size, floor(n / 3))
  # Calculate resulting number of windows
  resulting_windows <- floor((n - suggested_size) / suggested_size) + 1
  return(list(
    window_size = suggested_size,
    resulting_windows = resulting_windows,
    data_length = n,
    recommendation = paste("Use window_size =", suggested_size,
                           "for approximately", resulting_windows, "windows")
  ))
}

#' Distance Matrix Summary
#'
#' Analyzes the properties of a distance matrix and provides summary statistics
#' and recommendations for network construction methods based on the distance scale.
#'
#' @param distance_matrix A square numeric distance matrix.
#'
#' @return A list containing:
#'   \item{statistics}{A list of basic statistics of the distances, including `matrix_size`,
#'     `min_distance`, `max_distance`, `mean_distance`, `median_distance`,
#'     `sd_distance`, `q25` (25th percentile), `q75` (75th percentile),
#'     `zero_distances` (count of zero distances), and `range` (max - min).}
#'   \item{scale_category}{A character string categorizing the scale of distances
#'     ("very_large", "large", "medium", "small").}
#'   \item{recommended_methods}{A character vector of network construction methods
#'     recommended for the detected distance scale.}
#'   \item{summary_text}{A human-readable summary of the distance matrix properties.}
#' @importFrom stats median quantile sd
#' @importFrom utils globalVariables
#' @name ts-utils-package

utils::globalVariables(c("ts_data", "min_length", "max_windows", "min_window_size", "distance_matrix"))
ts_distance_summary <- function(distance_matrix) {
  if (!is.matrix(distance_matrix) || nrow(distance_matrix) != ncol(distance_matrix)) {
    stop("Input must be a square matrix")
  }
  # Extract upper triangular distances (excluding diagonal)
  distances <- distance_matrix[upper.tri(distance_matrix)]
  # Calculate statistics
  stats <- list(
    matrix_size = nrow(distance_matrix),
    min_distance = min(distances),
    max_distance = max(distances),
    mean_distance = mean(distances),
    median_distance = median(distances),
    sd = sd(distances),
    q25 = quantile(distances, 0.25),
    q75 = quantile(distances, 0.75),
    zero_distances = sum(distances == 0),
    range = max(distances) - min(distances)
  )
  # Determine scale and recommend methods
  if (stats$max_distance > 1000) {
    scale <- "very_large"
    recommended_methods <- c("knn", "percentile")
  } else if (stats$max_distance > 100) {
    scale <- "large"
    recommended_methods <- c("knn", "percentile", "threshold")
  } else if (stats$max_distance > 10) {
    scale <- "medium"
    recommended_methods <- c("threshold", "percentile")
  } else {
    scale <- "small"
    recommended_methods <- c("epsilon", "threshold")
  }
  return(list(
    statistics = stats,
    scale_category = scale,
    recommended_methods = recommended_methods,
    summary_text = paste("Distance range:", round(stats$min_distance, 2), "to",
                         round(stats$max_distance, 2), "- Scale:", scale)
  ))
}

#' Package Information
#'
#' Displays an overview of the main functions and capabilities provided by the
#' `tsnetworks` package. This serves as a quick reference for users.
#'
#' @return A list of categorized function names, printed to the console.
#' @export
ts_info <- function() {
  functions_list <- list(
    distance_functions = c("ts_distance", "ts_distance_pair", "ts_distance_methods"),
    network_functions = c("ts_network", "ts_network_auto", "ts_network_methods", "ts_network_recommend"),
    plotting_functions = c("plot_ts_network", "plot_timeseries"),
    utility_functions = c("ts_check", "ts_window_size", "ts_distance_summary", "ts_info"),
    discretization_functions = c("stna", "compute_quantile_states"),
    regime_detection_functions = c("rolling_measures", "detect_regime"),
    conversion_functions = c("distance_to_similarity", "normalize_matrix", "prepare_for_qgraph", "quick_qgraph", "show_conversion_examples")
  )

  cat("=== Time Series Network Analysis Package ===\n\n")
  cat("This package provides tools for analyzing time series data through network-based approaches.\n")
  cat("Key functionalities include Visibility Graphs, Dynamic Complexity and Regime Detection,\n")
  cat("and Time Series Proximity Networks.\n\n")

  cat("MAIN FUNCTION CATEGORIES:\n")
  for (category_name in names(functions_list)) {
    cat(paste0("  ", toupper(gsub("_", " ", category_name)), ":\n"))
    for (func_name in functions_list[[category_name]]) {
      cat(paste0("    - ", func_name, "\n"))
    }
    cat("\n")
  }

  cat("For detailed help on any function, use `?function_name` (e.g., `?ts_distance`).\n")
  cat("For a full workflow example, refer to the package vignettes (e.g., `vignette(\"proximity-networks\")`).\n")

  invisible(functions_list)
}

#' Generate sliding windows from a time series
#'
#' This internal helper function creates a list of sliding windows from a given time series.
#'
#' @param series Numeric vector. The input time series.
#' @param window_size Integer. The size of each window.
#' @param step_size Integer. The step size between the start of consecutive windows.
#' @param drop_last_if_shorter Logical. If TRUE, drops the last window if it's shorter than window_size.
#' @return A list containing:
#'   \item{windows}{A list of numeric vectors, each representing a window.}
#'   \item{start_indices}{A vector of starting indices for each window.}
#'   \item{actual_lengths}{A vector of actual lengths for each window.}
#' @keywords internal
.generate_windows <- function(series, window_size, step_size, drop_last_if_shorter = FALSE) {
  n <- length(series)
  if (n < window_size) {
    warning("Series length is less than window size. No windows generated.")
    return(list(windows = list(), start_indices = integer(0), actual_lengths = integer(0)))
  }
  
  starts <- seq(1, n - window_size + 1, by = step_size)
  windows <- list()
  actual_lengths <- integer(length(starts))
  
  for (i in seq_along(starts)) {
    start_idx <- starts[i]
    end_idx <- min(start_idx + window_size - 1, n)
    
    current_window <- series[start_idx:end_idx]
    current_length <- length(current_window)
    
    if (drop_last_if_shorter && current_length < window_size) {
      # Skip this window if it's shorter and we're dropping partials
      next
    }
    
    windows[[length(windows) + 1]] <- current_window
    actual_lengths[length(windows)] <- current_length
  }
  
  # Filter out any NAs that might have resulted from dropping windows
  starts_filtered <- starts[1:length(windows)]
  actual_lengths_filtered <- actual_lengths[1:length(windows)]
  
  return(list(windows = windows, start_indices = starts_filtered, actual_lengths = actual_lengths_filtered))
}
