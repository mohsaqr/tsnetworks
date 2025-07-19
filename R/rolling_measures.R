#' Calculate Rolling Window Measures for Time Series Data
#'
#' Computes various rolling window measures including dynamic complexity,
#' fluctuation, distribution, autocorrelation, and basic statistics for
#' univariate or multivariate time series data. Supports both grouped
#' (local) and ungrouped (global) calculations.
#'
#' @param data A data frame containing time series data.
#' @param id_col Character string specifying the column name for grouping
#'   identifiers. If \code{NULL}, global calculations are performed.
#'   Default is \code{NULL}.
#' @param ts_cols Character vector specifying column names containing time
#'   series data. If \code{NULL}, all numeric columns (excluding \code{id_col})
#'   are used. Default is \code{NULL}.
#' @param measures Character vector specifying which measures to calculate.
#'   Options: "complexity", "fluctuation", "distribution", "autocorrelation",
#'   "max", "min", "variance". If \code{NULL}, all measures are calculated.
#'   Default is \code{NULL}.
#' @param window_width Positive integer specifying the rolling window size.
#'   Must be at least 2. Default is 7.
#' @param scale_min Numeric value for minimum scale used in complexity-based
#'   measures. If \code{NULL}, computed from data. Default is \code{NULL}.
#' @param scale_max Numeric value for maximum scale used in complexity-based
#'   measures. If \code{NULL}, computed from data. Default is \code{NULL}.
#' @param rescale Logical indicating whether to rescale complexity-based
#'   measures to original data range. Default is \code{FALSE}.
#' @param na_action Character string specifying how to handle missing values.
#'   Options: "omit" (remove), "fail" (stop on NA), "pass" (keep).
#'   Default is "pass".
#'
#' @return A data frame containing the original data with additional columns
#'   for each calculated measure. Column names follow the pattern
#'   \code{<original_column>_<measure>}. For multivariate data, additional
#'   columns \code{mean_<measure>} and \code{median_<measure>} contain
#'   row-wise statistics across variables.
#'
#' @details
#' The function calculates the following measures:
#' \itemize{
#'   \item \strong{complexity}: Product of fluctuation and distribution measures
#'   \item \strong{fluctuation}: Root mean square of successive differences
#'   \item \strong{distribution}: Deviation from uniform distribution
#'   \item \strong{autocorrelation}: Lag-1 autocorrelation coefficient
#'   \item \strong{max}: Rolling maximum
#'   \item \strong{min}: Rolling minimum
#'   \item \strong{variance}: Rolling variance
#' }
#'
#' Scale-dependent measures (complexity, fluctuation, distribution) require
#' \code{scale_min} and \code{scale_max} parameters. If not provided, they
#' are computed from the data range.
#'
#' @examples
#' # Create sample data
#' set.seed(123)
#' n <- 100
#' data <- data.frame(
#'   id = rep(c("A", "B"), each = n/2),
#'   time = rep(1:(n/2), 2),
#'   ts1 = c(sin(seq(0, 4*pi, length.out = n/2)) + rnorm(n/2, 0, 0.1),
#'           cos(seq(0, 4*pi, length.out = n/2)) + rnorm(n/2, 0, 0.1)),
#'   ts2 = rnorm(n)
#' )
#'
#' # Calculate all measures with grouping
#' result <- rolling_measures(data, id_col = "id", window_width = 10)
#'
#' # Calculate specific measures globally
#' result_global <- rolling_measures(
#'   data,
#'   measures = c("complexity", "variance"),
#'   window_width = 10
#' )
#'
#' @export
#' @importFrom zoo rollapply
#' @importFrom stats acf var median

rolling_measures <- function(data,
                            id_col = NULL,
                            ts_cols = NULL,
                            measures = NULL,
                            window_width = 7L,
                            scale_min = NULL,
                            scale_max = NULL,
                            rescale = FALSE,
                            na_action = c("pass", "omit", "fail")) {
  
  # Input validation
  .validate_inputs(data, id_col, ts_cols, measures, window_width,
                   scale_min, scale_max, rescale, na_action)
  
  na_action <- match.arg(na_action)
  
  # Determine time series columns
  ts_cols <- .determine_ts_columns(data, id_col, ts_cols)
  
  # Determine measures to calculate
  measures <- .determine_measures(measures)
  
  # Handle missing values
  data <- .handle_missing_values(data, na_action)
  
  # Determine calculation scope
  use_grouping <- !is.null(id_col) && id_col %in% names(data)
  
  if (use_grouping) {
    result <- .calculate_grouped_measures(data, id_col, ts_cols, measures,
                                         window_width, scale_min, scale_max,
                                         rescale)
  } else {
    result <- .calculate_global_measures(data, ts_cols, measures,
                                        window_width, scale_min, scale_max,
                                        rescale)
  }
  
  return(result)
}

.validate_inputs <- function(data, id_col, ts_cols, measures, window_width,
                            scale_min, scale_max, rescale, na_action) {
  
  # Check data
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame", call. = FALSE)
  }
  
  if (nrow(data) == 0) {
    stop("'data' cannot be empty", call. = FALSE)
  }
  
  # Check id_col
  if (!is.null(id_col)) {
    if (!is.character(id_col) || length(id_col) != 1) {
      stop("'id_col' must be a single character string", call. = FALSE)
    }
    if (!id_col %in% names(data)) {
      stop("'id_col' must be a column name in 'data'", call. = FALSE)
    }
  }
  
  # Check ts_cols
  if (!is.null(ts_cols)) {
    if (!is.character(ts_cols)) {
      stop("'ts_cols' must be a character vector", call. = FALSE)
    }
    if (!all(ts_cols %in% names(data))) {
      missing_cols <- ts_cols[!ts_cols %in% names(data)]
      stop("Column(s) not found in data: ", paste(missing_cols, collapse = ", "),
           call. = FALSE)
    }
    if (!all(sapply(data[ts_cols], is.numeric))) {
      non_numeric <- ts_cols[!sapply(data[ts_cols], is.numeric)]
      stop("Non-numeric column(s) in 'ts_cols': ", paste(non_numeric, collapse = ", "),
           call. = FALSE)
    }
  }
  
  # Check measures
  valid_measures <- c("complexity", "fluctuation", "distribution",
                     "autocorrelation", "max", "min", "variance")
  if (!is.null(measures)) {
    if (!is.character(measures)) {
      stop("'measures' must be a character vector", call. = FALSE)
    }
    invalid_measures <- measures[!measures %in% valid_measures]
    if (length(invalid_measures) > 0) {
      stop("Invalid measure(s): ", paste(invalid_measures, collapse = ", "),
           "\nValid options: ", paste(valid_measures, collapse = ", "),
           call. = FALSE)
    }
  }
  
  # Check window_width
  if (!is.numeric(window_width) || length(window_width) != 1) {
    stop("'window_width' must be a single numeric value", call. = FALSE)
  }
  if (window_width < 2 || window_width != as.integer(window_width)) {
    stop("'window_width' must be an integer >= 2", call. = FALSE)
  }
  
  # Check scales
  if (!is.null(scale_min)) {
    if (!is.numeric(scale_min) || length(scale_min) != 1 || !is.finite(scale_min)) {
      stop("'scale_min' must be a single finite numeric value", call. = FALSE)
    }
  }
  
  if (!is.null(scale_max)) {
    if (!is.numeric(scale_max) || length(scale_max) != 1 || !is.finite(scale_max)) {
      stop("'scale_max' must be a single finite numeric value", call. = FALSE)
    }
  }
  
  if (!is.null(scale_min) && !is.null(scale_max)) {
    if (scale_min >= scale_max) {
      stop("'scale_min' must be less than 'scale_max'", call. = FALSE)
    }
  }
  
  # Check rescale
  if (!is.logical(rescale) || length(rescale) != 1) {
    stop("'rescale' must be a single logical value", call. = FALSE)
  }
  
  # Check na_action
  if (length(na_action) == 1) {
    if (!na_action %in% c("pass", "omit", "fail")) {
      stop("'na_action' must be one of: 'pass', 'omit', 'fail'", call. = FALSE)
    }
  }
}

.determine_ts_columns <- function(data, id_col, ts_cols) {
  if (is.null(ts_cols)) {
    # Auto-detect numeric columns
    exclude_cols <- if (!is.null(id_col)) id_col else character(0)
    candidate_cols <- setdiff(names(data), exclude_cols)
    ts_cols <- candidate_cols[sapply(data[candidate_cols], is.numeric)]
    
    if (length(ts_cols) == 0) {
      stop("No numeric columns found for time series analysis", call. = FALSE)
    }
    
    message("Using numeric columns: ", paste(ts_cols, collapse = ", "))
  } else {
    # Remove id_col if accidentally included
    if (!is.null(id_col) && id_col %in% ts_cols) {
      ts_cols <- setdiff(ts_cols, id_col)
      warning("Removed 'id_col' from 'ts_cols'", call. = FALSE)
    }
    
    if (length(ts_cols) == 0) {
      stop("No valid time series columns remaining", call. = FALSE)
    }
  }
  
  return(ts_cols)
}

.determine_measures <- function(measures) {
  if (is.null(measures)) {
    measures <- c("complexity", "fluctuation", "distribution",
                 "autocorrelation", "max", "min", "variance")
  }
  return(measures)
}

.handle_missing_values <- function(data, na_action) {
  if (na_action == "fail" && any(is.na(data))) {
    stop("Missing values found in data and na_action = 'fail'", call. = FALSE)
  }
  
  if (na_action == "omit") {
    original_rows <- nrow(data)
    data <- na.omit(data)
    removed_rows <- original_rows - nrow(data)
    if (removed_rows > 0) {
      message("Removed ", removed_rows, " rows with missing values")
    }
  }
  
  return(data)
}

.calculate_grouped_measures <- function(data, id_col, ts_cols, measures,
                                       window_width, scale_min, scale_max,
                                       rescale) {
  
  # Split data by groups
  groups <- split(data, data[[id_col]], drop = TRUE)
  
  results <- lapply(names(groups), function(group_name) {
    group_data <- groups[[group_name]]
    
    if (nrow(group_data) < window_width) {
      warning("Group '", group_name, "' has fewer rows (", nrow(group_data),
              ") than window_width (", window_width, "). Returning NAs.",
              call. = FALSE)
      return(.create_na_result(group_data, ts_cols, measures))
    }
    
    ts_data <- group_data[ts_cols]
    measure_results <- .calculate_measures_for_data(ts_data, measures,
                                                    window_width, scale_min,
                                                    scale_max, rescale)
    
    # Combine with original data
    cbind(group_data, measure_results)
  })
  
  # Combine all groups
  do.call(rbind, results)
}

.calculate_global_measures <- function(data, ts_cols, measures, window_width,
                                      scale_min, scale_max, rescale) {
  
  if (nrow(data) < window_width) {
    warning("Data has fewer rows (", nrow(data), ") than window_width (",
            window_width, "). Returning NAs.", call. = FALSE)
    return(.create_na_result(data, ts_cols, measures))
  }
  
  ts_data <- data[ts_cols]
  measure_results <- .calculate_measures_for_data(ts_data, measures,
                                                  window_width, scale_min,
                                                  scale_max, rescale)
  
  # Combine with original data
  cbind(data, measure_results)
}

.calculate_measures_for_data <- function(ts_data, measures, window_width,
                                        scale_min, scale_max, rescale) {
  
  # Determine scales for complexity-based measures
  scale_dependent <- c("complexity", "fluctuation", "distribution")
  needs_scales <- any(measures %in% scale_dependent)
  
  if (needs_scales) {
    if (is.null(scale_min)) {
      scale_min <- min(ts_data, na.rm = TRUE)
    }
    if (is.null(scale_max)) {
      scale_max <- max(ts_data, na.rm = TRUE)
    }
    
    if (!is.finite(scale_min) || !is.finite(scale_max) || scale_min >= scale_max) {
      warning("Invalid scale range computed from data. Using NAs for scale-dependent measures.",
              call. = FALSE)
      scale_min <- scale_max <- NA_real_
    }
  }
  
  # Calculate each measure
  results <- list()
  
  for (measure in measures) {
    measure_result <- .calculate_single_measure(ts_data, measure, window_width,
                                               scale_min, scale_max, rescale)
    
    # Add column names
    if (is.vector(measure_result)) {
      col_name <- paste0(names(ts_data)[1], "_", measure)
      results[[col_name]] <- measure_result
    } else {
      # Data frame result
      col_names <- paste0(names(ts_data), "_", measure)
      names(measure_result)[seq_along(col_names)] <- col_names
      
      # Add summary statistics for multivariate data
      if (ncol(ts_data) > 1) {
        measure_result[[paste0("mean_", measure)]] <- 
          rowMeans(measure_result[col_names], na.rm = TRUE)
        measure_result[[paste0("median_", measure)]] <- 
          apply(measure_result[col_names], 1, median, na.rm = TRUE)
        
        # Set to NA where all values are NA
        all_na <- apply(measure_result[col_names], 1, function(x) all(is.na(x)))
        measure_result[[paste0("mean_", measure)]][all_na] <- NA_real_
        measure_result[[paste0("median_", measure)]][all_na] <- NA_real_
      }
      
      for (col in names(measure_result)) {
        results[[col]] <- measure_result[[col]]
      }
    }
  }
  
  as.data.frame(results)
}

.calculate_single_measure <- function(ts_data, measure, window_width,
                                     scale_min, scale_max, rescale) {
  
  if (ncol(ts_data) == 1) {
    # Univariate case
    return(.calculate_measure_univariate(ts_data[[1]], measure, window_width,
                                        scale_min, scale_max, rescale))
  } else {
    # Multivariate case
    results <- lapply(ts_data, function(x) {
      .calculate_measure_univariate(x, measure, window_width,
                                   scale_min, scale_max, rescale)
    })
    return(as.data.frame(results))
  }
}

.calculate_measure_univariate <- function(x, measure, window_width,
                                         scale_min, scale_max, rescale) {
  
  switch(measure,
    "complexity" = .rolling_complexity(x, window_width, scale_min, scale_max, rescale),
    "fluctuation" = .rolling_fluctuation(x, window_width, scale_min, scale_max, rescale),
    "distribution" = .rolling_distribution(x, window_width, scale_min, scale_max, rescale),
    "autocorrelation" = .rolling_autocorrelation(x, window_width),
    "max" = .rolling_max(x, window_width),
    "min" = .rolling_min(x, window_width),
    "variance" = .rolling_variance(x, window_width),
    stop("Unknown measure: ", measure, call. = FALSE)
  )
}

.create_na_result <- function(data, ts_cols, measures) {
  n_rows <- nrow(data)
  results <- list()
  
  for (measure in measures) {
    for (col in ts_cols) {
      results[[paste0(col, "_", measure)]] <- rep(NA_real_, n_rows)
    }
    
    if (length(ts_cols) > 1) {
      results[[paste0("mean_", measure)]] <- rep(NA_real_, n_rows)
      results[[paste0("median_", measure)]] <- rep(NA_real_, n_rows)
    }
  }
  
  cbind(data, as.data.frame(results))
}

.rmssd <- function(x) {
  x_clean <- x[!is.na(x)]
  n <- length(x_clean)
  
  if (n < 2) {
    return(NA_real_)
  }
  
  diffs <- diff(x_clean)
  sqrt(mean(diffs^2))
}

.adjust_to_scale <- function(x, target_min, target_max) {
  if (!is.finite(target_min) || !is.finite(target_max) || target_min >= target_max) {
    return(x)
  }
  
  if (all(is.na(x))) {
    return(x)
  }
  
  x_range <- range(x, na.rm = TRUE, finite = TRUE)
  if (diff(x_range) < .Machine$double.eps^0.5) {
    return(rep(mean(c(target_min, target_max)), length(x)))
  }
  
  (x - x_range[1]) / diff(x_range) * (target_max - target_min) + target_min
}

.fluctuation_degree <- function(x, scale_min, scale_max) {
  if (!is.finite(scale_min) || !is.finite(scale_max) || scale_min >= scale_max) {
    return(NA_real_)
  }
  
  if (all(is.na(x)) || length(x[!is.na(x)]) < 2) {
    return(NA_real_)
  }
  
  f_max <- scale_max - scale_min
  f_obs <- .rmssd(x)
  
  if (is.na(f_obs)) {
    return(NA_real_)
  }
  
  f <- f_obs / f_max
  max(0, min(1, f))
}

.distribution_degree <- function(x, scale_min, scale_max) {
  if (!is.finite(scale_min) || !is.finite(scale_max) || scale_min >= scale_max) {
    return(NA_real_)
  }
  
  x_clean <- x[!is.na(x)]
  n <- length(x_clean)
  
  if (n < 2) {
    return(NA_real_)
  }
  
  uniform <- seq(from = scale_min, to = scale_max, length.out = n)
  empirical <- sort(x_clean)
  
  uni_diff <- diff(uniform)
  if (any(uni_diff <= 0)) {
    return(NA_real_)
  }
  
  emp_diff <- diff(empirical)
  deviation <- uni_diff - emp_diff
  dev_h <- deviation * (sign(deviation) + 1) / 2
  div_diff <- dev_h / uni_diff
  div_diff[is.infinite(div_diff)] <- NA
  
  D <- 1 - mean(div_diff, na.rm = TRUE)
  max(0, min(1, D))
}

.autocorr_lag1 <- function(x) {
  x_clean <- x[!is.na(x)]
  n <- length(x_clean)
  
  if (n < 2 || length(unique(x_clean)) < 2) {
    return(NA_real_)
  }
  
  tryCatch({
    acf_result <- stats::acf(x_clean, lag.max = 1, plot = FALSE, na.action = na.pass)
    if (length(acf_result$acf) >= 2) {
      return(acf_result$acf[2])
    } else {
      return(NA_real_)
    }
  }, error = function(e) {
    return(NA_real_)
  })
}

.rolling_complexity <- function(x, window_width, scale_min, scale_max, rescale) {
  if (is.na(scale_min) || is.na(scale_max)) {
    return(rep(NA_real_, length(x)))
  }
  
  fluctuation <- zoo::rollapply(
    x, 
    width = window_width,
    FUN = .fluctuation_degree,
    scale_min = scale_min,
    scale_max = scale_max,
    partial = FALSE,
    fill = NA,
    align = "right"
  )
  
  distribution <- zoo::rollapply(
    x,
    width = window_width,
    FUN = .distribution_degree,
    scale_min = scale_min,
    scale_max = scale_max,
    partial = FALSE,
    fill = NA,
    align = "right"
  )
  
  complexity <- fluctuation * distribution
  
  if (rescale) {
    complexity <- .adjust_to_scale(complexity, scale_min, scale_max)
  }
  
  complexity
}

.rolling_fluctuation <- function(x, window_width, scale_min, scale_max, rescale) {
  if (is.na(scale_min) || is.na(scale_max)) {
    return(rep(NA_real_, length(x)))
  }
  
  result <- zoo::rollapply(
    x,
    width = window_width,
    FUN = .fluctuation_degree,
    scale_min = scale_min,
    scale_max = scale_max,
    partial = FALSE,
    fill = NA,
    align = "right"
  )
  
  if (rescale) {
    result <- .adjust_to_scale(result, scale_min, scale_max)
  }
  
  result
}

.rolling_distribution <- function(x, window_width, scale_min, scale_max, rescale) {
  if (is.na(scale_min) || is.na(scale_max)) {
    return(rep(NA_real_, length(x)))
  }
  
  result <- zoo::rollapply(
    x,
    width = window_width,
    FUN = .distribution_degree,
    scale_min = scale_min,
    scale_max = scale_max,
    partial = FALSE,
    fill = NA,
    align = "right"
  )
  
  if (rescale) {
    result <- .adjust_to_scale(result, scale_min, scale_max)
  }
  
  result
}

.rolling_autocorrelation <- function(x, window_width) {
  zoo::rollapply(
    x,
    width = window_width,
    FUN = .autocorr_lag1,
    partial = FALSE,
    fill = NA,
    align = "right"
  )
}

.rolling_max <- function(x, window_width) {
  zoo::rollapply(
    x,
    width = window_width,
    FUN = function(vals) {
      if (all(is.na(vals))) return(NA_real_)
      result <- max(vals, na.rm = TRUE)
      if (is.infinite(result) && result < 0) return(NA_real_)
      result
    },
    partial = FALSE,
    fill = NA,
    align = "right"
  )
}

.rolling_min <- function(x, window_width) {
  zoo::rollapply(
    x,
    width = window_width,
    FUN = function(vals) {
      if (all(is.na(vals))) return(NA_real_)
      result <- min(vals, na.rm = TRUE)
      if (is.infinite(result) && result > 0) return(NA_real_)
      result
    },
    partial = FALSE,
    fill = NA,
    align = "right"
  )
}

#' Calculate Rolling Variance (Internal Helper)
#'
#' This internal helper function calculates the rolling variance for a time series.
#'
#' @param x A numeric vector representing the time series.
#' @param window_width Positive integer specifying the rolling window size.
#' @return A numeric vector with the calculated rolling variance values.
#' @importFrom zoo rollapply
#' @keywords internal
.rolling_variance <- function(x, window_width) {
  zoo::rollapply(
    x,
    width = window_width,
    FUN = function(vals) stats::var(vals, na.rm = TRUE),
    partial = FALSE,
    fill = NA,
    align = "right"
  )
}