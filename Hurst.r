################################################################################
# COMPLETE HURST ANALYSIS WITH EARLY WARNING SIGNAL DETECTION
# Version: 2.0
# Description: Comprehensive toolkit for Hurst exponent analysis with
#              time-varying states and graded early warning signals
################################################################################

# Main function for Hurst analysis
calculate_hurst <- function(data = NULL,
                            ts_col = NULL,
                            ts_data = NULL,
                            method = c("dfa", "rs", "mfdfa"),
                            scaling = c("none", "center", "standardize", "minmax", "iqr"),
                            min_scale = 4,
                            max_scale = NULL,
                            n_scales = 10,
                            q_values = seq(-5, 5, 0.5),
                            window_size = 50,
                            step_size = 1,
                            add_state = TRUE) {
  #' Calculate Hurst Exponent and Time-Varying Persistence States
  #'
  #' @description
  #' This function provides a comprehensive toolkit for analyzing the long-range
  #' dependence (long-term memory) of a time series via the Hurst exponent (H).
  #' It can calculate a single, global Hurst exponent for the entire series or,
  #' more powerfully, compute a time-varying Hurst exponent to identify how the
  #' persistence characteristics of the series change over time.
  #'
  #' The analysis can be performed using one of three established methods:
  #' - **DFA:** Detrended Fluctuation Analysis (robust to non-stationarity).
  #' - **R/S:** Rescaled Range Analysis (the classical method developed by Hurst).
  #' - **MFDFA:** Multifractal Detrended Fluctuation Analysis (an extension
  #'   of DFA that captures how scaling behavior varies).
  #'
  #' @section Hurst Exponent Interpretation:
  #' The Hurst exponent quantifies the nature of a time series. It is a value
  #' between 0 and 1, indicating whether the series is trending, mean-reverting,
  #' or behaving randomly.
  #'
  #' - **`H > 0.5 (Persistent)`**: A value in this range indicates a persistent
  #'   or trend-following behavior. A high value is followed by another high
  #'   value, and a low value by another low value. This is characteristic of
  #'   series with "momentum." The closer H is to 1.0, the stronger the trend.
  #'   *Example*: Stock prices that are trending upwards tend to continue
  #'   upwards in the short term.
  #'
  #' - **`H = 0.5 (Random Walk)`**: This indicates a completely random series,
  #'   like a geometric Brownian motion. There is no correlation between past
  #'   and future values. The series has no "memory."
  #'   *Example*: The path of a pollen grain in water (Brownian motion) or
  #'   the theoretical behavior of a perfectly efficient market.
  #'
  #' - **`H < 0.5 (Anti-persistent)`**: This indicates anti-persistent or
  #'   mean-reverting behavior. A high value is more likely to be followed by a
  #'   low value, and vice-versa. The series tends to revert to its long-term mean.
  #'   The closer H is to 0, the stronger the mean reversion.
  #'   *Example*: The volatility of a financial asset often exhibits
  #'   mean-reversion, where periods of high volatility are followed by
  #'   periods of low volatility.

  # Match method and scaling arguments
  method <- match.arg(method)
  scaling <- match.arg(scaling)

  # Handle input data
  use_dataframe <- FALSE
  original_data <- NULL
  ts_vector <- NULL

  if (!is.null(data) && !is.null(ts_col)) {
    # Using data frame input
    if (!is.data.frame(data)) {
      stop("'data' must be a data frame")
    }
    if (!ts_col %in% names(data)) {
      stop(paste("Column", ts_col, "not found in data"))
    }
    ts_vector <- data[[ts_col]]
    original_data <- data
    use_dataframe <- TRUE
  } else if (!is.null(ts_data)) {
    # Using vector input
    ts_vector <- ts_data
    use_dataframe <- FALSE
  } else {
    stop("Either provide 'data' and 'ts_col', or 'ts_data'")
  }

  # Validate time series
  if (!is.numeric(ts_vector)) {
    stop("Time series must be numeric")
  }

  # Apply scaling based on user input
  if (scaling != "none") {
    cat(sprintf("Applying '%s' scaling to the time series...\n", scaling))
    if (scaling == "center") {
      ts_vector <- ts_vector - mean(ts_vector, na.rm = TRUE)
    } else if (scaling == "standardize") {
      ts_vector <- scale(ts_vector, center = TRUE, scale = TRUE)[,1]
    } else if (scaling == "minmax") {
      min_val <- min(ts_vector, na.rm = TRUE)
      max_val <- max(ts_vector, na.rm = TRUE)
      ts_vector <- (ts_vector - min_val) / (max_val - min_val)
    } else if (scaling == "iqr") {
      median_val <- median(ts_vector, na.rm = TRUE)
      iqr_val <- IQR(ts_vector, na.rm = TRUE)
      if (iqr_val > 0) {
        ts_vector <- (ts_vector - median_val) / iqr_val
      } else {
        # Fallback to centering if IQR is zero
        ts_vector <- ts_vector - median_val
      }
    }
  }

  n <- length(ts_vector)

  if (!add_state) {
    # Calculate global Hurst only
    if (method == "dfa") {
      return(.calculate_hurst_dfa(ts_vector, min_scale, max_scale, n_scales))
    } else if (method == "rs") {
      return(.calculate_hurst_rs(ts_vector, min_scale, max_scale, n_scales))
    } else if (method == "mfdfa") {
      return(.calculate_multifractal_hurst(ts_vector, q_values, min_scale, max_scale, n_scales))
    }
  }

  # Calculate time-varying Hurst states

  # Validate window parameters
  if (window_size < 2 * min_scale) {
    stop(paste("Window size must be at least", 2 * min_scale))
  }
  if (window_size > n) {
    stop("Window size cannot exceed time series length")
  }
  if (step_size < 1) {
    stop("Step size must be at least 1")
  }

  # Initialize state vectors
  hurst_state <- rep(NA, n)
  r_squared_state <- rep(NA, n)

  if (method == "mfdfa") {
    # For MFDFA, also store multifractal width
    mf_width_state <- rep(NA, n)
    # Store the index of q=2 for classical Hurst
    idx_q2 <- which.min(abs(q_values - 2))
    if (length(idx_q2) == 0) {
      # If q=2 not in q_values, add it
      q_values <- sort(c(q_values, 2))
      idx_q2 <- which(q_values == 2)
    }
  }

  # Calculate states using rolling window
  # We'll calculate at positions where we have enough data
  half_window <- floor(window_size / 2)

  # Progress tracking for long series
  total_windows <- floor((n - window_size) / step_size) + 1
  progress_step <- max(1, floor(total_windows / 10))

  cat("Calculating time-varying Hurst states...\n")

  window_count <- 0
  for (i in seq(half_window + 1, n - half_window, by = step_size)) {
    # Extract window
    start_idx <- i - half_window
    end_idx <- i + half_window - 1

    if (end_idx > n) break

    window_data <- ts_vector[start_idx:end_idx]

    # Skip if too many NAs in window
    if (sum(!is.na(window_data)) < 0.8 * window_size) {
      next
    }

    # Calculate Hurst for this window
    tryCatch({
      if (method == "dfa") {
        result <- .calculate_hurst_dfa(window_data, min_scale,
                                       min(max_scale, floor(window_size/4)),
                                       min(n_scales, 5))
        hurst_state[i] <- result$hurst
        r_squared_state[i] <- result$r_squared

      } else if (method == "rs") {
        result <- .calculate_hurst_rs(window_data, min_scale,
                                      min(max_scale, floor(window_size/4)),
                                      min(n_scales, 5))
        hurst_state[i] <- result$hurst
        r_squared_state[i] <- result$r_squared

      } else if (method == "mfdfa") {
        result <- .calculate_multifractal_hurst(window_data, q_values, min_scale,
                                                min(max_scale, floor(window_size/4)),
                                                min(n_scales, 5))
        hurst_state[i] <- result$hurst  # H(q=2)
        r_squared_state[i] <- result$r_squared
        mf_width_state[i] <- result$multifractal_width
      }
    }, error = function(e) {
      # Skip this window if calculation fails
      NULL
    })

    # Progress indicator
    window_count <- window_count + 1
    if (window_count %% progress_step == 0) {
      cat(".", sep = "")
    }
  }
  cat(" Done!\n")

  # Interpolate states to fill gaps
  # First, forward/backward fill at edges
  first_valid <- which(!is.na(hurst_state))[1]
  last_valid <- tail(which(!is.na(hurst_state)), 1)

  if (!is.na(first_valid)) {
    hurst_state[1:first_valid] <- hurst_state[first_valid]
    r_squared_state[1:first_valid] <- r_squared_state[first_valid]

    hurst_state[last_valid:n] <- hurst_state[last_valid]
    r_squared_state[last_valid:n] <- r_squared_state[last_valid]

    if (method == "mfdfa") {
      mf_width_state[1:first_valid] <- mf_width_state[first_valid]
      mf_width_state[last_valid:n] <- mf_width_state[last_valid]
    }
  }

  # Linear interpolation for ALL gaps (including those from step_size > 1)
  na_pos <- which(is.na(hurst_state))
  if (length(na_pos) > 0) {
    valid_pos <- which(!is.na(hurst_state))
    if (length(valid_pos) >= 2) {
      hurst_state <- approx(valid_pos, hurst_state[valid_pos], xout = 1:n, rule = 2)$y
      r_squared_state <- approx(valid_pos, r_squared_state[valid_pos], xout = 1:n, rule = 2)$y

      if (method == "mfdfa") {
        mf_width_state <- approx(valid_pos, mf_width_state[valid_pos], xout = 1:n, rule = 2)$y
      }
    }
  }

  # Create state categories using cut() function
  breaks <- c(-Inf, 0.4, 0.5, 0.6, 0.7, Inf)
  labels <- c("strong_antipersistent", "antipersistent", "random_walk", "persistent", "strong_persistent")
  hurst_category <- cut(hurst_state, breaks = breaks, labels = labels, right = FALSE, include.lowest = TRUE)

  # Calculate state transitions
  state_changes <- c(0, diff(as.numeric(factor(hurst_category))))

  # Prepare output
  if (use_dataframe) {
    # Add columns to original data frame
    result_df <- original_data
    result_df[[paste0("hurst_", method)]] <- hurst_state
    result_df[[paste0("hurst_", method, "_r2")]] <- r_squared_state
    result_df[[paste0(method, ".state")]] <- hurst_category  # Standardized name
    result_df[[paste0("hurst_", method, "_transition")]] <- state_changes

    if (method == "mfdfa") {
      result_df[[paste0("mf_width")]] <- mf_width_state

      # Add multifractal category
      mf_category <- character(n)
      mf_category[mf_width_state < 0.1] <- "monofractal"
      mf_category[mf_width_state >= 0.1 & mf_width_state < 0.3] <- "weak_multifractal"
      mf_category[mf_width_state >= 0.3 & mf_width_state < 0.5] <- "moderate_multifractal"
      mf_category[mf_width_state >= 0.5] <- "strong_multifractal"
      result_df[[paste0("mf_category")]] <- mf_category
    }

    # Add summary attributes
    attr(result_df, "hurst_summary") <- list(
      method = method,
      window_size = window_size,
      step_size = step_size,
      mean_hurst = mean(hurst_state, na.rm = TRUE),
      sd_hurst = sd(hurst_state, na.rm = TRUE),
      min_hurst = min(hurst_state, na.rm = TRUE),
      max_hurst = max(hurst_state, na.rm = TRUE),
      n_transitions = sum(abs(state_changes) > 0)
    )

    return(result_df)

  } else {
    # Return a data frame with states
    result_df <- data.frame(
      index = 1:n,
      value = ts_vector,
      hurst = hurst_state,
      hurst_r2 = r_squared_state,
      .state = hurst_category,  # Standardized name
      hurst_transition = state_changes
    )

    if (method == "mfdfa") {
      result_df$mf_width <- mf_width_state

      mf_category <- character(n)
      mf_category[mf_width_state < 0.1] <- "monofractal"
      mf_category[mf_width_state >= 0.1 & mf_width_state < 0.3] <- "weak_multifractal"
      mf_category[mf_width_state >= 0.3 & mf_width_state < 0.5] <- "moderate_multifractal"
      mf_category[mf_width_state >= 0.5] <- "strong_multifractal"
      result_df$mf_category <- mf_category
    }

    # Add column names based on method
    names(result_df)[3:6] <- c(
      paste0("hurst_", method),
      paste0("hurst_", method, "_r2"),
      paste0(method, ".state"),
      paste0("hurst_", method, "_transition")
    )

    # Add summary attributes
    attr(result_df, "hurst_summary") <- list(
      method = method,
      window_size = window_size,
      step_size = step_size,
      mean_hurst = mean(hurst_state, na.rm = TRUE),
      sd_hurst = sd(hurst_state, na.rm = TRUE),
      min_hurst = min(hurst_state, na.rm = TRUE),
      max_hurst = max(hurst_state, na.rm = TRUE),
      n_transitions = sum(abs(state_changes) > 0)
    )

    return(result_df)
  }
}

# Internal DFA function
.calculate_hurst_dfa <- function(ts_data, min_scale = 4, max_scale = NULL, n_scales = 10) {
  # Start timing
  start_time <- Sys.time()

  # Input validation
  if (!is.numeric(ts_data)) {
    stop("Input ts_data must be a numeric vector")
  }

  if (any(is.na(ts_data))) {
    ts_data <- ts_data[!is.na(ts_data)]
  }

  n <- length(ts_data)

  if (n < 2 * min_scale) {
    stop(paste("Time series too short. Need at least", 2 * min_scale, "observations"))
  }

  if (min_scale < 4) {
    min_scale <- 4
  }

  if (is.null(max_scale)) {
    max_scale <- floor(n / 4)
  } else if (max_scale > n / 2) {
    max_scale <- floor(n / 4)
  }

  if (max_scale <= min_scale) {
    stop("max_scale must be greater than min_scale")
  }

  # Create cumulative sum
  ts_mean <- mean(ts_data)
  ts_centered <- ts_data - ts_mean
  profile <- cumsum(ts_centered)

  # Create logarithmically spaced scales
  log_scales <- seq(log(min_scale), log(max_scale), length.out = n_scales)
  scales <- unique(round(exp(log_scales)))
  scales <- scales[scales >= min_scale & scales <= max_scale]
  n_scales_actual <- length(scales)

  fluctuations <- numeric(n_scales_actual)

  # Calculate fluctuations at each scale
  for (i in 1:n_scales_actual) {
    scale <- scales[i]
    n_segments <- floor(n / scale)

    if (n_segments < 2) {
      fluctuations[i] <- NA
      next
    }

    segment_fluct <- numeric(2 * n_segments)

    # Forward direction
    for (j in 1:n_segments) {
      start_idx <- (j - 1) * scale + 1
      end_idx <- j * scale

      segment <- profile[start_idx:end_idx]
      x_vals <- 1:scale

      x_mean <- mean(x_vals)
      y_mean <- mean(segment)

      xy_cov <- sum((x_vals - x_mean) * (segment - y_mean))
      xx_var <- sum((x_vals - x_mean)^2)

      if (xx_var > 0) {
        slope <- xy_cov / xx_var
        intercept <- y_mean - slope * x_mean

        trend <- intercept + slope * x_vals
        detrended <- segment - trend

        segment_fluct[j] <- sqrt(mean(detrended^2))
      } else {
        segment_fluct[j] <- 0
      }
    }

    # Backward direction
    for (j in 1:n_segments) {
      start_idx <- n - j * scale + 1
      end_idx <- n - (j - 1) * scale

      segment <- profile[start_idx:end_idx]
      x_vals <- 1:scale

      x_mean <- mean(x_vals)
      y_mean <- mean(segment)

      xy_cov <- sum((x_vals - x_mean) * (segment - y_mean))
      xx_var <- sum((x_vals - x_mean)^2)

      if (xx_var > 0) {
        slope <- xy_cov / xx_var
        intercept <- y_mean - slope * x_mean

        trend <- intercept + slope * x_vals
        detrended <- segment - trend

        segment_fluct[n_segments + j] <- sqrt(mean(detrended^2))
      } else {
        segment_fluct[n_segments + j] <- 0
      }
    }

    fluctuations[i] <- sqrt(mean(segment_fluct^2))
  }

  # Remove NA values
  valid_idx <- !is.na(fluctuations) & fluctuations > 0
  scales_valid <- scales[valid_idx]
  fluctuations_valid <- fluctuations[valid_idx]

  if (length(scales_valid) < 3) {
    stop("Not enough valid scales")
  }

  # Calculate Hurst exponent
  log_scales <- log(scales_valid)
  log_fluct <- log(fluctuations_valid)

  x_mean <- mean(log_scales)
  y_mean <- mean(log_fluct)

  xy_cov <- sum((log_scales - x_mean) * (log_fluct - y_mean))
  xx_var <- sum((log_scales - x_mean)^2)

  slope <- xy_cov / xx_var
  intercept <- y_mean - slope * x_mean

  y_pred <- intercept + slope * log_scales
  ss_tot <- sum((log_fluct - y_mean)^2)
  ss_res <- sum((log_fluct - y_pred)^2)
  r_squared <- 1 - (ss_res / ss_tot)

  hurst <- slope

  end_time <- Sys.time()
  computation_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  output <- list(
    hurst = hurst,
    r_squared = r_squared,
    scales = scales_valid,
    fluctuations = fluctuations_valid,
    slope = slope,
    intercept = intercept,
    computation_time = computation_time,
    method = "DFA",
    n_observations = n,
    n_scales_used = length(scales_valid),
    scale_range = c(min(scales_valid), max(scales_valid))
  )

  class(output) <- c("hurst_result", "list")
  return(output)
}

# Internal R/S function
.calculate_hurst_rs <- function(ts_data, min_scale = 4, max_scale = NULL, n_scales = 10) {
  start_time <- Sys.time()

  if (!is.numeric(ts_data)) {
    stop("Input ts_data must be a numeric vector")
  }

  if (any(is.na(ts_data))) {
    ts_data <- ts_data[!is.na(ts_data)]
  }

  n <- length(ts_data)

  if (n < 2 * min_scale) {
    stop(paste("Time series too short. Need at least", 2 * min_scale, "observations"))
  }

  if (min_scale < 4) {
    min_scale <- 4
  }

  if (is.null(max_scale)) {
    max_scale <- floor(n / 4)
  } else if (max_scale > n / 2) {
    max_scale <- floor(n / 4)
  }

  if (max_scale <= min_scale) {
    stop("max_scale must be greater than min_scale")
  }

  log_scales <- seq(log(min_scale), log(max_scale), length.out = n_scales)
  scales <- unique(round(exp(log_scales)))
  scales <- scales[scales >= min_scale & scales <= max_scale]
  n_scales_actual <- length(scales)

  rs_values <- numeric(n_scales_actual)
  expected_rs <- numeric(n_scales_actual)

  for (i in 1:n_scales_actual) {
    scale <- scales[i]
    n_subseries <- floor(n / scale)

    if (n_subseries < 1) {
      rs_values[i] <- NA
      next
    }

    subseries_rs <- numeric(n_subseries)

    for (j in 1:n_subseries) {
      start_idx <- (j - 1) * scale + 1
      end_idx <- j * scale

      subseries <- ts_data[start_idx:end_idx]
      mean_val <- mean(subseries)
      adjusted <- subseries - mean_val
      cumsum_series <- cumsum(adjusted)

      R <- max(cumsum_series) - min(cumsum_series)
      S <- sqrt(mean((subseries - mean_val)^2))

      if (S > 0) {
        subseries_rs[j] <- R / S
      } else {
        subseries_rs[j] <- 0
      }
    }

    valid_rs <- subseries_rs[subseries_rs > 0]

    if (length(valid_rs) > 0) {
      rs_values[i] <- mean(valid_rs)
    } else {
      rs_values[i] <- NA
    }

    if (scale <= 340) {
      expected_rs[i] <- sqrt(scale * pi / 2) * (1 - 0.5 / scale - 0.125 / scale^2)
    } else {
      expected_rs[i] <- sqrt(scale * pi / 2)
    }
  }

  valid_idx <- !is.na(rs_values) & rs_values > 0
  scales_valid <- scales[valid_idx]
  rs_values_valid <- rs_values[valid_idx]
  expected_rs_valid <- expected_rs[valid_idx]

  if (length(scales_valid) < 3) {
    stop("Not enough valid scales")
  }

  log_scales <- log(scales_valid)
  log_rs <- log(rs_values_valid)

  x_mean <- mean(log_scales)
  y_mean <- mean(log_rs)

  xy_cov <- sum((log_scales - x_mean) * (log_rs - y_mean))
  xx_var <- sum((log_scales - x_mean)^2)

  slope <- xy_cov / xx_var
  intercept <- y_mean - slope * x_mean

  y_pred <- intercept + slope * log_scales
  ss_tot <- sum((log_rs - y_mean)^2)
  ss_res <- sum((log_rs - y_pred)^2)
  r_squared <- 1 - (ss_res / ss_tot)

  hurst <- slope

  if (hurst > 1) {
    hurst <- 1.0
  } else if (hurst < 0) {
    hurst <- 0.0
  }

  end_time <- Sys.time()
  computation_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  output <- list(
    hurst = hurst,
    r_squared = r_squared,
    scales = scales_valid,
    rs_values = rs_values_valid,
    expected_rs = expected_rs_valid,
    slope = slope,
    intercept = intercept,
    computation_time = computation_time,
    method = "R/S",
    n_observations = n,
    n_scales_used = length(scales_valid),
    scale_range = c(min(scales_valid), max(scales_valid))
  )

  class(output) <- c("hurst_result", "list")
  return(output)
}

# Internal Multifractal DFA function
.calculate_multifractal_hurst <- function(ts_data, q_values = seq(-5, 5, 0.5),
                                          min_scale = 4, max_scale = NULL, n_scales = 10) {
  start_time <- Sys.time()

  if (!is.numeric(ts_data)) {
    stop("Input ts_data must be a numeric vector")
  }

  if (any(is.na(ts_data))) {
    ts_data <- ts_data[!is.na(ts_data)]
  }

  n <- length(ts_data)

  if (n < 2 * min_scale) {
    stop(paste("Time series too short. Need at least", 2 * min_scale, "observations"))
  }

  if (min_scale < 4) {
    min_scale <- 4
  }

  if (!is.numeric(q_values)) {
    stop("q_values must be a numeric vector")
  }

  if (0 %in% q_values) {
    q_values <- q_values[q_values != 0]
  }

  if (is.null(max_scale)) {
    max_scale <- floor(n / 4)
  } else if (max_scale > n / 2) {
    max_scale <- floor(n / 4)
  }

  if (max_scale <= min_scale) {
    stop("max_scale must be greater than min_scale")
  }

  ts_mean <- mean(ts_data)
  ts_centered <- ts_data - ts_mean
  profile <- cumsum(ts_centered)

  log_scales <- seq(log(min_scale), log(max_scale), length.out = n_scales)
  scales <- unique(round(exp(log_scales)))
  scales <- scales[scales >= min_scale & scales <= max_scale]
  n_scales_actual <- length(scales)

  n_q <- length(q_values)
  fluctuations <- matrix(NA, nrow = n_scales_actual, ncol = n_q)

  for (i in 1:n_scales_actual) {
    scale <- scales[i]
    n_segments <- floor(n / scale)

    if (n_segments < 2) {
      next
    }

    segment_fluct_sq <- numeric(2 * n_segments)

    # Forward direction
    for (j in 1:n_segments) {
      start_idx <- (j - 1) * scale + 1
      end_idx <- j * scale
      
      segment <- profile[start_idx:end_idx]
      x_vals <- 1:scale

      x_mean <- mean(x_vals)
      y_mean <- mean(segment)

      xy_cov <- sum((x_vals - x_mean) * (segment - y_mean))
      xx_var <- sum((x_vals - x_mean)^2)

      if (xx_var > 0) {
        slope <- xy_cov / xx_var
        intercept <- y_mean - slope * x_mean

        trend <- intercept + slope * x_vals
        detrended <- segment - trend

        segment_fluct_sq[j] <- mean(detrended^2)
      } else {
        segment_fluct_sq[j] <- 0
      }
    }

    # Backward direction
    for (j in 1:n_segments) {
      start_idx <- n - j * scale + 1
      end_idx <- n - (j - 1) * scale

      segment <- profile[start_idx:end_idx]
      x_vals <- 1:scale

      x_mean <- mean(x_vals)
      y_mean <- mean(segment)

      xy_cov <- sum((x_vals - x_mean) * (segment - y_mean))
      xx_var <- sum((x_vals - x_mean)^2)

      if (xx_var > 0) {
        slope <- xy_cov / xx_var
        intercept <- y_mean - slope * x_mean

        trend <- intercept + slope * x_vals
        detrended <- segment - trend

        segment_fluct_sq[n_segments + j] <- mean(detrended^2)
      } else {
        segment_fluct_sq[n_segments + j] <- 0
      }
    }

    segment_fluct_sq <- segment_fluct_sq[segment_fluct_sq > 0]

    if (length(segment_fluct_sq) == 0) {
      next
    }

    for (k in 1:n_q) {
      q <- q_values[k]

      if (q > 0) {
        fluctuations[i, k] <- (mean(segment_fluct_sq^(q/2)))^(1/q)
      } else {
        fluctuations[i, k] <- exp(mean(log(sqrt(segment_fluct_sq))))
      }
    }
  }

  h_q <- numeric(n_q)
  r_squared_q <- numeric(n_q)

  for (k in 1:n_q) {
    fluct_q <- fluctuations[, k]
    valid_idx <- !is.na(fluct_q) & fluct_q > 0

    if (sum(valid_idx) < 3) {
      h_q[k] <- NA
      r_squared_q[k] <- NA
      next
    }

    scales_valid <- scales[valid_idx]
    fluct_valid <- fluct_q[valid_idx]

    log_scales <- log(scales_valid)
    log_fluct <- log(fluct_valid)

    x_mean <- mean(log_scales)
    y_mean <- mean(log_fluct)

    xy_cov <- sum((log_scales - x_mean) * (log_fluct - y_mean))
    xx_var <- sum((log_scales - x_mean)^2)

    slope <- xy_cov / xx_var

    y_pred <- y_mean - x_mean * slope + slope * log_scales
    ss_tot <- sum((log_fluct - y_mean)^2)
    ss_res <- sum((log_fluct - y_pred)^2)
    r_squared_q[k] <- 1 - (ss_res / ss_tot)

    h_q[k] <- slope
  }

  valid_h <- !is.na(h_q)
  q_valid <- q_values[valid_h]
  h_q_valid <- h_q[valid_h]
  r_squared_valid <- r_squared_q[valid_h]

  if (length(h_q_valid) < 3) {
    stop("Not enough valid H(q) values calculated")
  }

  tau_q <- q_valid * h_q_valid - 1

  n_h <- length(q_valid)
  h_alpha <- numeric(n_h - 2)

  for (i in 2:(n_h - 1)) {
    h_alpha[i-1] <- (tau_q[i+1] - tau_q[i-1]) / (q_valid[i+1] - q_valid[i-1])
  }

  q_alpha <- q_valid[2:(n_h-1)]
  D_alpha <- q_alpha * h_alpha - tau_q[2:(n_h-1)]

  multifractal_width <- max(h_q_valid) - min(h_q_valid)

  idx_q2 <- which.min(abs(q_valid - 2))
  hurst_classic <- h_q_valid[idx_q2]
  r_squared_classic <- r_squared_valid[idx_q2]

  end_time <- Sys.time()
  computation_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  output <- list(
    hurst = hurst_classic,
    r_squared = r_squared_classic,
    h_q = h_q_valid,
    q_values = q_valid,
    r_squared_q = r_squared_valid,
    tau_q = tau_q,
    h_alpha = h_alpha,
    D_alpha = D_alpha,
    multifractal_width = multifractal_width,
    scales = scales,
    computation_time = computation_time,
    method = "MFDFA",
    n_observations = n,
    n_scales_used = n_scales_actual,
    scale_range = c(min(scales), max(scales))
  )

  class(output) <- c("mfdfa_result", "hurst_result", "list")
  return(output)
}

################################################################################
# EARLY WARNING SIGNAL DETECTION FUNCTIONS
################################################################################

# Advanced Early Warning Signal Detection with Graded Warnings
detect_hurst_ews <- function(result_df,
                             method = NULL,
                             # Threshold parameters
                             extreme_percentiles = c(0.05, 0.95),
                             trend_percentile = 0.90,
                             volatility_percentile = 0.90,
                             # Window parameters
                             trend_window = 30,
                             volatility_window = 30,
                             flicker_window = 20,
                             # Advanced parameters
                             use_multivariate = TRUE,
                             use_variance_ratio = TRUE,
                             use_spectral = TRUE,
                             append_to_data = TRUE) {
  #' Advanced Early Warning Signal Detection with Graded Risk Levels
  #'
  #' @description
  #' Detects early warning signals using multiple sophisticated indicators
  #' and assigns graded warning levels (0-4) to each time point:
  #' - 0: None (no warning)
  #' - 1: Low (minor indicators present)
  #' - 2: Moderate (multiple indicators or moderate severity)
  #' - 3: High (strong indicators present)
  #' - 4: Critical (multiple strong indicators)

  # Auto-detect method
  if (is.null(method)) {
    state_cols <- grep("\\.state$", names(result_df), value = TRUE)
    if (length(state_cols) > 0) {
      method <- gsub("\\.state$", "", state_cols[1])
    } else {
      stop("No .state columns found in data")
    }
  }

  hurst_col <- paste0("hurst_", method)
  state_col <- paste0(method, ".state")

  if (!hurst_col %in% names(result_df)) {
    stop(paste("Column", hurst_col, "not found"))
  }

  hurst_values <- result_df[[hurst_col]]
  n <- length(hurst_values)

  # Initialize indicator matrices
  indicators <- matrix(0, nrow = n, ncol = 10)
  colnames(indicators) <- c("extreme_low", "extreme_high", "trend_up", "trend_down",
                            "high_volatility", "flickering", "variance_ratio",
                            "spectral_shift", "autocorr_increase", "state_persistence")

  # Calculate thresholds based on data
  thresh_low <- quantile(hurst_values, extreme_percentiles[1], na.rm = TRUE)
  thresh_high <- quantile(hurst_values, extreme_percentiles[2], na.rm = TRUE)

  # 1. Extreme values (separated into high and low)
  indicators[, "extreme_low"] <- hurst_values < thresh_low
  indicators[, "extreme_high"] <- hurst_values > thresh_high

  # 2. Trend detection with statistical significance
  if (n >= trend_window) {
    for (i in trend_window:n) {
      window_idx <- (i - trend_window + 1):i
      x <- 1:trend_window
      y <- hurst_values[window_idx]

      if (sum(!is.na(y)) >= trend_window * 0.8) {
        # Robust regression
        lm_fit <- lm(y ~ x)
        slope <- coef(lm_fit)[2]
        se <- summary(lm_fit)$coefficients[2, 2]
        t_stat <- abs(slope) / se

        # Significant trends
        if (!is.na(t_stat) && t_stat > 2) {
          trend_thresh <- quantile(abs(hurst_values - mean(hurst_values, na.rm = TRUE)),
                                   trend_percentile, na.rm = TRUE) / trend_window

          if (slope > trend_thresh) {
            indicators[i, "trend_up"] <- 1
          } else if (slope < -trend_thresh) {
            indicators[i, "trend_down"] <- 1
          }
        }
      }
    }
  }

  # 3. Volatility clustering (GARCH-like)
  if (n >= volatility_window) {
    # Calculate rolling volatility
    rolling_vol <- rep(NA, n)

    for (i in volatility_window:n) {
      window_idx <- (i - volatility_window + 1):i
      # Use squared deviations from local mean
      local_mean <- mean(hurst_values[window_idx], na.rm = TRUE)
      rolling_vol[i] <- sqrt(mean((hurst_values[window_idx] - local_mean)^2, na.rm = TRUE))
    }

    # Detect volatility clusters
    vol_thresh <- quantile(rolling_vol, volatility_percentile, na.rm = TRUE)
    indicators[, "high_volatility"] <- rolling_vol > vol_thresh
  }

  # 4. Flickering (state transitions)
  if (state_col %in% names(result_df) && n >= flicker_window) {
    states <- as.numeric(factor(result_df[[state_col]]))

    for (i in flicker_window:n) {
      window_idx <- (i - flicker_window + 1):i
      # Count unique states in window
      n_unique <- length(unique(states[window_idx]))
      # Count transitions
      n_trans <- sum(diff(states[window_idx]) != 0)

      # High flickering if many transitions relative to window size
      indicators[i, "flickering"] <- (n_trans > flicker_window * 0.2) | (n_unique >= 3)
    }
  }

  # 5. Variance ratio test (for random walk deviation)
  if (use_variance_ratio && n >= 60) {
    for (i in 60:n) {
      window_idx <- (i - 59):i
      h_vals <- hurst_values[window_idx]

      if (sum(!is.na(h_vals)) >= 50) {
        # Calculate variance ratios
        var_2 <- var(diff(h_vals), na.rm = TRUE)
        var_4 <- var(diff(h_vals, lag = 2), na.rm = TRUE) / 2
        var_8 <- var(diff(h_vals, lag = 4), na.rm = TRUE) / 4

        # Variance ratio should be ~1 for random walk
        vr_2 <- var_4 / var_2
        vr_4 <- var_8 / var_2

        # Significant deviation from 1
        indicators[i, "variance_ratio"] <- (abs(vr_2 - 1) > 0.5) | (abs(vr_4 - 1) > 0.5)
      }
    }
  }

  # 6. Spectral indicators (dominant frequency shift)
  if (use_spectral && n >= 100) {
    for (i in 100:n) {
      window_idx <- (i - 99):i
      h_vals <- hurst_values[window_idx]

      if (sum(!is.na(h_vals)) >= 90) {
        # Simple spectral analysis
        spec <- spectrum(h_vals, method = "pgram", plot = FALSE)
        # Find dominant frequency
        dom_freq_idx <- which.max(spec$spec)
        dom_period <- 1 / spec$freq[dom_freq_idx]

        # Detect shift to low frequencies (longer periods)
        indicators[i, "spectral_shift"] <- dom_period > 20
      }
    }
  }

  # 7. Autocorrelation increase
  if (n >= 50) {
    for (i in 50:n) {
      if (i >= 100) {
        # Compare recent vs past autocorrelation
        recent_idx <- (i - 49):i
        past_idx <- (i - 99):(i - 50)

        recent_acf <- acf(hurst_values[recent_idx], lag.max = 5, plot = FALSE)$acf[2:6]
        past_acf <- acf(hurst_values[past_idx], lag.max = 5, plot = FALSE)$acf[2:6]

        # Significant increase in autocorrelation
        indicators[i, "autocorr_increase"] <- mean(recent_acf) > mean(past_acf) + 0.2
      }
    }
  }

  # 8. State persistence
  if (state_col %in% names(result_df)) {
    states <- result_df[[state_col]]
    state_runs <- rle(as.character(states))

    # Mark long persistent states
    run_end_idx <- cumsum(state_runs$lengths)
    run_start_idx <- c(1, run_end_idx[-length(run_end_idx)] + 1)

    for (i in 1:length(state_runs$lengths)) {
      if (state_runs$lengths[i] > quantile(state_runs$lengths, 0.9)) {
        indicators[run_start_idx[i]:run_end_idx[i], "state_persistence"] <- 1
      }
    }
  }

  # Calculate weighted warning scores
  weights <- c(
    extreme_low = 2.0,
    extreme_high = 2.0,
    trend_up = 1.5,
    trend_down = 1.5,
    high_volatility = 1.5,
    flickering = 2.0,
    variance_ratio = 1.0,
    spectral_shift = 1.0,
    autocorr_increase = 1.5,
    state_persistence = 1.0
  )

  # Calculate raw scores
  raw_scores <- indicators %*% weights

  # Add proximity effects (warnings cluster)
  proximity_window <- 10
  enhanced_scores <- raw_scores

  for (i in (proximity_window + 1):(n - proximity_window)) {
    window_idx <- (i - proximity_window):(i + proximity_window)
    # Enhance score if neighbors also have warnings
    neighbor_effect <- mean(raw_scores[window_idx] > 0)
    enhanced_scores[i] <- raw_scores[i] * (1 + neighbor_effect)
  }

  # Convert to 5-level warning system
  warning_levels <- rep(0, n)

  # Dynamic thresholds based on score distribution
  if (any(enhanced_scores > 0)) {
    score_quantiles <- quantile(enhanced_scores[enhanced_scores > 0],
                                c(0.4, 0.6, 0.8, 0.95), na.rm = TRUE)

    warning_levels[enhanced_scores > 0 & enhanced_scores <= score_quantiles[1]] <- 1  # Low
    warning_levels[enhanced_scores > score_quantiles[1] & enhanced_scores <= score_quantiles[2]] <- 2  # Moderate
    warning_levels[enhanced_scores > score_quantiles[2] & enhanced_scores <= score_quantiles[3]] <- 3  # High
    warning_levels[enhanced_scores > score_quantiles[3]] <- 4  # Critical
  }

  # Create warning categories
  warning_categories <- factor(warning_levels,
                               levels = 0:4,
                               labels = c("none", "low", "moderate", "high", "critical"))

  # Create binary warning
  warning_binary <- warning_levels > 0

  # Calculate confidence scores (0-100)
  max_score <- max(enhanced_scores, na.rm = TRUE)
  if (max_score > 0) {
    confidence_scores <- pmin(100, (enhanced_scores / max_score) * 100)
  } else {
    confidence_scores <- rep(0, n)
  }

  # Prepare results
  if (append_to_data) {
    # Append all results to original data
    result_df$ews_level <- warning_levels
    result_df$ews_category <- warning_categories
    result_df$ews_binary <- warning_binary
    result_df$ews_confidence <- round(confidence_scores, 1)
    result_df$ews_score <- round(enhanced_scores, 3)

    # Add individual indicators
    for (i in 1:ncol(indicators)) {
      result_df[[paste0("ews_", colnames(indicators)[i])]] <- indicators[, i]
    }

    # Add attributes with summary
    attr(result_df, "ews_summary") <- list(
      method = method,
      warning_distribution = table(warning_categories),
      current_level = tail(warning_levels, 1),
      current_category = tail(warning_categories, 1),
      high_risk_periods = which(warning_levels >= 3),
      indicators_used = colnames(indicators),
      weights = weights
    )

    return(result_df)

  } else {
    # Return separate results
    results <- list(
      warning_levels = warning_levels,
      warning_categories = warning_categories,
      warning_binary = warning_binary,
      confidence_scores = confidence_scores,
      raw_scores = raw_scores,
      enhanced_scores = enhanced_scores,
      indicators = as.data.frame(indicators),
      summary = list(
        distribution = table(warning_categories),
        prop_warned = mean(warning_binary),
        mean_confidence = mean(confidence_scores[warning_binary]),
        high_risk_periods = which(warning_levels >= 3)
      )
    )

    class(results) <- c("hurst_ews", "list")
    return(results)
  }
}

################################################################################
# PLOTTING FUNCTIONS
################################################################################

# Flexible plotting function for time series and network
plot_tsn <- function(result_df, type = c("both", "series", "network"),
                     ts_col = "value", time_col = NULL, method = NULL, ...) {
  #' Flexible plotting for time series and/or network visualization

  type <- match.arg(type)

  # Auto-detect method and state column if not provided
  if (is.null(method)) {
    state_cols <- grep("\\.state$", names(result_df), value = TRUE)
    if (length(state_cols) > 0) {
      method <- gsub("\\.state$", "", state_cols[1])
    } else {
      stop("No .state columns found in data")
    }
  }

  # Get column names
  hurst_col <- paste0("hurst_", method)
  state_col <- paste0(method, ".state")

  if (!hurst_col %in% names(result_df)) {
    stop(paste("Column", hurst_col, "not found"))
  }

  if (!state_col %in% names(result_df)) {
    stop(paste("Column", state_col, "not found"))
  }

  # Determine time axis
  if (!is.null(time_col)) {
    if (!time_col %in% names(result_df)) {
      stop(paste("Time column", time_col, "not found"))
    }
    time_axis <- result_df[[time_col]]
  } else {
    time_axis <- 1:nrow(result_df)
  }

  # Handle date axes
  is_date_axis <- inherits(time_axis, "Date") || inherits(time_axis, "POSIXt")
  date_axis_params <- list()
  if (is_date_axis) {
      time_range_days <- as.numeric(difftime(max(time_axis, na.rm=TRUE),
                                             min(time_axis, na.rm=TRUE), units = "days"))
      interval <- "months"
      date_format <- "%b %Y"
      if (time_range_days > 365 * 3) {
          interval <- "6 months"
      } else if (time_range_days > 365) {
          interval <- "3 months"
      } else if (time_range_days < 90) {
          date_format <- "%d %b"
      }

      date_axis_params$at <- seq(min(time_axis, na.rm=TRUE),
                                 max(time_axis, na.rm=TRUE), by = interval)
      date_axis_params$format <- date_format
  }

  if (type == "both") {
    # Set up 2x2 layout with series plots on left, network on right
    layout(matrix(c(1,3,2,3), 2, 2, byrow = TRUE), widths = c(2,1))
    par(mar = c(4, 4, 2, 1))
  }

  if (type %in% c("both", "series")) {
    # 1. Time series with colored background by state
    plot(time_axis, result_df[[ts_col]], type = "l",
         main = "Time Series with Hurst States",
         xlab = if(is_date_axis) "" else ifelse(!is.null(time_col), time_col, "Index"),
         ylab = ts_col, xaxt = if(is_date_axis) "n" else "s", ...)

    if (is_date_axis) {
      axis.Date(1, at = date_axis_params$at, format = date_axis_params$format)
    }

    # Add colored regions
    colors <- c(
      "strong_antipersistent" = "#FF6B6B",
      "antipersistent" = "#FFA06B",
      "random_walk" = "#FFD93D",
      "persistent" = "#6BCF7F",
      "strong_persistent" = "#4ECDC4"
    )

    # Find state change points
    state_numeric <- as.numeric(factor(result_df[[state_col]]))
    change_points <- which(diff(state_numeric) != 0)
    change_points <- c(0, change_points, nrow(result_df))

    # Add colored rectangles
    for (i in 1:(length(change_points) - 1)) {
      start <- change_points[i] + 1
      end <- change_points[i + 1]
      state <- as.character(result_df[[state_col]][start])

      if (!is.na(state) && state %in% names(colors)) {
        rect(time_axis[start], par("usr")[3], time_axis[end], par("usr")[4],
             col = adjustcolor(colors[state], alpha = 0.3),
             border = NA)
      }
    }

    # Redraw the time series
    lines(time_axis, result_df[[ts_col]], col = "black", lwd = 1.5)

    # 2. Hurst exponent over time
    plot(time_axis, result_df[[hurst_col]], type = "l",
         main = paste("Hurst Exponent (", toupper(method), ")", sep = ""),
         xlab = if(is_date_axis) "" else ifelse(!is.null(time_col), time_col, "Index"),
         ylab = "Hurst Exponent",
         ylim = c(0, 1), col = "blue", lwd = 2,
         xaxt = if(is_date_axis) "n" else "s", ...)

    if (is_date_axis) {
      axis.Date(1, at = date_axis_params$at, format = date_axis_params$format)
    }

    # Add reference lines
    abline(h = 0.5, lty = 2, col = "gray")
    abline(h = c(0.4, 0.6), lty = 3, col = "gray")

    # Add state regions as background
    for (i in 1:(length(change_points) - 1)) {
      start <- change_points[i] + 1
      end <- change_points[i + 1]
      state <- as.character(result_df[[state_col]][start])

      if (!is.na(state) && state %in% names(colors)) {
        rect(time_axis[start], 0, time_axis[end], 1,
             col = adjustcolor(colors[state], alpha = 0.15),
             border = NA)
      }
    }

    # Redraw the Hurst line
    lines(time_axis, result_df[[hurst_col]], col = "blue", lwd = 2)
  }

  if (type %in% c("both", "network")) {
    if (type == "network") {
      par(mar = c(2, 2, 2, 2))
    }

    # 3. TNA Network plot
    if (requireNamespace("tna", quietly = TRUE)) {
      # Prepare data for TNA
      tna_data <- tna::prepare_data(result_df, action = state_col)
      tna_series <- tna::tna(tna_data)

      # Create the network plot
      plot(tna_series, scale_nodes = "Instrength",
           main = "State Transition Network")
    } else {
      # Fallback to transition matrix heatmap if tna not available
      state_transitions <- table(
        From = result_df[[state_col]][-nrow(result_df)],
        To = result_df[[state_col]][-1]
      )

      # Normalize to probabilities
      trans_prob <- prop.table(state_transitions, margin = 1)

      # Create heatmap-style plot
      image(as.matrix(trans_prob),
            main = "Transition Probabilities\n(install 'tna' for network plot)",
            xlab = "From State", ylab = "To State",
            col = heat.colors(20))

      # Add grid
      grid(nx = nrow(trans_prob), ny = ncol(trans_prob), col = "white", lty = 1)
    }
  }

  if (type == "both") {
    layout(1)  # Reset layout
  }
  par(mar = c(5, 4, 4, 2) + 0.1)  # Reset margins
}

# Enhanced plotting function for graded warnings
plot_ews_graded <- function(result_df, time_col = NULL,
                            show_indicators = TRUE, ...) {
  #' Plot graded early warning signals

  # Check for EWS columns
  if (!"ews_level" %in% names(result_df)) {
    stop("No EWS data found. Run detect_hurst_ews() with append_to_data = TRUE")
  }

  # Determine time axis
  if (!is.null(time_col) && time_col %in% names(result_df)) {
    time_axis <- result_df[[time_col]]
  } else {
    time_axis <- 1:nrow(result_df)
  }

  # Color scheme for warning levels
  warning_colors <- c(
    "0" = "#2ECC71",  # Green - None
    "1" = "#F1C40F",  # Yellow - Low
    "2" = "#E67E22",  # Orange - Moderate
    "3" = "#E74C3C",  # Red - High
    "4" = "#8E44AD"   # Purple - Critical
  )

  # Set up plot layout
  if (show_indicators) {
    layout(matrix(c(1,1,2,2,3,3,3), ncol = 1))
    par(mar = c(2, 4, 2, 2))
  } else {
    par(mfrow = c(2, 1))
    par(mar = c(4, 4, 2, 2))
  }

  # 1. Time series with warning level background
  # Find Hurst column
  hurst_cols <- grep("^hurst_", names(result_df), value = TRUE)
  hurst_col <- hurst_cols[!grepl("_r2|_transition", hurst_cols)][1]

  plot(time_axis, result_df[[hurst_col]], type = "n",
       main = "Hurst Exponent with Graded Warning Levels",
       xlab = "", ylab = "Hurst Exponent", ylim = c(0, 1))

  # Add colored background for warning levels
  warning_changes <- c(0, which(diff(result_df$ews_level) != 0), nrow(result_df))

  for (i in 1:(length(warning_changes) - 1)) {
    start_idx <- warning_changes[i] + 1
    end_idx <- warning_changes[i + 1]
    level <- as.character(result_df$ews_level[start_idx])

    rect(time_axis[start_idx], 0, time_axis[end_idx], 1,
         col = adjustcolor(warning_colors[level], alpha = 0.3),
         border = NA)
  }

  # Add Hurst line
  lines(time_axis, result_df[[hurst_col]], col = "black", lwd = 2)
  abline(h = 0.5, lty = 2, col = "gray")

  # 2. Warning level time series
  plot(time_axis, result_df$ews_level, type = "s", lwd = 2,
       col = "black", ylim = c(-0.5, 4.5),
       main = "Early Warning Signal Level",
       xlab = ifelse(show_indicators, "",
                      ifelse(!is.null(time_col), time_col, "Index")),
       ylab = "Warning Level", yaxt = "n")

  # Add colored points
  points(time_axis, result_df$ews_level,
         col = warning_colors[as.character(result_df$ews_level)],
         pch = 19, cex = 0.8)

  # Custom y-axis
  axis(2, at = 0:4, labels = c("None", "Low", "Moderate", "High", "Critical"))

  # Add confidence as background
  if ("ews_confidence" %in% names(result_df)) {
    # Add semi-transparent bars for confidence
    for (i in 1:nrow(result_df)) {
      if (result_df$ews_level[i] > 0) {
        rect(time_axis[i] - 0.4, -0.5, time_axis[i] + 0.4,
             result_df$ews_confidence[i]/100 * 4.5 - 0.5,
             col = adjustcolor("gray", alpha = 0.2), border = NA)
      }
    }
  }

  # 3. Individual indicators (if requested)
  if (show_indicators) {
    # Get indicator columns
    indicator_cols <- grep("^ews_", names(result_df), value = TRUE)
    indicator_cols <- indicator_cols[!indicator_cols %in%
                                     c("ews_level", "ews_category", "ews_binary",
                                       "ews_confidence", "ews_score")]

    if (length(indicator_cols) > 0) {
      # Create heatmap of indicators
      indicator_matrix <- t(as.matrix(result_df[, indicator_cols]))

      # Replace column names for better display
      nice_names <- gsub("ews_", "", indicator_cols)
      nice_names <- gsub("_", " ", nice_names)
      nice_names <- tools::toTitleCase(nice_names)
      
      image(x = time_axis, y = 1:length(indicator_cols),
            z = t(indicator_matrix),
            col = c("white", "darkred"),
            main = "Individual Warning Indicators",
            xlab = ifelse(!is.null(time_col), time_col, "Index"),
            ylab = "", yaxt = "n")

      axis(2, at = 1:length(indicator_cols), labels = nice_names, las = 2)

      # Add grid
      abline(h = seq(0.5, length(indicator_cols) + 0.5, 1),
             col = "gray", lty = 3)
    }
  }

  # Reset layout
  par(mfrow = c(1, 1))
  layout(1)
}

# Separate time series plotting function
plot_ts <- function(data, ts_col = "value", time_col = NULL,
                    state_col = NULL, main = NULL, ...) {
  #' Plot time series with optional state coloring

  # Validate inputs
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }

  if (!ts_col %in% names(data)) {
    stop(paste("Column", ts_col, "not found in data"))
  }

  # Determine time axis
  if (!is.null(time_col)) {
    if (!time_col %in% names(data)) {
      stop(paste("Time column", time_col, "not found"))
    }
    time_axis <- data[[time_col]]
    xlab <- time_col
  } else {
    time_axis <- 1:nrow(data)
    xlab <- "Index"
  }

  is_date_axis <- inherits(time_axis, "Date") || inherits(time_axis, "POSIXt")

  # Set default title
  if (is.null(main)) {
    main <- paste("Time Series:", ts_col)
  }

  # Create base plot
  plot(time_axis, data[[ts_col]], type = "l",
       main = main, xlab = if(is_date_axis) "" else xlab, ylab = ts_col,
       xaxt = if(is_date_axis) "n" else "s", ...)

  # Add custom date axis if needed
  if (is_date_axis) {
      time_range_days <- as.numeric(difftime(max(time_axis, na.rm=TRUE),
                                             min(time_axis, na.rm=TRUE), units = "days"))
      interval <- "months"
      date_format <- "%b %Y"
      if (time_range_days > 365 * 3) {
          interval <- "6 months"
      } else if (time_range_days > 365) {
          interval <- "3 months"
      } else if (time_range_days < 90) {
          date_format <- "%d %b %Y"
      }

      axis.Date(1, at = seq(min(time_axis, na.rm=TRUE),
                           max(time_axis, na.rm=TRUE), by = interval),
                format = date_format)
  }

  # Add state coloring if specified
  if (!is.null(state_col) && state_col %in% names(data)) {
    # Define colors
    colors <- c(
      "strong_antipersistent" = "#FF6B6B",
      "antipersistent" = "#FFA06B",
      "random_walk" = "#FFD93D",
      "persistent" = "#6BCF7F",
      "strong_persistent" = "#4ECDC4",
      "monofractal" = "#E8E8E8",
      "weak_multifractal" = "#C8E6C9",
      "moderate_multifractal" = "#81C784",
      "strong_multifractal" = "#4CAF50"
    )

    # Find state change points
    state_factor <- factor(data[[state_col]])
    state_numeric <- as.numeric(state_factor)
    change_points <- which(diff(state_numeric) != 0)
    change_points <- c(0, change_points, nrow(data))

    # Add colored regions
    for (i in 1:(length(change_points) - 1)) {
      start <- change_points[i] + 1
      end <- change_points[i + 1]
      state <- as.character(data[[state_col]][start])

      if (state %in% names(colors)) {
        rect(time_axis[start], par("usr")[3], time_axis[end], par("usr")[4],
             col = adjustcolor(colors[state], alpha = 0.3),
             border = NA)
      }
    }

    # Redraw the time series
    lines(time_axis, data[[ts_col]], col = "black", lwd = 1.5)

    # Add legend
    unique_states <- unique(as.character(data[[state_col]]))
    legend_colors <- colors[unique_states]
    legend_colors <- legend_colors[!is.na(legend_colors)]

    if (length(legend_colors) > 0) {
      legend("topright",
             legend = names(legend_colors),
             fill = adjustcolor(legend_colors, alpha = 0.3),
             cex = 0.8)
    }
  }
}

# Plotting functions for global Hurst results
plot_hurst_scaling <- function(result, main = NULL, ...) {
  if (!inherits(result, "hurst_result")) {
    stop("Input must be a hurst_result object")
  }

  if (is.null(main)) {
    main <- paste(result$method, ": Hurst =", round(result$hurst, 3),
                  ", R² =", round(result$r_squared, 3))
  }

  if (result$method == "DFA") {
    xlab <- "log(Scale)"
    ylab <- "log(Fluctuation)"
    x_data <- log(result$scales)
    y_data <- log(result$fluctuations)
  } else if (result$method == "R/S") {
    xlab <- "log(Scale)"
    ylab <- "log(R/S)"
    x_data <- log(result$scales)
    y_data <- log(result$rs_values)
  } else {
    stop("Plotting not implemented for MFDFA. Use plot_mfdfa() instead.")
  }

  plot(x_data, y_data,
       xlab = xlab, ylab = ylab,
       main = main,
       pch = 19, col = "blue", ...)

  if (result$method == "DFA") {
    abline(a = result$intercept, b = result$slope, col = "red", lwd = 2)
  } else if (result$method == "R/S") {
    abline(a = result$intercept, b = result$slope, col = "red", lwd = 2)

    if (!is.null(result$expected_rs)) {
      points(log(result$scales), log(result$expected_rs),
             pch = 4, col = "green")
      legend("bottomright",
             legend = c("Observed", "Fitted", "Expected (H=0.5)"),
             col = c("blue", "red", "green"),
             pch = c(19, NA, 4),
             lty = c(NA, 1, NA),
             lwd = c(NA, 2, NA))
    }
  }

  grid()
}

plot_mfdfa <- function(result, type = c("hq", "spectrum", "both"), ...) {
  if (!inherits(result, "mfdfa_result")) {
    stop("Input must be an mfdfa_result object")
  }

  type <- match.arg(type)

  if (type == "both") {
    par(mfrow = c(1, 2))
  }

  if (type %in% c("hq", "both")) {
    plot(result$q_values, result$h_q,
         xlab = "q", ylab = "H(q)",
         main = paste("Generalized Hurst Exponents\nWidth =",
                      round(result$multifractal_width, 3)),
         pch = 19, col = "blue", type = "b", ...)

    abline(h = result$hurst, col = "red", lty = 2)
    text(min(result$q_values), result$hurst + 0.02,
         "H(q=2)", pos = 4, col = "red")

    grid()
  }

  if (type %in% c("spectrum", "both")) {
    plot(result$h_alpha, result$D_alpha,
         xlab = expression(alpha), ylab = expression(f(alpha)),
         main = "Singularity Spectrum",
         pch = 19, col = "blue", type = "b", ...)

    max_idx <- which.max(result$D_alpha)
    points(result$h_alpha[max_idx], result$D_alpha[max_idx],
           col = "red", pch = 19, cex = 1.5)

    grid()
  }

  if (type == "both") {
    par(mfrow = c(1, 1))
  }
}

################################################################################
# PRINT METHODS
################################################################################

print.hurst_result <- function(x, ...) {
  cat("Hurst Exponent Analysis (", x$method, ")\n", sep = "")
  cat("=====================================\n")
  cat("Hurst exponent: ", sprintf("%.4f", x$hurst), "\n", sep = "")
  cat("R-squared:      ", sprintf("%.4f", x$r_squared), "\n", sep = "")
  cat("Interpretation: ")

  if (x$hurst < 0.4) {
    cat("Strong anti-persistence (mean-reverting)\n")
  } else if (x$hurst < 0.5) {
    cat("Anti-persistent (mean-reverting)\n")
  } else if (x$hurst < 0.6) {
    cat("Random walk (no memory)\n")
  } else if (x$hurst < 0.7) {
    cat("Persistent (trending)\n")
  } else {
    cat("Strong persistence (trending)\n")
  }

  cat("\nAnalysis details:\n")
  cat("- Time series length: ", x$n_observations, "\n", sep = "")
  cat("- Scales analyzed:    ", x$n_scales_used, "\n", sep = "")
  cat("- Scale range:        [", x$scale_range[1], ", ", x$scale_range[2], "]\n", sep = "")
  cat("- Computation time:   ", sprintf("%.3f", x$computation_time), " seconds\n", sep = "")
}

print.mfdfa_result <- function(x, ...) {
  cat("Multifractal DFA Analysis\n")
  cat("=====================================\n")
  cat("Classical Hurst (q=2): ", sprintf("%.4f", x$hurst), "\n", sep = "")
  cat("R-squared (q=2):       ", sprintf("%.4f", x$r_squared), "\n", sep = "")
  cat("Multifractal width:    ", sprintf("%.4f", x$multifractal_width), "\n", sep = "")

  cat("\nMultifractal properties:\n")
  if (x$multifractal_width < 0.1) {
    cat("- Monofractal (uniform scaling)\n")
  } else if (x$multifractal_width < 0.3) {
    cat("- Weakly multifractal\n")
  } else if (x$multifractal_width < 0.5) {
    cat("- Moderately multifractal\n")
  } else {
    cat("- Strongly multifractal\n")
  }

  cat("\nH(q) range: [", sprintf("%.3f", min(x$h_q)), ", ",
      sprintf("%.3f", max(x$h_q)), "]\n", sep = "")
  cat("q range:    [", min(x$q_values), ", ", max(x$q_values), "]\n", sep = "")

  cat("\nAnalysis details:\n")
  cat("- Time series length: ", x$n_observations, "\n", sep = "")
  cat("- Scales analyzed:    ", x$n_scales_used, "\n", sep = "")
  cat("- Scale range:        [", x$scale_range[1], ", ", x$scale_range[2], "]\n", sep = "")
  cat("- Computation time:   ", sprintf("%.3f", x$computation_time), " seconds\n", sep = "")
}

################################################################################
# UTILITY FUNCTIONS
################################################################################

# Legacy wrapper for backward compatibility
plot_hurst_states <- function(result_df, ts_col = "value", time_col = NULL, method = NULL, ...) {
  plot_tsn(result_df, type = "both", ts_col = ts_col, time_col = time_col, method = method, ...)
}

# Get state summary statistics
get_hurst_summary <- function(result_df, method = NULL) {
  #' Extract summary statistics for Hurst states

  # Auto-detect method if not provided
  if (is.null(method)) {
    state_cols <- grep("\\.state$", names(result_df), value = TRUE)
    if (length(state_cols) > 0) {
      method <- gsub("\\.state$", "", state_cols[1])
    } else {
      stop("No .state columns found in data")
    }
  }

  hurst_col <- paste0("hurst_", method)
  state_col <- paste0(method, ".state")

  if (!hurst_col %in% names(result_df)) {
    stop(paste("Column", hurst_col, "not found"))
  }

  # Calculate state durations
  state_runs <- rle(as.character(result_df[[state_col]]))

  # State transition matrix
  state_transitions <- table(
    From = result_df[[state_col]][-nrow(result_df)],
    To = result_df[[state_col]][-1]
  )

  # Normalize to get transition probabilities
  trans_prob <- prop.table(state_transitions, margin = 1)

  summary_list <- list(
    method = method,
    basic_stats = list(
      mean_hurst = mean(result_df[[hurst_col]], na.rm = TRUE),
      sd_hurst = sd(result_df[[hurst_col]], na.rm = TRUE),
      min_hurst = min(result_df[[hurst_col]], na.rm = TRUE),
      max_hurst = max(result_df[[hurst_col]], na.rm = TRUE),
      median_hurst = median(result_df[[hurst_col]], na.rm = TRUE)
    ),
    state_distribution = table(result_df[[state_col]]),
    state_proportions = prop.table(table(result_df[[state_col]])),
    state_durations = list(
      mean_duration = mean(state_runs$lengths),
      max_duration = max(state_runs$lengths),
      min_duration = min(state_runs$lengths),
      duration_by_state = tapply(state_runs$lengths, state_runs$values, mean)
    ),
    transition_matrix = state_transitions,
    transition_probabilities = trans_prob,
    n_transitions = sum(diff(as.numeric(factor(result_df[[state_col]]))) != 0)
  )

  class(summary_list) <- c("hurst_summary", "list")
  return(summary_list)
}

# Print method for hurst_summary
print.hurst_summary <- function(x, ...) {
  cat("Hurst State Analysis Summary (", toupper(x$method), ")\n", sep = "")
  cat("=========================================\n\n")

  cat("Basic Statistics:\n")
  cat(sprintf("  Mean Hurst:   %.4f\n", x$basic_stats$mean_hurst))
  cat(sprintf("  SD Hurst:     %.4f\n", x$basic_stats$sd_hurst))
  cat(sprintf("  Range:        [%.4f, %.4f]\n", x$basic_stats$min_hurst, x$basic_stats$max_hurst))

  cat("\nState Distribution:\n")
  print(x$state_proportions)

  cat("\nAverage State Durations:\n")
  print(round(x$state_durations$duration_by_state, 1))

  cat("\nTransition Probabilities:\n")
  print(round(x$transition_probabilities, 3))

  cat(sprintf("\nTotal Transitions: %d\n", x$n_transitions))
}

# Export state sequences for further analysis
export_state_sequence <- function(result_df, method = NULL, format = c("vector", "rle", "dataframe")) {
  #' Export state sequence in various formats

  format <- match.arg(format)

  # Auto-detect method if not provided
  if (is.null(method)) {
    state_cols <- grep("\\.state$", names(result_df), value = TRUE)
    if (length(state_cols) > 0) {
      method <- gsub("\\.state$", "", state_cols[1])
    } else {
      stop("No .state columns found in data")
    }
  }

  state_col <- paste0(method, ".state")

  if (!state_col %in% names(result_df)) {
    stop(paste("Column", state_col, "not found"))
  }

  states <- result_df[[state_col]]

  if (format == "vector") {
    return(states)
  } else if (format == "rle") {
    return(rle(as.character(states)))
  } else if (format == "dataframe") {
    state_runs <- rle(as.character(states))
    df <- data.frame(
      state = state_runs$values,
      start = cumsum(c(1, state_runs$lengths[-length(state_runs$lengths)])),
      end = cumsum(state_runs$lengths),
      duration = state_runs$lengths
    )
    return(df)
  }
}

# Compare multiple methods
compare_hurst_methods <- function(data, ts_col, methods = c("dfa", "rs"),
                                  window_size = 50, step_size = 1) {
  #' Compare Hurst estimates across multiple methods

  result_df <- data

  for (method in methods) {
    cat(sprintf("Calculating %s...\n", toupper(method)))

    if (method == "mfdfa") {
      temp_result <- calculate_hurst(data = data, ts_col = ts_col,
                                     method = method,
                                     window_size = window_size,
                                     step_size = step_size,
                                     q_values = seq(-3, 3, 1))
    } else {
      temp_result <- calculate_hurst(data = data, ts_col = ts_col,
                                     method = method,
                                     window_size = window_size,
                                     step_size = step_size)
    }

    # Extract only the new columns
    new_cols <- setdiff(names(temp_result), names(result_df))
    result_df <- cbind(result_df, temp_result[, new_cols, drop = FALSE])
  }

  return(result_df)
}

# Summary function for EWS results
summarize_ews <- function(result_df) {
  #' Summarize early warning signal results

  if (!"ews_level" %in% names(result_df)) {
    stop("No EWS data found in result_df")
  }

  # Get EWS summary from attributes if available
  ews_summary <- attr(result_df, "ews_summary")

  # Calculate additional summaries
  summary_list <- list(
    # Distribution
    level_distribution = table(result_df$ews_level),
    category_distribution = table(result_df$ews_category),

    # Time in each level
    proportion_by_level = prop.table(table(result_df$ews_level)),

    # Transitions
    level_transitions = table(
      From = result_df$ews_level[-nrow(result_df)],
      To = result_df$ews_level[-1]
    ),

    # High risk analysis
    n_high_risk = sum(result_df$ews_level >= 3),
    prop_high_risk = mean(result_df$ews_level >= 3),

    # Current status
    current_level = tail(result_df$ews_level, 1),
    current_category = tail(result_df$ews_category, 1),
    current_confidence = tail(result_df$ews_confidence, 1),

    # Trend
    recent_trend = ifelse(
      mean(tail(result_df$ews_level, 20)) > mean(head(tail(result_df$ews_level, 40), 20)),
      "increasing", "stable/decreasing"
    )
  )

  # Add method info if available
  if (!is.null(ews_summary)) {
    summary_list$method_info = ews_summary
  }

  class(summary_list) <- c("ews_summary", "list")
  return(summary_list)
}

# Print method for EWS summary
print.ews_summary <- function(x, ...) {
  cat("Early Warning Signal Summary\n")
  cat("===========================\n\n")

  cat("Current Status:\n")
  cat(sprintf("  Level: %d (%s)\n", x$current_level, x$current_category))
  cat(sprintf("  Confidence: %.1f%%\n", x$current_confidence))
  cat(sprintf("  Recent trend: %s\n", x$recent_trend))

  cat("\nWarning Level Distribution:\n")
  print(x$proportion_by_level)

  cat(sprintf("\nHigh Risk Periods: %d (%.1f%% of time)\n",
              x$n_high_risk, x$prop_high_risk * 100))

  cat("\nLevel Transition Matrix:\n")
  print(x$level_transitions)
}

################################################################################
# EXAMPLE USAGE FUNCTIONS
################################################################################

# Example with EWS detection
example_complete_analysis <- function() {
  # Generate test data with regime changes
  set.seed(123)
  n <- 1000

  # Create regime-switching data
  regime1 <- arima.sim(list(ar = 0.5), n = 300)
  transition <- arima.sim(list(ar = c(0.8, -0.3)), n = 200)  # Unstable
  regime2 <- arima.sim(list(ar = -0.4), n = 300)
  critical <- rnorm(200, sd = 2)  # High volatility

  ts_data <- cumsum(c(regime1, transition, regime2, critical))

  df <- data.frame(
    date = seq(as.Date("2020-01-01"), length.out = n, by = "day"),
    value = ts_data
  )

  # Step 1: Calculate Hurst states
  cat("Step 1: Calculating Hurst states...\n")
  df_hurst <- calculate_hurst(data = df, ts_col = "value", method = "dfa",
                              window_size = 100, step_size = 5)

  # Step 2: Detect early warning signals
  cat("\nStep 2: Detecting early warning signals...\n")
  df_ews <- detect_hurst_ews(df_hurst, append_to_data = TRUE)

  # Step 3: Summarize results
  cat("\nStep 3: Summary of results\n")
  hurst_summary <- get_hurst_summary(df_ews)
  print(hurst_summary)

  ews_summary <- summarize_ews(df_ews)
  print(ews_summary)

  # Step 4: Visualize
  cat("\nStep 4: Creating visualizations...\n")

  # Plot Hurst states with network
  plot_tsn(df_ews, type = "both", ts_col = "value", time_col = "date")

  # Plot EWS with indicators
  plot_ews_graded(df_ews, time_col = "date", show_indicators = TRUE)

  # Show critical periods
  critical_periods <- which(df_ews$ews_level == 4)
  if (length(critical_periods) > 0) {
    cat("\nCritical warning periods detected at:\n")
    print(df_ews$date[critical_periods])
  }

  return(df_ews)
}

# To run the example and see the toolkit in action:
# result <- example_complete_analysis()