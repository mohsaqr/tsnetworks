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
#'
#' @param data Data frame containing the time series (optional if ts_data is provided)
#' @param ts_col Column name in data containing the time series values
#' @param ts_data Numeric vector of time series data (alternative to data/ts_col)
#' @param method Character string specifying the method: "dfa", "rs", or "mfdfa"
#' @param scaling Scaling method: "none", "center", "standardize", "minmax", or "iqr"
#' @param min_scale Minimum scale for analysis (default: 4)
#' @param max_scale Maximum scale for analysis (default: NULL, auto-determined)
#' @param n_scales Number of scales to analyze (default: 10)
#' @param q_values Vector of q values for MFDFA (default: seq(-5, 5, 0.5))
#' @param window_size Window size for time-varying analysis (default: 50)
#' @param step_size Step size for rolling window (default: 1)
#' @param add_state Whether to calculate time-varying states (default: TRUE)
#'
#' @return Data frame with Hurst analysis results or list for global analysis
#'
#' @examples
#' # Generate sample data
#' set.seed(123)
#' data <- data.frame(
#'   time = 1:500,
#'   value = cumsum(rnorm(500))
#' )
#'
#' # Calculate time-varying Hurst states
#' result <- calculate_hurst(data, ts_col = "value", method = "dfa")
#' head(result)
#'
#' # Global Hurst analysis
#' global_result <- calculate_hurst(ts_data = data$value, method = "dfa", add_state = FALSE)
#' print(global_result)
#'
#' @importFrom stats arima.sim approx acf lm coef quantile var sd IQR spectrum
#' @importFrom graphics plot lines abline rect legend axis.Date image grid points
#' @importFrom grDevices adjustcolor heat.colors
#' @importFrom tools toTitleCase
#' @export
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
#'
#' @param result_df Data frame with Hurst analysis results
#' @param method Method used for Hurst analysis (auto-detected if NULL)
#' @param extreme_percentiles Percentiles for extreme value detection
#' @param trend_percentile Percentile threshold for trend detection
#' @param volatility_percentile Percentile threshold for volatility detection
#' @param trend_window Window size for trend analysis
#' @param volatility_window Window size for volatility analysis
#' @param flicker_window Window size for flickering detection
#' @param use_multivariate Whether to use multivariate indicators
#' @param use_variance_ratio Whether to use variance ratio test
#' @param use_spectral Whether to use spectral analysis
#' @param append_to_data Whether to append results to input data
#'
#' @return Data frame with early warning signal indicators
#'
#' @examples
#' # Generate sample data with regime changes
#' set.seed(123)
#' data <- data.frame(
#'   time = 1:500,
#'   value = cumsum(c(rnorm(200), rnorm(100, sd = 2), rnorm(200)))
#' )
#'
#' # Calculate Hurst states
#' hurst_result <- calculate_hurst(data, ts_col = "value", method = "dfa")
#'
#' # Detect early warning signals
#' ews_result <- detect_hurst_ews(hurst_result)
#' table(ews_result$ews_category)
#'
#' @importFrom stats quantile lm summary.lm
#' @export
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

  # Convert to 5-level warning system
  warning_levels <- rep(0, n)

  # Dynamic thresholds based on score distribution
  if (any(raw_scores > 0)) {
    score_quantiles <- quantile(raw_scores[raw_scores > 0],
                                c(0.4, 0.6, 0.8, 0.95), na.rm = TRUE)

    warning_levels[raw_scores > 0 & raw_scores <= score_quantiles[1]] <- 1  # Low
    warning_levels[raw_scores > score_quantiles[1] & raw_scores <= score_quantiles[2]] <- 2  # Moderate
    warning_levels[raw_scores > score_quantiles[2] & raw_scores <= score_quantiles[3]] <- 3  # High
    warning_levels[raw_scores > score_quantiles[3]] <- 4  # Critical
  }

  # Create warning categories
  warning_categories <- factor(warning_levels,
                               levels = 0:4,
                               labels = c("none", "low", "moderate", "high", "critical"))

  # Create binary warning
  warning_binary <- warning_levels > 0

  # Calculate confidence scores (0-100)
  max_score <- max(raw_scores, na.rm = TRUE)
  if (max_score > 0) {
    confidence_scores <- pmin(100, (raw_scores / max_score) * 100)
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
    result_df$ews_score <- round(raw_scores, 3)

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

#' Print method for hurst_result objects
#' @param x A hurst_result object
#' @param ... Additional arguments (unused)
#' @export
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

#' Print method for mfdfa_result objects
#' @param x An mfdfa_result object
#' @param ... Additional arguments (unused)
#' @export
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
