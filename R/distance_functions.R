# Native Distance Function Implementations
# Updated to provide native implementations without external dependencies

# Time Series Distance Functions
# Individual distance calculation methods

#' Jitter spike train for significance testing
#' @param spike_train Numeric vector of spike times
#' @param jitter_val Numeric. The range for jittering.
#' @param series_length_for_jitter Numeric. The maximum possible time point in the series.
#' @return Jittered spike train
#' @keywords internal
.jitter_spike_train <- function(spike_train, jitter_val, series_length_for_jitter = NULL) {
  if (length(spike_train) == 0) return(numeric(0))
  
  jittered_train <- spike_train + stats::runif(length(spike_train), -jitter_val, jitter_val)
  
  # Ensure jittered spikes stay within reasonable bounds if series_length_for_jitter is provided
  if (!is.null(series_length_for_jitter)) {
    jittered_train <- pmax(0, jittered_train) # Spikes cannot be negative
    jittered_train <- pmin(series_length_for_jitter, jittered_train) # Spikes cannot exceed series length
  }
  
  return(sort(jittered_train))
}

dist_cor <- function(series1, series2, method = "pearson", type = "absolute", significance_test = TRUE) {
  series1 <- validate_time_series(series1, "series1")
  series2 <- validate_time_series(series2, "series2")
  .ensure_same_length(series1, series2, "dist_cor: ")
  
  if(length(series1) < 3) {
    warning("dist_cor: Fewer than 3 observations. Returning NA.")
    return(list(correlation = NA, p_value = if(significance_test) NA else NULL))
  }
  
  test_result <- stats::cor.test(series1, series2, method = method, exact = FALSE, na.action = "na.omit")
  correlation <- test_result$estimate
  
  result_cor <- switch(type, 
                      absolute = abs(correlation), 
                      positive = max(0, correlation, na.rm = TRUE),
                      negative = min(0, correlation, na.rm = TRUE), 
                      raw = correlation)
  
  output <- list(correlation = result_cor)
  if (significance_test) output$p_value <- test_result$p_value
  
  return(output)
}

#' Calculate cross-correlation function between time series
#' @param series1 First time series
#' @param series2 Second time series
#' @param type Type of CCF to return ("max_absolute", "max_positive", "min_negative", "full_sequence")
#' @param max_lag Maximum lag to consider
#' @importFrom stats ccf
#' @return List with CCF value and lag
dist_ccf <- function(series1, series2, type = "max_absolute", max_lag = NULL) {
  series1 <- validate_time_series(series1, "series1")
  series2 <- validate_time_series(series2, "series2")
  .ensure_same_length(series1, series2, "dist_ccf: ")
  
  if(length(series1) == 0) {
    warning("dist_ccf: Empty series provided.")
    if (type == "full_sequence") return(NULL)
    return(list(ccf_value = NA, lag = NA))
  }
  
  ccf_result <- stats::ccf(series1, series2, lag.max = max_lag, plot = FALSE, na.action = na.pass)
  
  if (type == "full_sequence") return(ccf_result)
  
  values <- ccf_result$acf[,,1]
  lags <- ccf_result$lag[,,1]
  
  selected_val <- NA
  selected_lag <- NA
  
  if (length(values) > 0 && !all(is.na(values))) {
    if (type == "max_absolute") {
      idx <- which.max(abs(values))
      selected_val <- values[idx]
      selected_lag <- lags[idx]
    } else if (type == "max_positive") {
      pos_vals <- values[values > 0 & !is.na(values)]
      if (length(pos_vals) > 0) {
        idx_in_pos <- which.max(pos_vals)
        orig_idx <- which(values == pos_vals[idx_in_pos] & values > 0 & !is.na(values))[1]
        selected_val <- values[orig_idx]
        selected_lag <- lags[orig_idx]
      } else {
        selected_val <- 0
        selected_lag <- NA
      }
    } else if (type == "min_negative") {
      neg_vals <- values[values < 0 & !is.na(values)]
      if (length(neg_vals) > 0) {
        idx_in_neg <- which.min(neg_vals)
        orig_idx <- which(values == neg_vals[idx_in_neg] & values < 0 & !is.na(values))[1]
        selected_val <- values[orig_idx]
        selected_lag <- lags[orig_idx]
      } else {
        selected_val <- 0
        selected_lag <- NA
      }
    }
  }
  
  return(list(ccf_value = selected_val, lag = selected_lag))
}

#' Calculate Dynamic Time Warping distance
#' @param series1 First time series
#' @param series2 Second time series
#' @param ... Additional arguments for DTW
#' @importFrom dtw dtw
#' @return List with distance
dist_dtw <- function(series1, series2, ...) {
  check_required_packages("dtw", "dist_dtw")
  
  series1 <- validate_time_series(series1, "series1")
  series2 <- validate_time_series(series2, "series2")
  
  if(length(series1) == 0 || length(series2) == 0) {
    warning("dist_dtw: Empty series provided. Returning NA.")
    return(list(distance = NA))
  }
  
  distance_val <- tryCatch({
    dtw::dtw(series1, series2, ...)$distance
  }, error = function(e) {
    warning(paste("dist_dtw: Error in DTW calculation:", e$message))
    return(NA)
  })
  
  return(list(distance = distance_val))
}

#' Calculate Normalized Mutual Information
#' @param series1 First time series
#' @param series2 Second time series
#' @param n_bins Number of bins for discretization
#' @param discretize_method Discretization method
#' @param already_discrete Whether series are already discrete
#' @importFrom aricode NMI
#' @return List with NMI value
dist_nmi <- function(series1, series2, n_bins = 10, discretize_method = "quantile", already_discrete = FALSE) {
  check_required_packages(c("aricode"), "dist_nmi")
  
  series1 <- validate_time_series(series1, "series1")
  series2 <- validate_time_series(series2, "series2")
  .ensure_same_length(series1, series2, "dist_nmi: ")
  
  if(length(series1) == 0) {
    warning("dist_nmi: Empty series. Returning NA.")
    return(list(NMI = NA))
  }
  
  ts1_discrete <- series1
  ts2_discrete <- series2
  
  if (!already_discrete) {
    is_continuous1 <- length(unique(series1)) > n_bins * 1.5 && !all(floor(series1) == series1, na.rm = TRUE)
    is_continuous2 <- length(unique(series2)) > n_bins * 1.5 && !all(floor(series2) == series2, na.rm = TRUE)
    
    if (is_continuous1 || is_continuous2) {
      discretized <- .discretize_series_pair(series1, series2, n_bins = n_bins, method = discretize_method)
      ts1_discrete <- discretized$d1
      ts2_discrete <- discretized$d2
    } else {
      ts1_discrete <- as.integer(series1)
      ts2_discrete <- as.integer(series2)
    }
  } else {
    ts1_discrete <- as.integer(series1)
    ts2_discrete <- as.integer(series2)
  }
  
  if (any(is.na(ts1_discrete)) || any(is.na(ts2_discrete))) {
    warning("dist_nmi: NAs present in discrete series. Result may be NA/NaN.")
  }
  
  nmi_value <- tryCatch({
    aricode::NMI(ts1_discrete, ts2_discrete)
  }, warning = function(w) NA, error = function(e) NA)
  
  return(list(NMI = nmi_value))
}

#' Calculate Variation of Information
#' @param series1 First time series
#' @param series2 Second time series
#' @param n_bins Number of bins for discretization
#' @param discretize_method Discretization method
#' @param already_discrete Whether series are already discrete
#' @importFrom aricode NVI
#' @return List with VI value
dist_voi <- function(series1, series2, n_bins = 10, discretize_method = "quantile", already_discrete = FALSE) {
  check_required_packages(c("aricode"), "dist_voi")
  
  series1 <- validate_time_series(series1, "series1")
  series2 <- validate_time_series(series2, "series2")
  .ensure_same_length(series1, series2, "dist_voi: ")
  
  if(length(series1) == 0) {
    warning("dist_voi: Empty series. Returning NA.")
    return(list(VI = NA))
  }
  
  ts1_discrete <- series1
  ts2_discrete <- series2
  
  if (!already_discrete) {
    is_continuous1 <- length(unique(series1)) > n_bins * 1.5 && !all(floor(series1) == series1, na.rm = TRUE)
    is_continuous2 <- length(unique(series2)) > n_bins * 1.5 && !all(floor(series2) == series2, na.rm = TRUE)
    
    if (is_continuous1 || is_continuous2) {
      discretized <- .discretize_series_pair(series1, series2, n_bins = n_bins, method = discretize_method)
      ts1_discrete <- discretized$d1
      ts2_discrete <- discretized$d2
    } else {
      ts1_discrete <- as.integer(series1)
      ts2_discrete <- as.integer(series2)
    }
  } else {
    ts1_discrete <- as.integer(series1)
    ts2_discrete <- as.integer(series2)
  }
  
  if (any(is.na(ts1_discrete)) || any(is.na(ts2_discrete))) {
    warning("dist_voi: NAs present in discrete series. Result may be NA/NaN.")
  }
  
  vi_value <- tryCatch({
          aricode::NVI(ts1_discrete, ts2_discrete)
  }, warning = function(w) NA, error = function(e) NA)
  
  return(list(VI = vi_value))
}

#' Calculate Maximal Information Coefficient
#' @param series1 First time series
#' @param series2 Second time series
#' @param ... Additional arguments for MIC
#' @importFrom minerva mine
#' @return List with MIC value
dist_mic <- function(series1, series2, ...) {
  check_required_packages("minerva", "dist_mic")
  
  series1 <- validate_time_series(series1, "series1")
  series2 <- validate_time_series(series2, "series2")
  .ensure_same_length(series1, series2, "dist_mic: ")
  
  if(length(series1) < 4 || anyNA(series1) || anyNA(series2)) {
    warning("dist_mic: Need at least 4 non-NA observations. Returning NA.")
    return(list(MIC = NA))
  }
  
  mic_value <- NA
  tryCatch({
    mine_result <- minerva::mine(series1, series2, ...)
    mic_value <- mine_result$MIC
  }, error = function(e) {
    warning(paste("dist_mic: Error in minerva::mine:", e$message))
  })
  
  return(list(MIC = mic_value))
}

#' Calculate Event Synchronization for spike trains
#' @param series1 First spike train (event times)
#' @param series2 Second spike train (event times)
#' @param tau Time tolerance for synchronization
#' @param significance_test Whether to perform significance test
#' @param n_surrogates Number of surrogates for testing
#' @param jitter_range Range for jittering in significance test
#' @param series_length_for_jitter Maximum length for boundary checking
#' @importFrom stats runif
#' @return List with synchronization measure and p-value
dist_es <- function(series1, series2, tau, significance_test = TRUE, n_surrogates = 100,
                    jitter_range = NULL, series_length_for_jitter = NULL) {
  
  series1 <- validate_time_series(series1, "series1")
  series2 <- validate_time_series(series2, "series2")
  
  # Count synchronized spikes
  count_sync_spikes <- function(train1, train2, tolerance) {
    count <- 0
    if (length(train1) == 0 || length(train2) == 0) return(0)
    
    idx_train2 <- 1
    for (t1 in train1) {
      while (idx_train2 <= length(train2) && train2[idx_train2] < t1 - tolerance) {
        idx_train2 <- idx_train2 + 1
      }
      if (idx_train2 <= length(train2) && train2[idx_train2] <= t1 + tolerance) {
        count <- count + 1
      }
    }
    return(count)
  }
  
  if (length(series1) == 0 || length(series2) == 0) {
    return(list(C = 0, p_value = if(significance_test) 1.0 else NA))
  }
  
  series1_sorted <- sort(series1)
  series2_sorted <- sort(series2)
  
  c12 <- count_sync_spikes(series1_sorted, series2_sorted, tau)
  c21 <- count_sync_spikes(series2_sorted, series1_sorted, tau)
  
  N1 <- length(series1_sorted)
  N2 <- length(series2_sorted)
  
  C12_norm <- if (N1 > 0) c12 / N1 else 0
  C21_norm <- if (N2 > 0) c21 / N2 else 0
  C_observed <- (C12_norm + C21_norm) / 2
  
  if(is.nan(C_observed)) C_observed <- 0
  
  output <- list(C = C_observed)
  
  if (significance_test) {
    jitter_val <- jitter_range
    if(is.null(jitter_val)) {
      warning("dist_es: jitter_range not set for significance test, using 5*tau.")
      jitter_val <- 5 * tau
    }
    
    C_surrogate_values <- numeric(n_surrogates)
    for (k in 1:n_surrogates) {
      s1_surrogate <- .jitter_spike_train(series1_sorted, jitter_val = jitter_val, 
                                         series_length_for_jitter = series_length_for_jitter)
      c12_surr <- count_sync_spikes(s1_surrogate, series2_sorted, tau)
      c21_surr <- count_sync_spikes(series2_sorted, s1_surrogate, tau)
      N1_surr <- length(s1_surrogate)
      C12_surr_norm <- if(N1_surr > 0) c12_surr / N1_surr else 0
      C21_surr_norm <- if(N2 > 0) c21_surr / N2 else 0
      C_surrogate_values[k] <- (C12_surr_norm + C21_surr_norm) / 2
      if(is.nan(C_surrogate_values[k])) C_surrogate_values[k] <- 0
    }
    
    p_value <- sum(C_surrogate_values >= C_observed, na.rm = TRUE) / sum(!is.na(C_surrogate_values))
    output$p_value <- p_value
  }
  
  return(output)
}

#' Calculate van Rossum distance for spike trains
#' @param series1 First spike train (event times)
#' @param series2 Second spike train (event times)
#' @param tau_vr Time constant for exponential kernel
#' @param dt Time step for discretization
#' @param t_max Maximum time for analysis
#' @param significance_test Whether to perform significance test
#' @param n_surrogates Number of surrogates for testing
#' @param jitter_range Range for jittering in significance test
#' @param series_length_for_jitter Maximum length for boundary checking
#' @return List with van Rossum distance and p-value
#' @importFrom stats sd
dist_vr <- function(series1, series2, tau_vr, dt = NULL, t_max = NULL,
                    significance_test = TRUE, n_surrogates = 100,
                    jitter_range = NULL, series_length_for_jitter = NULL) {
  
  series1 <- validate_time_series(series1, "series1")
  series2 <- validate_time_series(series2, "series2")
  
  if (tau_vr <= 0) stop("dist_vr: tau_vr must be positive.")
  
  current_dt <- dt
  if(is.null(current_dt)) current_dt <- tau_vr / 10
  if (current_dt <= 0) stop("dist_vr: dt must be positive.")
  
  # Create filtered spike trains
  create_filtered_train <- function(spike_times, time_vec, tau_constant) {
    filtered_signal <- numeric(length(time_vec))
    if (length(spike_times) == 0 || length(time_vec) == 0) return(filtered_signal)
    
    for (spike_time in spike_times) {
      if(is.na(spike_time)) next
      effective_indices <- which(time_vec >= spike_time)
      if (length(effective_indices) > 0) {
        filtered_signal[effective_indices] <- filtered_signal[effective_indices] + 
          exp(-(time_vec[effective_indices] - spike_time) / tau_constant)
      }
    }
    return(filtered_signal)
  }
  
  calc_t_max <- t_max
  if(is.null(calc_t_max)){
    max_s1 <- if(length(series1) > 0) max(series1, na.rm = TRUE) else 0
    max_s2 <- if(length(series2) > 0) max(series2, na.rm = TRUE) else 0
    calc_t_max <- max(max_s1, max_s2, 0, na.rm = TRUE) + 5 * tau_vr
    if (calc_t_max == 0 && (length(series1) > 0 || length(series2) > 0)) calc_t_max <- 5 * tau_vr
    if (calc_t_max == 0) calc_t_max <- current_dt
  }
  if(calc_t_max <= 0 && !(length(series1) == 0 && length(series2) == 0)) calc_t_max <- max(5 * tau_vr, current_dt)
  else if (calc_t_max <= 0) calc_t_max <- current_dt
  
  min_time_in_pair <- 0
  if(length(series1) > 0 || length(series2) > 0) {
    min_time_in_pair <- min(series1, series2, Inf, na.rm = TRUE)
    if(is.infinite(min_time_in_pair)) min_time_in_pair <- 0
  }
  
  time_vector <- seq(min_time_in_pair, calc_t_max, by = current_dt)
  if(length(time_vector) == 0 && calc_t_max >= min_time_in_pair && current_dt > 0) {
    time_vector <- c(min_time_in_pair, calc_t_max)
  }
  if(length(time_vector) == 0) time_vector <- c(min_time_in_pair)
  
  f1_filtered <- create_filtered_train(series1, time_vector, tau_vr)
  f2_filtered <- create_filtered_train(series2, time_vector, tau_vr)
  
  distance_observed <- NA
  if(length(f1_filtered) == 0 && length(f2_filtered) == 0) {
    distance_observed <- 0
  } else if (length(f1_filtered) != length(f2_filtered)) {
    warning("dist_vr: Filtered trains have different lengths. Returning NA.")
    distance_observed <- NA
  } else {
    sum_squared_diff <- sum((f1_filtered - f2_filtered)^2, na.rm = TRUE)
    distance_observed <- if (tau_vr == 0) NA else sqrt(sum_squared_diff * current_dt / tau_vr)
    if(tau_vr == 0) warning("dist_vr: tau_vr is zero. Distance ill-defined.")
  }
  
  output <- list(van_rossum_distance = distance_observed)
  
  if (significance_test && !is.na(distance_observed)) {
    jitter_val <- jitter_range
    if(is.null(jitter_val)) {
      warning("dist_vr: jitter_range not set for significance test, using 5*tau_vr.")
      jitter_val <- 5 * tau_vr
    }
    
    distance_squared_surrogates <- numeric(n_surrogates)
    for (k in 1:n_surrogates) {
      s1_surrogate <- .jitter_spike_train(series1, jitter_val, series_length_for_jitter = series_length_for_jitter)
      f1_surrogate_filtered <- create_filtered_train(s1_surrogate, time_vector, tau_vr)
      if(length(f1_surrogate_filtered) != length(f2_filtered)) {
        distance_squared_surrogates[k] <- NA
        next
      }
      sum_squared_diff_surrogate <- sum((f1_surrogate_filtered - f2_filtered)^2, na.rm = TRUE)
      distance_squared_surrogates[k] <- if(tau_vr == 0) NA else sum_squared_diff_surrogate * current_dt / tau_vr
    }
    
    valid_surrogates <- sqrt(distance_squared_surrogates[!is.na(distance_squared_surrogates)])
    output$p_value <- if(length(valid_surrogates) > 0) {
      sum(valid_surrogates <= distance_observed, na.rm = TRUE) / length(valid_surrogates)
    } else {
      NA
    }
  } else if (significance_test) {
    output$p_value <- NA
  }
  
  return(output)
}

# Native Distance Functions
# Implementation of common distance measures without external dependencies

#' Calculate Euclidean distance between two vectors
#' @param x First vector
#' @param y Second vector
#' @return Euclidean distance
euclidean_distance <- function(x, y) {
  if (length(x) != length(y)) {
    stop("Vectors must have the same length")
  }
  sqrt(sum((x - y)^2))
}

#' Calculate Manhattan (L1) distance between two vectors
#' @param x First vector
#' @param y Second vector
#' @return Manhattan distance
manhattan_distance <- function(x, y) {
  if (length(x) != length(y)) {
    stop("Vectors must have the same length")
  }
  sum(abs(x - y))
}

#' Calculate Cosine distance between two vectors
#' @param x First vector
#' @param y Second vector
#' @return Cosine distance (1 - cosine similarity)
cosine_distance <- function(x, y) {
  if (length(x) != length(y)) {
    stop("Vectors must have the same length")
  }
  
  dot_product <- sum(x * y)
  norm_x <- sqrt(sum(x^2))
  norm_y <- sqrt(sum(y^2))
  
  if (norm_x == 0 || norm_y == 0) {
    return(1)  # Maximum distance for zero vectors
  }
  
  cosine_sim <- dot_product / (norm_x * norm_y)
  return(1 - cosine_sim)
}

#' Calculate Minkowski distance between two vectors
#' @param x First vector
#' @param y Second vector
#' @param p Power parameter (p=1: Manhattan, p=2: Euclidean)
#' @return Minkowski distance
minkowski_distance <- function(x, y, p = 2) {
  if (length(x) != length(y)) {
    stop("Vectors must have the same length")
  }
  (sum(abs(x - y)^p))^(1/p)
}

#' Calculate Chebyshev (maximum) distance between two vectors
#' @param x First vector
#' @param y Second vector
#' @return Chebyshev distance
chebyshev_distance <- function(x, y) {
  if (length(x) != length(y)) {
    stop("Vectors must have the same length")
  }
  max(abs(x - y))
}

#' Calculate Canberra distance between two vectors
#' @param x First vector
#' @param y Second vector
#' @return Canberra distance
canberra_distance <- function(x, y) {
  if (length(x) != length(y)) {
    stop("Vectors must have the same length")
  }
  
  numerator <- abs(x - y)
  denominator <- abs(x) + abs(y)
  
  # Handle division by zero
  valid_indices <- denominator != 0
  if (sum(valid_indices) == 0) {
    return(0)
  }
  
  sum(numerator[valid_indices] / denominator[valid_indices])
}

#' Calculate correlation-based distance
#' @param x First vector
#' @param y Second vector
#' @param method Correlation method ("pearson", "spearman", "kendall")
#' @return Correlation distance (1 - |correlation|)
#' @importFrom stats cor
correlation_distance <- function(x, y, method = "pearson") {
  if (length(x) != length(y)) {
    stop("Vectors must have the same length")
  }
  
  if (length(x) < 3) {
    warning("Correlation requires at least 3 observations")
    return(NA)
  }
  
  correlation <- cor(x, y, method = method, use = "complete.obs")
  if (is.na(correlation)) {
    return(1)  # Maximum distance for undefined correlation
  }
  
  return(1 - abs(correlation))
}

#' Native DTW implementation (simplified)
#' @param x First vector
#' @param y Second vector
#' @return DTW distance
dtw_distance <- function(x, y) {
  n <- length(x)
  m <- length(y)
  
  # Create distance matrix
  dist_matrix <- matrix(0, n, m)
  for (i in 1:n) {
    for (j in 1:m) {
      dist_matrix[i, j] <- abs(x[i] - y[j])
    }
  }
  
  # Initialize DTW matrix
  dtw_matrix <- matrix(Inf, n + 1, m + 1)
  dtw_matrix[1, 1] <- 0
  
  # Fill DTW matrix
  for (i in 2:(n + 1)) {
    for (j in 2:(m + 1)) {
      cost <- dist_matrix[i - 1, j - 1]
      dtw_matrix[i, j] <- cost + min(
        dtw_matrix[i - 1, j],     # insertion
        dtw_matrix[i, j - 1],     # deletion
        dtw_matrix[i - 1, j - 1]  # match
      )
    }
  }
  
  return(dtw_matrix[n + 1, m + 1])
}

#' Calculate all available native distances
#' @param x First vector
#' @param y Second vector
#' @param method Distance method name
#' @param ... Additional parameters
#' @return Distance value
#' @importFrom stats median
calculate_native_distance <- function(x, y, method, ...) {
  method <- tolower(method)
  
  result <- switch(method,
    "euclidean" = euclidean_distance(x, y),
    "manhattan" = manhattan_distance(x, y),
    "cosine" = cosine_distance(x, y),
    "minkowski" = minkowski_distance(x, y, ...),
    "chebyshev" = chebyshev_distance(x, y),
    "canberra" = canberra_distance(x, y),
    "correlation" = correlation_distance(x, y, ...),
    "dtw" = dtw_distance(x, y),
    "pearson" = correlation_distance(x, y, method = "pearson"),
    "spearman" = correlation_distance(x, y, method = "spearman"),
    stop(paste("Unknown distance method:", method))
  )
  
  return(result)
}

#' Get list of all available native distance methods
#' @return Vector of method names
get_native_methods <- function() {
  c("euclidean", "manhattan", "cosine", "minkowski", "chebyshev", 
    "canberra", "correlation", "dtw", "pearson", "spearman")
}

#' Check if method is a native implementation
#'
#' Checks if a given distance method name corresponds to a native implementation
#' within the package.
#'
#' @param method Character string. The name of the method to check.
#' @return Logical. `TRUE` if the method is native, `FALSE` otherwise.
#' @keywords internal
is_native_method <- function(method) {
  tolower(method) %in% get_native_methods()
}
