# Helper functions for detect_regime function

#' @importFrom graphics hist
#' @importFrom zoo rollapply na.locf
#' @importFrom stats sd lm coef var quantile pnorm na.pass
#' @importFrom utils globalVariables

utils::globalVariables(c("complexity_vector", "params", "min_change_size"))

.process_regime_input_v2 <- function(complexity_data, complexity_col) {
  if (is.numeric(complexity_data) && is.vector(complexity_data)) {
    if (length(complexity_data) == 0) stop("Input complexity vector is empty.", call. = FALSE)
    return(list(
      complexity = complexity_data,
      time_index = 1:length(complexity_data)
    ))
  }
  
  if (!is.data.frame(complexity_data)) {
    stop("complexity_data must be a data frame or numeric vector.", call. = FALSE)
  }
  if (nrow(complexity_data) == 0) stop("Input complexity data frame is empty.", call. = FALSE)
  
  complexity_col_internal <- complexity_col
  if (is.null(complexity_col_internal)) {
    complexity_col_internal <- .detect_complexity_column_v2(complexity_data) 
  } else {
    if (!is.character(complexity_col_internal) || length(complexity_col_internal) != 1) {
      stop("complexity_col must be a single character string.", call. = FALSE)
    }
    if (!complexity_col_internal %in% names(complexity_data)) {
      stop(paste("Specified complexity_col '", complexity_col_internal, "' not found in data frame."), call. = FALSE)
    }
    if (!is.numeric(complexity_data[[complexity_col_internal]])) {
      stop(paste("Column '", complexity_col_internal, "' must be numeric."), call. = FALSE)
    }
  }
  
  extracted_complexity <- complexity_data[[complexity_col_internal]]
  
  time_idx <- NULL
  if ("time_index" %in% names(complexity_data)) {
    time_idx_candidate <- complexity_data$time_index
    if(length(time_idx_candidate) != length(extracted_complexity)){
      warning("'time_index' length differs from complexity column length. Using 1:n for analysis time index.", call. = FALSE)
      time_idx <- 1:length(extracted_complexity)
    } else {
      time_idx <- time_idx_candidate
    }
  } else {
    time_idx <- 1:length(extracted_complexity)
  }
  
  return(list(
    complexity = extracted_complexity,
    time_index = time_idx
  ))
}

.detect_complexity_column_v2 <- function(data_frame) {
  potential_cols <- c("complexity", "complexity_score", "complexity_value", 
                     "complexity_score_value", "score", "value")
  
  for (col in potential_cols) {
    if (col %in% names(data_frame) && is.numeric(data_frame[[col]])) {
      return(col)
    }
  }
  
  numeric_cols <- names(data_frame)[sapply(data_frame, is.numeric)]
  if (length(numeric_cols) > 0) {
    return(numeric_cols[1])
  }
  
  stop("No suitable complexity column found in data frame.", call. = FALSE)
}

.get_sensitivity_params_v2 <- function(sensitivity_level, base_window_size,
                                       base_peak_threshold, base_cumulative_threshold) {
  params <- list(
    window_size = max(3, floor(base_window_size)), 
    peak_threshold = base_peak_threshold,
    cumulative_threshold = base_cumulative_threshold,
    change_test_threshold = 0.05, 
    variance_ratio_threshold = 2.0, 
    gradient_signif_threshold = 0.1, 
    entropy_nbins = 10,
    entropy_relative_change_threshold = 0.2,
    sensitivity = sensitivity_level 
  )
  
  if (sensitivity_level == "high") {
    params$window_size <- max(3, floor(params$window_size * 0.7)) 
    params$peak_threshold <- params$peak_threshold * 0.6 
    params$cumulative_threshold <- params$cumulative_threshold * 0.5 
    params$change_test_threshold <- params$change_test_threshold * 0.6
    params$variance_ratio_threshold <- 1 + (params$variance_ratio_threshold - 1) * 0.6 
    params$gradient_signif_threshold <- params$gradient_signif_threshold * 0.7
    params$entropy_nbins <- max(5, floor(params$entropy_nbins * 0.7))
    params$entropy_relative_change_threshold <- params$entropy_relative_change_threshold * 0.6
  } else if (sensitivity_level == "low") {
    params$window_size <- max(5, floor(params$window_size * 1.3)) 
    params$peak_threshold <- params$peak_threshold * 1.4 
    params$cumulative_threshold <- min(1, params$cumulative_threshold * 1.3) 
    params$change_test_threshold <- params$change_test_threshold * 1.4
    params$variance_ratio_threshold <- 1 + (params$variance_ratio_threshold - 1) * 1.4
    params$gradient_signif_threshold <- params$gradient_signif_threshold * 1.3
    params$entropy_nbins <- floor(params$entropy_nbins * 1.3)
    params$entropy_relative_change_threshold <- params$entropy_relative_change_threshold * 1.4
  }
  
  if(params$window_size %% 2 == 0) params$window_size <- params$window_size + 1
  params$entropy_nbins <- max(2, params$entropy_nbins)
  
  return(params)
}

.apply_min_change_constraint <- function(regime_changes, min_change_size) {
  if (length(regime_changes) == 0 || min_change_size <= 1) {
    return(regime_changes)
  }
  
  n <- length(regime_changes)
  constrained_changes <- regime_changes
  change_indices <- which(regime_changes)
  
  if (length(change_indices) <= 1) {
    return(constrained_changes)
  }
  
  for (i in 2:length(change_indices)) {
    distance <- change_indices[i] - change_indices[i-1]
    if (distance < min_change_size) {
      constrained_changes[change_indices[i]] <- FALSE
    }
  }
  
  return(constrained_changes)
}

.generate_regime_ids <- function(regime_changes) {
  if (length(regime_changes) == 0) {
    return(integer(0))
  }
  
  regime_ids <- integer(length(regime_changes))
  current_regime <- 1L
  regime_ids[1] <- current_regime
  
  for (i in 2:length(regime_changes)) {
    if (regime_changes[i]) {
      current_regime <- current_regime + 1L
    }
    regime_ids[i] <- current_regime
  }
  
  return(regime_ids)
}

.apply_single_method_v2 <- function(complexity_vector, time_index, method, params, min_change_size) {
  switch(method,
         "cumulative_peaks" = .detect_cumulative_peaks_v2(complexity_vector, time_index, params, min_change_size),
         "changepoint" = .detect_changepoints_v2(complexity_vector, time_index, params, min_change_size), 
         "threshold" = .detect_threshold_regimes_v2(complexity_vector, time_index, params, min_change_size),
         "variance_shift" = .detect_variance_shifts_v2(complexity_vector, time_index, params, min_change_size),
         "gradient" = .detect_gradient_changes_v2(complexity_vector, time_index, params, min_change_size),
         "entropy" = .detect_entropy_changes_v2(complexity_vector, time_index, params, min_change_size),
         stop("Unknown single method in .apply_single_method_v2: ", method, call. = FALSE)
  )
}

.apply_all_methods_v2 <- function(complexity_vector, time_index, params, min_change_size, single_methods_list) {
  all_results_list <- list()
  n <- length(complexity_vector)
  
  for (s_method in single_methods_list) {
    res <- tryCatch({
      .apply_single_method_v2(complexity_vector, time_index, s_method, params, min_change_size)
    }, error = function(e) {
      warning("Method '", s_method, "' failed in 'all' mode: ", e$message, call. = FALSE)
      list(regime_change = rep(FALSE, n), regime_id = rep(1L, n),
           change_type = rep("error", n), change_magnitude = rep(0, n), confidence = rep(0, n))
    })
    all_results_list[[s_method]] <- res$regime_change 
  }
  
  if (length(all_results_list) == 0) {
    warning("No methods successfully ran in 'all' mode.", call. = FALSE)
    return(list(regime_change = rep(FALSE, n), regime_id = rep(1L, n), change_type = rep("all_failed", n),
                change_magnitude = rep(0,n), confidence = rep(0,n)))
  }
  
  vote_matrix <- do.call(cbind, all_results_list)
  if(is.null(dim(vote_matrix)) && length(vote_matrix) == n ) vote_matrix <- matrix(vote_matrix, nrow = n) 
  else if (is.null(dim(vote_matrix)) && length(vote_matrix) != n ) {
    warning("Vote matrix malformed in 'all' mode. Returning no changes.", call. = FALSE)
    return(list(regime_change = rep(FALSE, n), regime_id = rep(1L, n), change_type = rep("vote_matrix_error", n),
                change_magnitude = rep(0,n), confidence = rep(0,n)))
  }
  
  ensemble_vote_threshold_prop <- 0.3 
  if (!is.null(params$sensitivity)) {
    if(params$sensitivity == "high") ensemble_vote_threshold_prop <- 0.15 
    if(params$sensitivity == "low") ensemble_vote_threshold_prop <- 0.5 
  }
  
  num_methods_voted = ncol(vote_matrix)
  if(is.null(num_methods_voted) || num_methods_voted == 0) { 
    vote_sums <- rep(0, n)
    consensus_changes <- rep(FALSE, n)
    num_methods_voted <- 0 
  } else {
    vote_sums <- rowSums(vote_matrix, na.rm = TRUE)
    min_votes_needed <- max(1, floor(ensemble_vote_threshold_prop * num_methods_voted))
    consensus_changes <- vote_sums >= min_votes_needed
  }
  
  consensus_changes <- .apply_min_change_constraint(consensus_changes, min_change_size)
  final_regime_id <- .generate_regime_ids(consensus_changes)
  
  final_change_type <- rep("none", n)
  final_change_magnitude <- vote_sums / max(1, num_methods_voted) 
  final_confidence <- vote_sums / max(1, num_methods_voted)
  
  final_change_type[consensus_changes] <- paste0("ensemble_vote_", round(final_change_magnitude[consensus_changes]*100), "%")
  
  return(list(
    regime_change = consensus_changes,
    regime_id = final_regime_id,
    change_type = final_change_type,
    change_magnitude = final_change_magnitude,
    confidence = final_confidence
  ))
}

.detect_cumulative_peaks_v2 <- function(complexity_vector, time_index, params, min_change_size) {
  n <- length(complexity_vector)
  window <- params$window_size
  z_score_thresh <- params$peak_threshold
  cumulative_prop_thresh <- params$cumulative_threshold
  
  if (n < window) { 
    return(list(regime_change = rep(FALSE, n), regime_id = rep(1L, n), change_type = rep("peaks_insufficient_data", n),
                change_magnitude = rep(0,n), confidence = rep(0,n)))
  }
  
  # Simple peak detection based on local maxima
  final_regime_changes <- rep(FALSE, n)
  if (n > 2) {
    for (i in 2:(n-1)) {
      if (complexity_vector[i] > complexity_vector[i-1] && complexity_vector[i] > complexity_vector[i+1]) {
        final_regime_changes[i] <- TRUE
      }
    }
  }
  
  final_regime_changes <- .apply_min_change_constraint(final_regime_changes, min_change_size)
  final_regime_id <- .generate_regime_ids(final_regime_changes)
  
  final_change_type <- rep("none", n)
  final_change_magnitude <- rep(0, n)
  final_confidence <- rep(0.5, n)
  
  final_change_type[final_regime_changes] <- "peak_detected"
  final_change_magnitude[final_regime_changes] <- abs(diff(c(0, complexity_vector)))[final_regime_changes]
  
  return(list(
    regime_change = final_regime_changes,
    regime_id = final_regime_id,
    change_type = final_change_type,
    change_magnitude = final_change_magnitude,
    confidence = final_confidence
  ))
}

.detect_changepoints_v2 <- function(complexity_vector, time_index, params, min_change_size) {
  n <- length(complexity_vector)
  window <- params$window_size
  change_thresh <- params$change_test_threshold
  
  if (n < window * 2) {
    return(list(regime_change = rep(FALSE, n), regime_id = rep(1L, n), change_type = rep("cp_insufficient_data", n),
                change_magnitude = rep(0,n), confidence = rep(0,n)))
  }
  
  final_regime_changes <- rep(FALSE, n)
  change_magnitudes <- rep(0, n)
  
  for (i in (window + 1):(n - window)) {
    left_window <- complexity_vector[(i - window):(i - 1)]
    right_window <- complexity_vector[(i + 1):(i + window)]
    
    left_window <- left_window[!is.na(left_window)]
    right_window <- right_window[!is.na(right_window)]
    
    if (length(left_window) >= 3 && length(right_window) >= 3) {
      mean_left <- mean(left_window)
      mean_right <- mean(right_window)
      
      pooled_sd <- sqrt((var(left_window) + var(right_window)) / 2)
      if (pooled_sd > 1e-9) {
        t_stat <- abs(mean_left - mean_right) / (pooled_sd * sqrt(2/window))
        p_value <- 2 * (1 - pnorm(abs(t_stat)))
        
        if (p_value < change_thresh) {
          final_regime_changes[i] <- TRUE
          change_magnitudes[i] <- abs(mean_left - mean_right)
        }
      }
    }
  }
  
  final_regime_changes <- .apply_min_change_constraint(final_regime_changes, min_change_size)
  final_regime_id <- .generate_regime_ids(final_regime_changes)
  
  final_change_type <- rep("none", n)
  final_confidence <- rep(0, n)
  
  final_change_type[final_regime_changes] <- "changepoint"
  final_confidence[final_regime_changes] <- 0.7
  
  return(list(
    regime_change = final_regime_changes,
    regime_id = final_regime_id,
    change_type = final_change_type,
    change_magnitude = change_magnitudes,
    confidence = final_confidence
  ))
}

.detect_threshold_regimes_v2 <- function(complexity_vector, time_index, params, min_change_size) {
  n <- length(complexity_vector)
  
  quartiles <- quantile(complexity_vector, probs = c(0.25, 0.75), na.rm = TRUE)
  
  regime_labels <- rep(2L, n)
  regime_labels[complexity_vector <= quartiles[1]] <- 1L
  regime_labels[complexity_vector >= quartiles[2]] <- 3L
  
  final_regime_changes <- rep(FALSE, n)
  if (n > 1) {
    final_regime_changes[2:n] <- diff(regime_labels) != 0
  }
  
  final_regime_changes <- .apply_min_change_constraint(final_regime_changes, min_change_size)
  final_regime_id <- .generate_regime_ids(final_regime_changes)
  
  final_change_type <- rep("none", n)
  final_change_magnitude <- rep(0, n)
  final_confidence <- rep(0.7, n)
  
  final_change_type[final_regime_changes] <- "threshold_change"
  final_change_magnitude[final_regime_changes] <- abs(diff(c(0, complexity_vector)))[final_regime_changes]
  
  return(list(
    regime_change = final_regime_changes,
    regime_id = final_regime_id,
    change_type = final_change_type,
    change_magnitude = final_change_magnitude,
    confidence = final_confidence
  ))
}

.detect_variance_shifts_v2 <- function(complexity_vector, time_index, params, min_change_size) {
  n <- length(complexity_vector)
  window <- params$window_size
  var_ratio_thresh <- params$variance_ratio_threshold
  
  if (n < window * 2) {
    return(list(regime_change = rep(FALSE, n), regime_id = rep(1L, n), change_type = rep("var_insufficient_data", n),
                change_magnitude = rep(0,n), confidence = rep(0,n)))
  }
  
  final_regime_changes <- rep(FALSE, n)
  variance_ratios <- rep(1, n)
  
  for (i in (window + 1):(n - window)) {
    left_window <- complexity_vector[(i - window):(i - 1)]
    right_window <- complexity_vector[(i + 1):(i + window)]
    
    left_window <- left_window[!is.na(left_window)]
    right_window <- right_window[!is.na(right_window)]
    
    if (length(left_window) >= 3 && length(right_window) >= 3) {
      var_left <- var(left_window)
      var_right <- var(right_window)
      
      if (var_left > 1e-9 && var_right > 1e-9) {
        var_ratio <- max(var_left, var_right) / min(var_left, var_right)
        variance_ratios[i] <- var_ratio
        
        if (var_ratio > var_ratio_thresh) {
          final_regime_changes[i] <- TRUE
        }
      }
    }
  }
  
  final_regime_changes <- .apply_min_change_constraint(final_regime_changes, min_change_size)
  final_regime_id <- .generate_regime_ids(final_regime_changes)
  
  final_change_type <- rep("none", n)
  final_confidence <- rep(0, n)
  
  final_change_type[final_regime_changes] <- "variance_shift"
  final_confidence[final_regime_changes] <- 0.6
  
  return(list(
    regime_change = final_regime_changes,
    regime_id = final_regime_id,
    change_type = final_change_type,
    change_magnitude = variance_ratios - 1,
    confidence = final_confidence
  ))
}

.detect_gradient_changes_v2 <- function(complexity_vector, time_index, params, min_change_size) {
  n <- length(complexity_vector)
  window <- params$window_size
  gradient_thresh <- params$gradient_signif_threshold
  
  if (n < window * 2) {
    return(list(regime_change = rep(FALSE, n), regime_id = rep(1L, n), change_type = rep("grad_insufficient_data", n),
                change_magnitude = rep(0,n), confidence = rep(0,n)))
  }
  
  final_regime_changes <- rep(FALSE, n)
  gradient_changes <- rep(0, n)
  
  # Simple gradient detection using differences
  if (n > 2) {
    gradients <- diff(complexity_vector)
    gradient_changes_raw <- abs(diff(gradients))
    
    for (i in 2:(n-1)) {
      if (gradient_changes_raw[i-1] > gradient_thresh) {
        final_regime_changes[i] <- TRUE
        gradient_changes[i] <- gradient_changes_raw[i-1]
      }
    }
  }
  
  final_regime_changes <- .apply_min_change_constraint(final_regime_changes, min_change_size)
  final_regime_id <- .generate_regime_ids(final_regime_changes)
  
  final_change_type <- rep("none", n)
  final_confidence <- rep(0, n)
  
  final_change_type[final_regime_changes] <- "gradient_change"
  final_confidence[final_regime_changes] <- 0.6
  
  return(list(
    regime_change = final_regime_changes,
    regime_id = final_regime_id,
    change_type = final_change_type,
    change_magnitude = gradient_changes,
    confidence = final_confidence
  ))
}

.detect_entropy_changes_v2 <- function(complexity_vector, time_index, params, min_change_size) {
  n <- length(complexity_vector)
  
  # Simple entropy-like measure using local variance
  final_regime_changes <- rep(FALSE, n)
  window <- params$window_size
  
  if (n < window) {
    return(list(regime_change = rep(FALSE, n), regime_id = rep(1L, n), change_type = rep("ent_insufficient_data", n),
                change_magnitude = rep(0,n), confidence = rep(0,n)))
  }
  
  # Calculate rolling variance as proxy for entropy
  rolling_vars <- rep(0, n)
  half_window <- floor(window / 2)
  
  for (i in 1:n) {
    start_idx <- max(1, i - half_window)
    end_idx <- min(n, i + half_window)
    window_data <- complexity_vector[start_idx:end_idx]
    rolling_vars[i] <- var(window_data, na.rm = TRUE)
  }
  
  # Detect changes in variance
  if (n > 1) {
    var_changes <- abs(diff(rolling_vars))
    threshold <- quantile(var_changes, 0.9, na.rm = TRUE)
    
    for (i in 2:n) {
      if (var_changes[i-1] > threshold) {
        final_regime_changes[i] <- TRUE
      }
    }
  }
  
  final_regime_changes <- .apply_min_change_constraint(final_regime_changes, min_change_size)
  final_regime_id <- .generate_regime_ids(final_regime_changes)
  
  final_change_type <- rep("none", n)
  final_change_magnitude <- rep(0, n)
  final_confidence <- rep(0, n)
  
  final_change_type[final_regime_changes] <- "entropy_change"
  final_change_magnitude[final_regime_changes] <- abs(diff(c(0, rolling_vars)))[final_regime_changes]
  final_confidence[final_regime_changes] <- 0.5
  
  return(list(
    regime_change = final_regime_changes,
    regime_id = final_regime_id,
    change_type = final_change_type,
    change_magnitude = final_change_magnitude,
    confidence = final_confidence
  ))
}

.detect_smart_combination_v2 <- function(complexity_vector, time_index, params, min_change_size) {
  gradient_result <- .detect_gradient_changes_v2(complexity_vector, time_index, params, min_change_size)
  peaks_result <- .detect_cumulative_peaks_v2(complexity_vector, time_index, params, min_change_size)
  changepoint_result <- .detect_changepoints_v2(complexity_vector, time_index, params, min_change_size)
  
  n <- length(complexity_vector)
  
  combined_changes <- gradient_result$regime_change | peaks_result$regime_change | changepoint_result$regime_change
  combined_changes <- .apply_min_change_constraint(combined_changes, min_change_size)
  final_regime_id <- .generate_regime_ids(combined_changes)
  
  final_change_type <- rep("none", n)
  final_change_magnitude <- rep(0, n)
  final_confidence <- rep(0, n)
  
  change_points_indices <- which(combined_changes)
  if(length(change_points_indices) > 0){
    for(cp_idx in change_points_indices){
      methods_detected <- c()
      confidences <- c()
      magnitudes <- c()
      
      if(gradient_result$regime_change[cp_idx]) {
        methods_detected <- c(methods_detected, "gradient")
        confidences <- c(confidences, gradient_result$confidence[cp_idx])
        magnitudes <- c(magnitudes, gradient_result$change_magnitude[cp_idx])
      }
      if(peaks_result$regime_change[cp_idx]) {
        methods_detected <- c(methods_detected, "peaks")
        confidences <- c(confidences, peaks_result$confidence[cp_idx])
        magnitudes <- c(magnitudes, peaks_result$change_magnitude[cp_idx])
      }
      if(changepoint_result$regime_change[cp_idx]) {
        methods_detected <- c(methods_detected, "changepoint")
        confidences <- c(confidences, changepoint_result$confidence[cp_idx])
        magnitudes <- c(magnitudes, changepoint_result$change_magnitude[cp_idx])
      }
      
      final_change_type[cp_idx] <- paste("smart", paste(methods_detected, collapse = "+"), sep = "_")
      final_confidence[cp_idx] <- mean(confidences, na.rm = TRUE)
      final_change_magnitude[cp_idx] <- mean(magnitudes, na.rm = TRUE)
    }
  }
  
  return(list(
    regime_change = combined_changes,
    regime_id = final_regime_id,
    change_type = final_change_type,
    change_magnitude = final_change_magnitude,
    confidence = final_confidence
  ))
}

.consolidate_similar_regimes_v2 <- function(result_regime_info, similarity_threshold, min_change_size) {
  regime_ids <- unique(result_regime_info$regime_id)
  n_regimes <- length(regime_ids)
  
  if (n_regimes <= 1) {
    return(list(result = result_regime_info, is_consolidated_event = rep(FALSE, nrow(result_regime_info))))
  }
  
  regime_means <- sapply(regime_ids, function(rid) {
    regime_data <- result_regime_info$complexity[result_regime_info$regime_id == rid]
    mean(regime_data, na.rm = TRUE)
  })
  
  consolidated_mapping <- regime_ids
  names(consolidated_mapping) <- regime_ids
  
  for (i in 1:(n_regimes-1)) {
    for (j in (i+1):n_regimes) {
      mean_diff <- abs(regime_means[i] - regime_means[j])
      pooled_mean <- mean(c(regime_means[i], regime_means[j]))
      
      if (pooled_mean > 1e-9) {
        normalized_diff <- mean_diff / pooled_mean
        if (normalized_diff < similarity_threshold) {
          consolidated_mapping[regime_ids[j]] <- regime_ids[i]
        }
      }
    }
  }
  
  result_regime_info$regime_id <- consolidated_mapping[as.character(result_regime_info$regime_id)]
  
  unique_consolidated <- unique(result_regime_info$regime_id)
  renumber_map <- seq_along(unique_consolidated)
  names(renumber_map) <- unique_consolidated
  result_regime_info$regime_id <- renumber_map[as.character(result_regime_info$regime_id)]
  
  is_consolidated <- rep(FALSE, nrow(result_regime_info))
  
  return(list(result = result_regime_info, is_consolidated_event = is_consolidated))
}

.calculate_stability_markings_v2 <- function(result_regime_info, min_change_size) {
  n <- nrow(result_regime_info)
  stability_score <- rep(0.5, n)
  regime_stability <- rep("Stable", n)
  
  if (n > 1) {
    regime_durations <- table(result_regime_info$regime_id)
    
    for (i in 1:n) {
      current_regime <- result_regime_info$regime_id[i]
      regime_duration <- regime_durations[as.character(current_regime)]
      
      stability_score[i] <- min(1, regime_duration / (min_change_size * 2))
      
      if (stability_score[i] > 0.7) {
        regime_stability[i] <- "Stable"
      } else if (stability_score[i] > 0.3) {
        regime_stability[i] <- "Transitional"
      } else {
        regime_stability[i] <- "Unstable"
      }
    }
  }
  
  return(list(
    stability_score = stability_score,
    regime_stability = regime_stability
  ))
}

.describe_complexity_patterns_v2 <- function(result_regime_info) {
  n <- nrow(result_regime_info)
  complexity_pattern <- rep("Stable", n)
  
  if (n > 1) {
    regime_ids <- unique(result_regime_info$regime_id)
    
    for (rid in regime_ids) {
      regime_indices <- which(result_regime_info$regime_id == rid)
      
      if (length(regime_indices) > 2) {
        regime_complexity <- result_regime_info$complexity[regime_indices]
        
        trend_lm <- lm(regime_complexity ~ seq_along(regime_complexity))
        slope <- coef(trend_lm)[2]
        
        if (!is.na(slope)) {
          if (slope > 0.01) {
            complexity_pattern[regime_indices] <- "Increasing"
          } else if (slope < -0.01) {
            complexity_pattern[regime_indices] <- "Decreasing"
          } else {
            complexity_var <- var(regime_complexity, na.rm = TRUE)
            complexity_mean <- mean(regime_complexity, na.rm = TRUE)
            
            if (complexity_mean > 1e-9 && complexity_var / complexity_mean > 0.1) {
              complexity_pattern[regime_indices] <- "Fluctuating"
            } else {
              complexity_pattern[regime_indices] <- "Stable"
            }
          }
        }
      }
    }
  }
  
  return(complexity_pattern)
}
