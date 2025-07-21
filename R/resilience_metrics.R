#' Calculate Comprehensive Resilience Metrics for TSnetworks
#'
#' @description
#' Calculates comprehensive resilience metrics across three capacity dimensions
#' (absorptive, restorative, adaptive) that are fully compatible with existing
#' TSnetworks data structures and workflows. This function adds resilience
#' analysis capabilities without modifying any existing TSnetworks functions.
#'
#' @param data Data frame containing time series data (compatible with rolling_measures(), stna(), etc.)
#' @param ts_cols Character vector of time series column names. If NULL, auto-detects numeric columns
#' @param id_col Character string for grouping column (same format as rolling_measures())
#' @param window_width Integer specifying rolling window size (default: 7)
#' @param scaling_method Character string specifying scaling method (uses existing apply_scaling())
#' @param capacity_types Character vector specifying which capacities to calculate
#' @param baseline_method Method for establishing baseline performance
#' @param na_action How to handle missing values (same as rolling_measures())
#'
#' @return Data frame with original data plus added resilience metric columns.
#'   New columns follow naming pattern: "{capacity}_{metric}" (e.g., "absorptive_vsi")
#'
#' @examples
#' # Basic resilience analysis
#' data(saqrsteps)
#' resilience_data <- calculate_resilience_metrics(saqrsteps, ts_cols = "Steps")
#' head(resilience_data)
#'
#' # Integrate with existing TSnetworks workflow
#' stna_result <- stna(saqrsteps, "Steps", num_states = 4, method = "quantile")
#' enhanced_data <- calculate_resilience_metrics(stna_result$data, ts_cols = "Steps")
#' 
#' # Use with existing plotting functions
#' plot_timeseries_enhanced(enhanced_data, ts_col = "Steps")
#'
#' @importFrom stats var sd median mad IQR acf lm coef
#' @importFrom zoo rollapply
#' @export
calculate_resilience_metrics <- function(data,
                                        ts_cols = NULL,
                                        id_col = NULL,
                                        window_width = 7L,
                                        scaling_method = "auto",
                                        capacity_types = c("absorptive", "restorative", "adaptive"),
                                        baseline_method = "auto",
                                        na_action = c("pass", "omit", "fail")) {
  
  # Input validation using TSnetworks patterns
  .validate_resilience_inputs(data, ts_cols, id_col, window_width, capacity_types)
  
  na_action <- match.arg(na_action)
  
  # Auto-detect time series columns using same logic as rolling_measures
  if (is.null(ts_cols)) {
    exclude_cols <- if (!is.null(id_col)) id_col else character(0)
    candidate_cols <- setdiff(names(data), exclude_cols)
    ts_cols <- candidate_cols[sapply(data[candidate_cols], is.numeric)]
    
    if (length(ts_cols) == 0) {
      stop("No numeric columns found for resilience analysis", call. = FALSE)
    }
    
    message("Auto-detected time series columns for resilience analysis: ", 
            paste(ts_cols, collapse = ", "))
  }
  
  # Apply scaling using existing TSnetworks scaling utilities
  if (scaling_method != "none") {
    if (scaling_method == "auto") {
      # Use existing auto-detection
      scaling_method <- auto_detect_scaling(data[[ts_cols[1]]])
      message("Auto-detected scaling method: ", scaling_method)
    }
    
    # Apply scaling to time series columns
    scaled_data <- scale_dataframe(data, cols = ts_cols, method = scaling_method)
  } else {
    scaled_data <- data
  }
  
  # Calculate resilience metrics by capacity type
  result_data <- scaled_data
  
  if ("absorptive" %in% capacity_types) {
    result_data <- .calculate_absorptive_metrics(result_data, ts_cols, id_col, 
                                                window_width, baseline_method, na_action)
  }
  
  if ("restorative" %in% capacity_types) {
    result_data <- .calculate_restorative_metrics(result_data, ts_cols, id_col, 
                                                 window_width, baseline_method, na_action)
  }
  
  if ("adaptive" %in% capacity_types) {
    result_data <- .calculate_adaptive_metrics(result_data, ts_cols, id_col, 
                                              window_width, baseline_method, na_action)
  }
  
  # Add metadata attributes for compatibility
  attr(result_data, "resilience_analysis") <- list(
    capacity_types = capacity_types,
    ts_cols = ts_cols,
    window_width = window_width,
    scaling_method = scaling_method,
    baseline_method = baseline_method,
    analysis_timestamp = Sys.time()
  )
  
  class(result_data) <- c("resilience_data", class(result_data))
  return(result_data)
}

#' Calculate Absorptive Capacity Metrics
#'
#' @description
#' Calculates absorptive capacity metrics that measure a system's ability
#' to absorb disturbances without changing state. Compatible with all
#' existing TSnetworks functions.
#'
#' @param data Data frame with time series data
#' @param ts_cols Time series column names
#' @param id_col Grouping column (optional)
#' @param window_width Rolling window size
#' @param baseline_method Method for baseline calculation
#' @param na_action How to handle missing values
#'
#' @return Data frame with absorptive capacity columns added
#'
#' @examples
#' data(saqrsteps)
#' absorptive_data <- calculate_absorptive_capacity(saqrsteps, ts_cols = "Steps")
#' head(absorptive_data)
#'
#' @export
calculate_absorptive_capacity <- function(data,
                                         ts_cols = NULL,
                                         id_col = NULL,
                                         window_width = 7L,
                                         baseline_method = "auto",
                                         na_action = c("pass", "omit", "fail")) {
  
  na_action <- match.arg(na_action)
  
  # Auto-detect columns if needed
  if (is.null(ts_cols)) {
    exclude_cols <- if (!is.null(id_col)) id_col else character(0)
    candidate_cols <- setdiff(names(data), exclude_cols)
    ts_cols <- candidate_cols[sapply(data[candidate_cols], is.numeric)]
  }
  
  result_data <- .calculate_absorptive_metrics(data, ts_cols, id_col, 
                                              window_width, baseline_method, na_action)
  
  attr(result_data, "capacity_type") <- "absorptive"
  return(result_data)
}

#' Calculate Restorative Capacity Metrics
#'
#' @description
#' Calculates restorative capacity metrics that measure a system's ability
#' to recover from disturbances. Compatible with all existing TSnetworks functions.
#'
#' @param data Data frame with time series data
#' @param ts_cols Time series column names
#' @param id_col Grouping column (optional)
#' @param window_width Rolling window size
#' @param baseline_method Method for baseline calculation
#' @param na_action How to handle missing values
#'
#' @return Data frame with restorative capacity columns added
#'
#' @examples
#' data(saqrsteps)
#' restorative_data <- calculate_restorative_capacity(saqrsteps, ts_cols = "Steps")
#' head(restorative_data)
#'
#' @export
calculate_restorative_capacity <- function(data,
                                          ts_cols = NULL,
                                          id_col = NULL,
                                          window_width = 7L,
                                          baseline_method = "auto",
                                          na_action = c("pass", "omit", "fail")) {
  
  na_action <- match.arg(na_action)
  
  # Auto-detect columns if needed
  if (is.null(ts_cols)) {
    exclude_cols <- if (!is.null(id_col)) id_col else character(0)
    candidate_cols <- setdiff(names(data), exclude_cols)
    ts_cols <- candidate_cols[sapply(data[candidate_cols], is.numeric)]
  }
  
  result_data <- .calculate_restorative_metrics(data, ts_cols, id_col, 
                                               window_width, baseline_method, na_action)
  
  attr(result_data, "capacity_type") <- "restorative"
  return(result_data)
}

#' Calculate Adaptive Capacity Metrics
#'
#' @description
#' Calculates adaptive capacity metrics that measure a system's ability
#' to adapt and evolve in response to disturbances. Compatible with all
#' existing TSnetworks functions.
#'
#' @param data Data frame with time series data
#' @param ts_cols Time series column names
#' @param id_col Grouping column (optional)
#' @param window_width Rolling window size
#' @param baseline_method Method for baseline calculation
#' @param na_action How to handle missing values
#'
#' @return Data frame with adaptive capacity columns added
#'
#' @examples
#' data(saqrsteps)
#' adaptive_data <- calculate_adaptive_capacity(saqrsteps, ts_cols = "Steps")
#' head(adaptive_data)
#'
#' @export
calculate_adaptive_capacity <- function(data,
                                       ts_cols = NULL,
                                       id_col = NULL,
                                       window_width = 7L,
                                       baseline_method = "auto",
                                       na_action = c("pass", "omit", "fail")) {
  
  na_action <- match.arg(na_action)
  
  # Auto-detect columns if needed
  if (is.null(ts_cols)) {
    exclude_cols <- if (!is.null(id_col)) id_col else character(0)
    candidate_cols <- setdiff(names(data), exclude_cols)
    ts_cols <- candidate_cols[sapply(data[candidate_cols], is.numeric)]
  }
  
  result_data <- .calculate_adaptive_metrics(data, ts_cols, id_col, 
                                            window_width, baseline_method, na_action)
  
  attr(result_data, "capacity_type") <- "adaptive"
  return(result_data)
}

# Internal validation function
.validate_resilience_inputs <- function(data, ts_cols, id_col, window_width, capacity_types) {
  
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
  
  # Check window_width
  if (!is.numeric(window_width) || length(window_width) != 1) {
    stop("'window_width' must be a single numeric value", call. = FALSE)
  }
  if (window_width < 2 || window_width != as.integer(window_width)) {
    stop("'window_width' must be an integer >= 2", call. = FALSE)
  }
  
  # Check capacity_types
  valid_capacities <- c("absorptive", "restorative", "adaptive")
  if (!all(capacity_types %in% valid_capacities)) {
    invalid <- capacity_types[!capacity_types %in% valid_capacities]
    stop("Invalid capacity type(s): ", paste(invalid, collapse = ", "),
         "\nValid options: ", paste(valid_capacities, collapse = ", "),
         call. = FALSE)
  }
}

# Internal function to calculate absorptive capacity metrics
.calculate_absorptive_metrics <- function(data, ts_cols, id_col, window_width, baseline_method, na_action) {
  
  result_data <- data
  
  # Determine calculation scope
  use_grouping <- !is.null(id_col) && id_col %in% names(data)
  
  for (col in ts_cols) {
    
    # Initialize result vectors
    vsi_values <- numeric(nrow(data))
    arch_values <- numeric(nrow(data))
    cv_values <- numeric(nrow(data))
    
    if (use_grouping) {
      # Grouped calculation
      groups <- split(data, data[[id_col]], drop = TRUE)
      
      for (group_name in names(groups)) {
        group_data <- groups[[group_name]]
        group_indices <- which(data[[id_col]] == group_name)
        
        if (nrow(group_data) >= window_width) {
          ts_values <- group_data[[col]]
          
          # Calculate metrics for this group
          vsi_values[group_indices] <- .calculate_variance_stability_index(ts_values, window_width, baseline_method)
          arch_values[group_indices] <- .calculate_arch_effects(ts_values, window_width)
          cv_values[group_indices] <- .calculate_adaptive_cv(ts_values, window_width)
        } else {
          # Not enough data for this group
          vsi_values[group_indices] <- NA_real_
          arch_values[group_indices] <- NA_real_
          cv_values[group_indices] <- NA_real_
        }
      }
      
    } else {
      # Global calculation
      ts_values <- data[[col]]
      
      if (length(ts_values) >= window_width) {
        vsi_values <- .calculate_variance_stability_index(ts_values, window_width, baseline_method)
        arch_values <- .calculate_arch_effects(ts_values, window_width)
        cv_values <- .calculate_adaptive_cv(ts_values, window_width)
      } else {
        vsi_values <- rep(NA_real_, length(ts_values))
        arch_values <- rep(NA_real_, length(ts_values))
        cv_values <- rep(NA_real_, length(ts_values))
      }
    }
    
    # Add columns to result
    result_data[[paste0("absorptive_vsi_", col)]] <- vsi_values
    result_data[[paste0("absorptive_arch_", col)]] <- arch_values
    result_data[[paste0("absorptive_cv_", col)]] <- cv_values
  }
  
  return(result_data)
}

# Internal function to calculate restorative capacity metrics
.calculate_restorative_metrics <- function(data, ts_cols, id_col, window_width, baseline_method, na_action) {
  
  result_data <- data
  
  # Determine calculation scope
  use_grouping <- !is.null(id_col) && id_col %in% names(data)
  
  for (col in ts_cols) {
    
    # Initialize result vectors
    recovery_time_values <- numeric(nrow(data))
    eng_resilience_values <- numeric(nrow(data))
    recovery_slope_values <- numeric(nrow(data))
    
    if (use_grouping) {
      # Grouped calculation
      groups <- split(data, data[[id_col]], drop = TRUE)
      
      for (group_name in names(groups)) {
        group_data <- groups[[group_name]]
        group_indices <- which(data[[id_col]] == group_name)
        
        if (nrow(group_data) >= window_width) {
          ts_values <- group_data[[col]]
          
          # Calculate metrics for this group
          recovery_time_values[group_indices] <- .calculate_recovery_time(ts_values, window_width)
          eng_resilience_values[group_indices] <- .calculate_engineering_resilience(ts_values, window_width)
          recovery_slope_values[group_indices] <- .calculate_recovery_slope(ts_values, window_width)
        } else {
          # Not enough data for this group
          recovery_time_values[group_indices] <- NA_real_
          eng_resilience_values[group_indices] <- NA_real_
          recovery_slope_values[group_indices] <- NA_real_
        }
      }
      
    } else {
      # Global calculation
      ts_values <- data[[col]]
      
      if (length(ts_values) >= window_width) {
        recovery_time_values <- .calculate_recovery_time(ts_values, window_width)
        eng_resilience_values <- .calculate_engineering_resilience(ts_values, window_width)
        recovery_slope_values <- .calculate_recovery_slope(ts_values, window_width)
      } else {
        recovery_time_values <- rep(NA_real_, length(ts_values))
        eng_resilience_values <- rep(NA_real_, length(ts_values))
        recovery_slope_values <- rep(NA_real_, length(ts_values))
      }
    }
    
    # Add columns to result
    result_data[[paste0("restorative_recovery_time_", col)]] <- recovery_time_values
    result_data[[paste0("restorative_eng_resilience_", col)]] <- eng_resilience_values
    result_data[[paste0("restorative_recovery_slope_", col)]] <- recovery_slope_values
  }
  
  return(result_data)
}

# Internal function to calculate adaptive capacity metrics
.calculate_adaptive_metrics <- function(data, ts_cols, id_col, window_width, baseline_method, na_action) {
  
  result_data <- data
  
  # Determine calculation scope
  use_grouping <- !is.null(id_col) && id_col %in% names(data)
  
  for (col in ts_cols) {
    
    # Initialize result vectors
    sample_entropy_values <- numeric(nrow(data))
    dfa_values <- numeric(nrow(data))
    capacity_ratio_values <- numeric(nrow(data))
    
    if (use_grouping) {
      # Grouped calculation
      groups <- split(data, data[[id_col]], drop = TRUE)
      
      for (group_name in names(groups)) {
        group_data <- groups[[group_name]]
        group_indices <- which(data[[id_col]] == group_name)
        
        if (nrow(group_data) >= window_width) {
          ts_values <- group_data[[col]]
          
          # Calculate metrics for this group
          sample_entropy_values[group_indices] <- .calculate_sample_entropy(ts_values, window_width)
          dfa_values[group_indices] <- .calculate_dfa_alpha(ts_values, window_width)
          capacity_ratio_values[group_indices] <- .calculate_adaptive_capacity_ratio(ts_values, window_width)
        } else {
          # Not enough data for this group
          sample_entropy_values[group_indices] <- NA_real_
          dfa_values[group_indices] <- NA_real_
          capacity_ratio_values[group_indices] <- NA_real_
        }
      }
      
    } else {
      # Global calculation
      ts_values <- data[[col]]
      
      if (length(ts_values) >= window_width) {
        sample_entropy_values <- .calculate_sample_entropy(ts_values, window_width)
        dfa_values <- .calculate_dfa_alpha(ts_values, window_width)
        capacity_ratio_values <- .calculate_adaptive_capacity_ratio(ts_values, window_width)
      } else {
        sample_entropy_values <- rep(NA_real_, length(ts_values))
        dfa_values <- rep(NA_real_, length(ts_values))
        capacity_ratio_values <- rep(NA_real_, length(ts_values))
      }
    }
    
    # Add columns to result
    result_data[[paste0("adaptive_sample_entropy_", col)]] <- sample_entropy_values
    result_data[[paste0("adaptive_dfa_", col)]] <- dfa_values
    result_data[[paste0("adaptive_capacity_ratio_", col)]] <- capacity_ratio_values
  }
  
  return(result_data)
}

# Core metric calculation functions

.calculate_variance_stability_index <- function(x, window_width, baseline_method) {
  
  # Calculate baseline variance
  if (baseline_method == "auto") {
    baseline_var <- var(x, na.rm = TRUE)
  } else {
    baseline_var <- var(x, na.rm = TRUE)  # Can be extended with other methods
  }
  
  # Calculate rolling variance
  rolling_var <- zoo::rollapply(
    x,
    width = window_width,
    FUN = function(vals) var(vals, na.rm = TRUE),
    partial = FALSE,
    fill = NA,
    align = "right"
  )
  
  # Calculate VSI
  vsi <- 1 - abs(rolling_var - baseline_var) / baseline_var
  vsi[vsi < 0] <- 0
  vsi[is.infinite(vsi) | is.nan(vsi)] <- NA
  
  return(vsi)
}

.calculate_arch_effects <- function(x, window_width) {
  
  # Calculate rolling ARCH test statistic
  arch_stats <- zoo::rollapply(
    x,
    width = window_width,
    FUN = function(vals) {
      if (length(vals) < 4 || all(is.na(vals))) return(NA_real_)
      
      # Fit AR(1) model
      tryCatch({
        ar_fit <- lm(vals[-1] ~ vals[-length(vals)])
        residuals <- residuals(ar_fit)
        
        # Test for heteroscedasticity
        squared_resid <- residuals^2
        if (length(squared_resid) < 3) return(NA_real_)
        
        # Simple ARCH test
        arch_test <- lm(squared_resid[-1] ~ squared_resid[-length(squared_resid)])
        arch_stat <- summary(arch_test)$r.squared
        
        return(arch_stat)
      }, error = function(e) NA_real_)
    },
    partial = FALSE,
    fill = NA,
    align = "right"
  )
  
  return(arch_stats)
}

.calculate_adaptive_cv <- function(x, window_width) {
  
  # Calculate rolling coefficient of variation
  cv_values <- zoo::rollapply(
    x,
    width = window_width,
    FUN = function(vals) {
      if (all(is.na(vals))) return(NA_real_)
      
      mean_val <- mean(vals, na.rm = TRUE)
      sd_val <- sd(vals, na.rm = TRUE)
      
      if (is.na(mean_val) || is.na(sd_val) || mean_val == 0) {
        return(NA_real_)
      }
      
      return(sd_val / abs(mean_val))
    },
    partial = FALSE,
    fill = NA,
    align = "right"
  )
  
  return(cv_values)
}

.calculate_recovery_time <- function(x, window_width) {
  
  # Calculate rolling recovery time based on AR(1) model
  recovery_times <- zoo::rollapply(
    x,
    width = window_width,
    FUN = function(vals) {
      if (length(vals) < 4 || all(is.na(vals))) return(NA_real_)
      
      tryCatch({
        # Fit AR(1) model
        ar_fit <- lm(vals[-1] ~ vals[-length(vals)])
        ar_coef <- coef(ar_fit)[2]
        
        if (is.na(ar_coef) || ar_coef >= 1 || ar_coef <= -1) {
          return(NA_real_)
        }
        
        # Calculate return time
        return_time <- -1 / log(abs(ar_coef))
        return(return_time)
      }, error = function(e) NA_real_)
    },
    partial = FALSE,
    fill = NA,
    align = "right"
  )
  
  return(recovery_times)
}

.calculate_engineering_resilience <- function(x, window_width) {
  
  # Calculate engineering resilience (1/recovery_time)
  recovery_times <- .calculate_recovery_time(x, window_width)
  
  eng_resilience <- ifelse(is.na(recovery_times) | recovery_times <= 0, 
                          NA_real_, 
                          1 / recovery_times)
  
  return(eng_resilience)
}

.calculate_recovery_slope <- function(x, window_width) {
  
  # Calculate rolling recovery slope
  recovery_slopes <- zoo::rollapply(
    x,
    width = window_width,
    FUN = function(vals) {
      if (length(vals) < 3 || all(is.na(vals))) return(NA_real_)
      
      tryCatch({
        # Detect potential disruption (largest negative change)
        diffs <- diff(vals)
        disruption_idx <- which.min(diffs)
        
        if (length(disruption_idx) == 0 || disruption_idx >= length(vals)) {
          return(NA_real_)
        }
        
        # Calculate recovery slope after disruption
        recovery_vals <- vals[(disruption_idx + 1):length(vals)]
        if (length(recovery_vals) < 2) return(NA_real_)
        
        time_points <- 1:length(recovery_vals)
        slope_fit <- lm(recovery_vals ~ time_points)
        slope <- coef(slope_fit)[2]
        
        return(slope)
      }, error = function(e) NA_real_)
    },
    partial = FALSE,
    fill = NA,
    align = "right"
  )
  
  return(recovery_slopes)
}

.calculate_sample_entropy <- function(x, window_width, pattern_length = 2) {
  
  # Calculate rolling sample entropy
  entropy_values <- zoo::rollapply(
    x,
    width = window_width,
    FUN = function(vals) {
      if (length(vals) < pattern_length + 1 || all(is.na(vals))) return(NA_real_)
      
      tryCatch({
        # Simple sample entropy calculation
        tolerance <- 0.2 * sd(vals, na.rm = TRUE)
        if (is.na(tolerance) || tolerance == 0) return(NA_real_)
        
        n <- length(vals)
        
        # Count pattern matches
        matches_m <- 0
        matches_m1 <- 0
        
        for (i in 1:(n - pattern_length)) {
          template <- vals[i:(i + pattern_length - 1)]
          
          for (j in 1:(n - pattern_length)) {
            if (i != j) {
              candidate <- vals[j:(j + pattern_length - 1)]
              
              if (max(abs(template - candidate), na.rm = TRUE) <= tolerance) {
                matches_m <- matches_m + 1
                
                # Check m+1 pattern
                if (i <= n - pattern_length - 1 && j <= n - pattern_length - 1) {
                  template_m1 <- vals[i:(i + pattern_length)]
                  candidate_m1 <- vals[j:(j + pattern_length)]
                  
                  if (max(abs(template_m1 - candidate_m1), na.rm = TRUE) <= tolerance) {
                    matches_m1 <- matches_m1 + 1
                  }
                }
              }
            }
          }
        }
        
        if (matches_m == 0 || matches_m1 == 0) return(NA_real_)
        
        sample_entropy <- -log(matches_m1 / matches_m)
        return(sample_entropy)
      }, error = function(e) NA_real_)
    },
    partial = FALSE,
    fill = NA,
    align = "right"
  )
  
  return(entropy_values)
}

.calculate_dfa_alpha <- function(x, window_width) {
  
  # Calculate rolling DFA alpha (simplified version)
  dfa_values <- zoo::rollapply(
    x,
    width = window_width,
    FUN = function(vals) {
      if (length(vals) < 8 || all(is.na(vals))) return(NA_real_)
      
      tryCatch({
        # Create cumulative sum profile
        profile <- cumsum(vals - mean(vals, na.rm = TRUE))
        
        # Use a few scales for simplified DFA
        scales <- c(4, 8, min(16, floor(length(vals)/4)))
        scales <- scales[scales <= length(vals)/2]
        
        if (length(scales) < 2) return(NA_real_)
        
        fluctuations <- numeric(length(scales))
        
        for (i in seq_along(scales)) {
          scale <- scales[i]
          n_segments <- floor(length(profile) / scale)
          
          if (n_segments < 1) {
            fluctuations[i] <- NA
            next
          }
          
          segment_flucts <- numeric(n_segments)
          
          for (j in 1:n_segments) {
            start_idx <- (j - 1) * scale + 1
            end_idx <- j * scale
            
            segment <- profile[start_idx:end_idx]
            x_vals <- 1:scale
            
            # Linear detrending
            fit <- lm(segment ~ x_vals)
            detrended <- residuals(fit)
            
            segment_flucts[j] <- sqrt(mean(detrended^2))
          }
          
          fluctuations[i] <- sqrt(mean(segment_flucts^2))
        }
        
        # Calculate scaling exponent
        valid_idx <- !is.na(fluctuations) & fluctuations > 0
        if (sum(valid_idx) < 2) return(NA_real_)
        
        log_scales <- log(scales[valid_idx])
        log_flucts <- log(fluctuations[valid_idx])
        
        fit <- lm(log_flucts ~ log_scales)
        alpha <- coef(fit)[2]
        
        return(alpha)
      }, error = function(e) NA_real_)
    },
    partial = FALSE,
    fill = NA,
    align = "right"
  )
  
  return(dfa_values)
}

.calculate_adaptive_capacity_ratio <- function(x, window_width) {
  
  # Calculate rolling adaptive capacity ratio
  capacity_ratios <- zoo::rollapply(
    x,
    width = window_width,
    FUN = function(vals) {
      if (length(vals) < 3 || all(is.na(vals))) return(NA_real_)
      
      tryCatch({
        # Use first differences as proxy for disturbance
        diffs <- diff(vals)
        
        # Calculate response variance (variance of the time series)
        response_var <- var(vals, na.rm = TRUE)
        
        # Calculate driver variance (variance of changes)
        driver_var <- var(diffs, na.rm = TRUE)
        
        if (is.na(response_var) || is.na(driver_var) || driver_var == 0) {
          return(NA_real_)
        }
        
        # Adaptive capacity ratio
        ac_ratio <- response_var / (driver_var + 1e-10)  # Small constant to avoid division by zero
        
        return(ac_ratio)
      }, error = function(e) NA_real_)
    },
    partial = FALSE,
    fill = NA,
    align = "right"
  )
  
  return(capacity_ratios)
}

#' Print method for resilience_data objects
#' @param x A resilience_data object
#' @param ... Additional arguments (unused)
#' @export
print.resilience_data <- function(x, ...) {
  cat("Resilience Analysis Data\n")
  cat("========================\n")
  
  # Get resilience analysis info
  resilience_info <- attr(x, "resilience_analysis")
  
  if (!is.null(resilience_info)) {
    cat("Capacity types analyzed:", paste(resilience_info$capacity_types, collapse = ", "), "\n")
    cat("Time series columns:", paste(resilience_info$ts_cols, collapse = ", "), "\n")
    cat("Window width:", resilience_info$window_width, "\n")
    cat("Scaling method:", resilience_info$scaling_method, "\n")
    cat("Analysis timestamp:", format(resilience_info$analysis_timestamp), "\n")
  }
  
  cat("Dimensions:", nrow(x), "x", ncol(x), "\n")
  
  # Identify resilience columns
  resilience_cols <- grep("^(absorptive|restorative|adaptive)_", names(x), value = TRUE)
  if (length(resilience_cols) > 0) {
    cat("Resilience metrics columns:", length(resilience_cols), "\n")
    if (length(resilience_cols) <= 10) {
      cat("  ", paste(resilience_cols, collapse = ", "), "\n")
    } else {
      cat("  ", paste(resilience_cols[1:10], collapse = ", "), "...\n")
    }
  }
  
  cat("\nFirst few rows:\n")
  print(head(x))
}
