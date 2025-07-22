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
  
  # CRITICAL FIX: Proper baseline variance calculation using robust statistics
  if (baseline_method == "auto") {
    # Use robust estimation for baseline - median absolute deviation scaled to variance
    baseline_var <- (1.4826 * mad(x, na.rm = TRUE))^2
    
    # If MAD is zero (constant series), handle gracefully
    if (baseline_var == 0 || is.na(baseline_var)) {
      # For constant series, set baseline variance to small positive value
      # This avoids division by zero and gives reasonable VSI = 1 for constant series
      baseline_var <- 1e-10
    }
  } else {
    baseline_var <- var(x, na.rm = TRUE)
    # Handle constant series in non-auto method too
    if (baseline_var == 0 || is.na(baseline_var)) {
      baseline_var <- 1e-10
    }
  }
  
  # CRITICAL CHECK: Ensure baseline variance is valid
  if (baseline_var <= 0 || is.na(baseline_var)) {
    warning("Invalid baseline variance for VSI calculation")
    return(rep(NA_real_, length(x)))
  }
  
  # Calculate rolling variance with proper edge handling
  rolling_var <- zoo::rollapply(
    x,
    width = window_width,
    FUN = function(vals) {
      if (length(vals) < max(3, window_width/4) || all(is.na(vals))) return(NA_real_)  # Reduced minimum requirement
      var(vals, na.rm = TRUE)
    },
    partial = TRUE,
    fill = NA,
    align = "center"  # Changed from "right" for better temporal resolution
  )
  
  # SCIENTIFIC CORRECTION: Proper VSI formula from literature
  # VSI = 1 - |σ²(t) - σ²_baseline| / (σ²_baseline + σ²(t))
  # This prevents division by zero and provides proper normalization
  vsi <- 1 - abs(rolling_var - baseline_var) / (baseline_var + rolling_var + 1e-10)
  
  # Apply bounds and handle edge cases
  vsi[vsi < 0] <- 0
  vsi[vsi > 1] <- 1
  vsi[is.infinite(vsi) | is.nan(vsi)] <- NA
  
  return(vsi)
}

.calculate_arch_effects <- function(x, window_width) {
  
  # Calculate rolling ARCH LM test statistic following Engle (1982)
  arch_stats <- zoo::rollapply(
    x,
    width = window_width,
    FUN = function(vals) {
      if (length(vals) < max(4, window_width/4) || all(is.na(vals))) return(NA_real_)  # Reduced from 6 to 4
      
      tryCatch({
        # SCIENTIFIC CORRECTION: Proper ARCH test procedure
        
        # Step 1: Fit mean equation (AR(1) or constant mean)
        n <- length(vals)
        vals_centered <- vals - mean(vals, na.rm = TRUE)
        
        # Try AR(1) first, fall back to constant mean if necessary
        if (n >= 4) {
          # AR(1) regression: x_t = φ*x_{t-1} + ε_t
          y <- vals_centered[-1]
          x_lag <- vals_centered[-n]
          
          if (sd(x_lag) > 1e-8) {  # Check for variation
            ar_fit <- lm(y ~ x_lag)
            residuals <- residuals(ar_fit)
          } else {
            residuals <- vals_centered
          }
        } else {
          residuals <- vals_centered
        }
        
        # Step 2: ARCH test on squared residuals
        squared_resid <- residuals^2
        n_resid <- length(squared_resid)
        
        if (n_resid < 4) return(NA_real_)
        
        # ARCH(1) test: σ²_t = α₀ + α₁*ε²_{t-1}
        y_arch <- squared_resid[-1]
        x_arch <- squared_resid[-n_resid]
        
        # Check for sufficient variation
        if (sd(x_arch) < 1e-8 || sd(y_arch) < 1e-8) {
          return(0.0)  # No heteroscedasticity if no variation
        }
        
        # Fit ARCH regression
        arch_reg <- lm(y_arch ~ x_arch)
        
        # LM test statistic: n*R² ~ χ²(1)
        r_squared <- summary(arch_reg)$r.squared
        lm_statistic <- (n_resid - 1) * r_squared
        
        # Return p-value (more interpretable than raw statistic)
        p_value <- 1 - pchisq(lm_statistic, df = 1)
        
        # Convert to arch effect measure (0 = no ARCH, 1 = strong ARCH)
        arch_effect <- 1 - p_value
        
        return(arch_effect)
        
      }, error = function(e) NA_real_)
    },
    partial = TRUE,
    fill = NA,
    align = "center"
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
  
  # Calculate rolling recovery time using proper AR methodology
  recovery_times <- zoo::rollapply(
    x,
    width = window_width,
    FUN = function(vals) {
      if (length(vals) < max(4, window_width/4) || all(is.na(vals))) return(NA_real_)  # Reduced from 8
      
      tryCatch({
        # SCIENTIFIC CORRECTION: Proper autoregressive analysis
        
        # Step 1: Center the series for stationarity
        vals_centered <- vals - mean(vals, na.rm = TRUE)
        n <- length(vals_centered)
        
        # Check for sufficient variation
        if (sd(vals_centered, na.rm = TRUE) < 1e-6) return(NA_real_)
        
        # Step 2: Fit AR(1) model with proper validation
        y <- vals_centered[-1]
        x_lag <- vals_centered[-n]
        
        # Check for variation in lagged values
        if (sd(x_lag, na.rm = TRUE) < 1e-6) return(NA_real_)
        
        # Fit AR(1) using robust methods
        ar_fit <- tryCatch({
          MASS::rlm(y ~ x_lag, method = "M", maxit = 50)
        }, warning = function(w) {
          # If convergence issues, try simpler method
          tryCatch({
            MASS::rlm(y ~ x_lag, method = "MM", maxit = 100)
          }, error = function(e) {
            # Final fallback to OLS
            lm(y ~ x_lag)
          })
        }, error = function(e) {
          # Final fallback to OLS
          lm(y ~ x_lag)
        })
        
        ar_coef <- coef(ar_fit)[2]  # Slope coefficient (φ)
        
        # Step 3: Validate AR coefficient
        if (is.na(ar_coef)) return(NA_real_)
        
        # Check for stationarity constraint: |φ| < 1
        if (abs(ar_coef) >= 0.99) {
          # Near unit root - extremely slow recovery
          return(Inf)
        }
        
        if (abs(ar_coef) < 0.01) {
          # Very fast recovery (white noise)
          return(1.0)
        }
        
        # Step 4: Calculate half-life (proper formula)
        # For AR(1): x_t = φ*x_{t-1} + ε_t
        # Half-life = ln(0.5) / ln(|φ|)
        half_life <- log(0.5) / log(abs(ar_coef))
        
        # Validate result and handle infinite values
        if (is.infinite(half_life) || is.nan(half_life) || half_life <= 0) {
          return(NA_real_)
        }
        
        # Cap at reasonable maximum (e.g., 5 times window width)
        max_recovery_time <- window_width * 5
        if (half_life > max_recovery_time) {
          return(max_recovery_time)
        }
        
        return(half_life)
        
      }, error = function(e) NA_real_)
    },
    partial = TRUE,
    fill = NA,
    align = "center"
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
  
  # Use the scientifically corrected recovery slopes calculation
  return(.calculate_recovery_slopes(x, window_width, threshold_sd = 2.0))
}

.calculate_sample_entropy <- function(x, window_width, pattern_length = 2) {
  
  # Calculate rolling sample entropy following Richman & Moorman (2000)
  entropy_values <- zoo::rollapply(
    x,
    width = window_width,
    FUN = function(vals) {
      if (length(vals) < pattern_length + max(4, window_width/4) || all(is.na(vals))) return(NA_real_)  # Reduced from 10
      
      tryCatch({
        # SCIENTIFIC CORRECTION: Proper sample entropy algorithm
        
        # Step 1: Standardize the data for scale-invariant analysis
        vals_std <- as.numeric(scale(vals))
        n <- length(vals_std)
        
        # Step 2: Calculate tolerance (critical parameter)
        # Use 0.2 * SD as recommended by Richman & Moorman (2000)
        tolerance <- 0.2 * sd(vals_std, na.rm = TRUE)
        
        # Alternative tolerance methods for robustness
        if (is.na(tolerance) || tolerance <= 0) {
          # Use median absolute deviation if SD fails
          tolerance <- 0.2 * 1.4826 * mad(vals_std, na.rm = TRUE)
          if (tolerance <= 0) return(NA_real_)
        }
        
        # Step 3: Count pattern matches for m and m+1
        matches_m <- 0
        matches_m1 <- 0
        total_comparisons <- 0
        
        # Iterate through all possible templates
        for (i in 1:(n - pattern_length)) {
          template_m <- vals_std[i:(i + pattern_length - 1)]
          
          # For m+1 patterns
          if (i <= n - pattern_length) {
            template_m1 <- vals_std[i:(i + pattern_length)]
          } else {
            next
          }
          
          # Compare with all other patterns (excluding self-matches)
          for (j in 1:(n - pattern_length)) {
            if (i == j) next  # Skip self-match
            
            candidate_m <- vals_std[j:(j + pattern_length - 1)]
            
            # Check if m-patterns match using Chebyshev distance
            max_diff_m <- max(abs(template_m - candidate_m))
            
            if (max_diff_m <= tolerance) {
              matches_m <- matches_m + 1
              
              # Check m+1 pattern if both templates are valid
              if (j <= n - pattern_length) {
                candidate_m1 <- vals_std[j:(j + pattern_length)]
                max_diff_m1 <- max(abs(template_m1 - candidate_m1))
                
                if (max_diff_m1 <= tolerance) {
                  matches_m1 <- matches_m1 + 1
                }
              }
            }
            
            total_comparisons <- total_comparisons + 1
          }
        }
        
        # Step 4: Calculate conditional probabilities
        if (matches_m == 0 || matches_m1 == 0 || total_comparisons == 0) {
          return(NA_real_)
        }
        
        # Relative template matching probabilities
        phi_m <- matches_m / total_comparisons
        phi_m1 <- matches_m1 / total_comparisons
        
        # Step 5: Sample entropy calculation
        if (phi_m1 > 0 && phi_m > 0) {
          sample_entropy <- -log(phi_m1 / phi_m)
          
          # Validate result (sample entropy should be non-negative)
          if (sample_entropy < 0 || is.infinite(sample_entropy)) {
            return(NA_real_)
          }
          
          return(sample_entropy)
        } else {
          return(NA_real_)
        }
        
      }, error = function(e) NA_real_)
    },
    partial = TRUE,
    fill = NA,
    align = "center"
  )
  
  return(entropy_values)
}

.calculate_dfa_alpha <- function(x, window_width) {
  
  # Calculate rolling DFA alpha using proper methodology (Peng et al. 1994)
  dfa_values <- zoo::rollapply(
    x,
    width = window_width,
    FUN = function(vals) {
      if (length(vals) < max(8, window_width/2) || all(is.na(vals))) return(NA_real_)  # Reduced from 16
      
      tryCatch({
        # SCIENTIFIC CORRECTION: Proper DFA implementation
        
        # Step 1: Remove mean and integrate the series
        vals_centered <- vals - mean(vals, na.rm = TRUE)
        integrated_series <- cumsum(vals_centered)
        n <- length(integrated_series)
        
        # Step 2: Define scale range (minimum 4 points per window)
        min_scale <- 4
        max_scale <- floor(n / 4)  # Ensure at least 4 windows
        
        if (max_scale <= min_scale) return(NA_real_)
        
        # Use logarithmically spaced scales for proper scaling analysis
        num_scales <- min(10, max_scale - min_scale + 1)  # Limit for efficiency
        scales <- unique(round(exp(seq(log(min_scale), log(max_scale), length.out = num_scales))))
        
        # Step 3: Calculate fluctuation function F(n) for each scale
        fluctuations <- numeric(length(scales))
        
        for (i in seq_along(scales)) {
          scale <- scales[i]
          num_windows <- floor(n / scale)
          
          if (num_windows < 2) {
            fluctuations[i] <- NA
            next
          }
          
          # Calculate local trends and detrended fluctuations
          window_variances <- numeric(num_windows)
          
          for (j in 1:num_windows) {
            start_idx <- (j - 1) * scale + 1
            end_idx <- j * scale
            window_data <- integrated_series[start_idx:end_idx]
            window_indices <- 1:scale
            
            # Linear detrending (DFA-1)
            if (scale >= 3) {
              linear_fit <- lm(window_data ~ window_indices)
              detrended <- residuals(linear_fit)
              window_variances[j] <- mean(detrended^2)
            } else {
              window_variances[j] <- var(window_data)
            }
          }
          
          # Root mean square fluctuation
          fluctuations[i] <- sqrt(mean(window_variances, na.rm = TRUE))
        }
        
        # Step 4: Calculate scaling exponent α from log-log regression
        valid_points <- !is.na(fluctuations) & fluctuations > 0
        
        if (sum(valid_points) < 3) return(NA_real_)
        
        log_scales <- log10(scales[valid_points])
        log_fluctuations <- log10(fluctuations[valid_points])
        
        # Robust linear regression for scaling exponent
        tryCatch({
          # Try robust regression first
          fit <- tryCatch({
            MASS::rlm(log_fluctuations ~ log_scales, method = "M", maxit = 50)
          }, warning = function(w) {
            # If convergence warning, try with different settings
            MASS::rlm(log_fluctuations ~ log_scales, method = "MM", maxit = 100)
          })
          
          alpha <- coef(fit)[2]  # Slope = scaling exponent
          
          # Validate alpha is in reasonable range
          if (alpha < 0.1 || alpha > 2.0) return(NA_real_)
          
          return(alpha)
        }, error = function(e) {
          # Fallback to ordinary least squares
          tryCatch({
            fit <- lm(log_fluctuations ~ log_scales)
            alpha <- coef(fit)[2]
            if (is.na(alpha) || alpha < 0.1 || alpha > 2.0) return(NA_real_)
            return(alpha)
          }, error = function(e2) NA_real_)
        })
        
      }, error = function(e) NA_real_)
    },
    partial = TRUE,
    fill = NA,
    align = "center"
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

.calculate_recovery_slopes <- function(x, window_width, threshold_sd = 2.0) {
  
  # Calculate rolling recovery slopes using proper resilience methodology
  recovery_slopes <- zoo::rollapply(
    x,
    width = window_width,
    FUN = function(vals) {
      if (length(vals) < max(6, window_width/3) || all(is.na(vals))) return(NA_real_)  # Reduced from 10
      
      tryCatch({
        # SCIENTIFIC CORRECTION: Proper recovery analysis
        
        # Step 1: Detect perturbations/shocks using robust statistics
        vals_centered <- vals - median(vals, na.rm = TRUE)
        mad_threshold <- threshold_sd * mad(vals_centered, na.rm = TRUE)
        
        if (mad_threshold == 0) {
          # Fall back to standard deviation if MAD is zero
          sd_threshold <- threshold_sd * sd(vals_centered, na.rm = TRUE)
          shock_indices <- which(abs(vals_centered) > sd_threshold)
        } else {
          shock_indices <- which(abs(vals_centered) > mad_threshold)
        }
        
        if (length(shock_indices) == 0) return(NA_real_)
        
        # Step 2: Analyze recovery trajectories
        recovery_slopes_vec <- numeric()
        
        for (shock_idx in shock_indices) {
          # Define recovery window (minimum 4 points after shock)
          recovery_start <- shock_idx + 1
          recovery_end <- min(shock_idx + 8, length(vals))  # Up to 8 points for recovery
          
          if (recovery_end - recovery_start < 3) next  # Need at least 4 points
          
          # Extract recovery trajectory
          recovery_window <- recovery_start:recovery_end
          recovery_values <- vals[recovery_window]
          time_points <- 1:length(recovery_values)
          
          # Check for valid recovery data
          if (all(is.na(recovery_values)) || sd(recovery_values) == 0) next
          
          # Step 3: Fit recovery model using robust regression
          # Model: exponential decay toward equilibrium
          baseline <- median(vals, na.rm = TRUE)
          initial_displacement <- vals[shock_idx] - baseline
          
          # Linear approximation of recovery in log space (if appropriate)
          if (all(recovery_values > 0) && initial_displacement != 0) {
            # Try exponential recovery model
            normalized_values <- abs(recovery_values - baseline) / abs(initial_displacement)
            normalized_values[normalized_values <= 0] <- 1e-6
            
            log_values <- log(normalized_values)
            
            # Robust linear regression on log-transformed data with fallbacks
            recovery_fit <- tryCatch({
              MASS::rlm(log_values ~ time_points, method = "M", maxit = 50)
            }, warning = function(w) {
              tryCatch({
                MASS::rlm(log_values ~ time_points, method = "MM", maxit = 100)
              }, error = function(e) lm(log_values ~ time_points))
            }, error = function(e) lm(log_values ~ time_points))
            
            recovery_slope <- coef(recovery_fit)[2]
          } else {
            # Simple linear recovery with fallbacks
            recovery_fit <- tryCatch({
              MASS::rlm(recovery_values ~ time_points, method = "M", maxit = 50)
            }, warning = function(w) {
              tryCatch({
                MASS::rlm(recovery_values ~ time_points, method = "MM", maxit = 100)
              }, error = function(e) lm(recovery_values ~ time_points))
            }, error = function(e) lm(recovery_values ~ time_points))
            
            recovery_slope <- coef(recovery_fit)[2]
          }
          
          # Validate slope (should indicate recovery toward baseline)
          if (is.finite(recovery_slope)) {
            recovery_slopes_vec <- c(recovery_slopes_vec, recovery_slope)
          }
        }
        
        # Step 4: Aggregate recovery slopes
        if (length(recovery_slopes_vec) == 0) return(NA_real_)
        
        # Use median recovery slope for robustness
        median_recovery_slope <- median(recovery_slopes_vec, na.rm = TRUE)
        
        # Ensure slope indicates recovery (negative for return to baseline)
        # Convert to resilience measure (faster recovery = higher resilience)
        resilience_score <- -median_recovery_slope  # Negative slope becomes positive resilience
        
        return(resilience_score)
        
      }, error = function(e) NA_real_)
    },
    partial = TRUE,
    fill = NA,
    align = "center"
  )
  
  return(recovery_slopes)
}

#' Print method for resilience_data objects
#' @param x resilience_data object
#' @param ... additional arguments
#' @export
print.resilience_data <- function(x, ...) {
  cat("🎯 Resilience Analysis Results\n")
  cat("==============================\n\n")
  
  # Get attributes
  analysis_info <- attr(x, "resilience_analysis")
  if (!is.null(analysis_info)) {
    cat("📋 Analysis Configuration:\n")
    cat("  • Capacity types: ", paste(analysis_info$capacity_types, collapse = ", "), "\n")
    cat("  • Time series columns: ", paste(analysis_info$ts_cols, collapse = ", "), "\n") 
    cat("  • Window width: ", analysis_info$window_width, "\n")
    cat("  • Scaling method: ", analysis_info$scaling_method, "\n\n")
  }
  
  # Data dimensions
  cat("📊 Data Summary:\n")
  cat("  • Rows: ", nrow(x), "\n")
  cat("  • Total columns: ", ncol(x), "\n")
  
  # Find resilience columns
  resilience_cols <- grep("^(absorptive|restorative|adaptive)_", names(x), value = TRUE)
  cat("  • Resilience metrics: ", length(resilience_cols), "\n\n")
  
  # Show metrics by category
  if (length(resilience_cols) > 0) {
    abs_cols <- grep("^absorptive_", resilience_cols, value = TRUE)
    rest_cols <- grep("^restorative_", resilience_cols, value = TRUE) 
    adapt_cols <- grep("^adaptive_", resilience_cols, value = TRUE)
    
    if (length(abs_cols) > 0) {
      cat("🛡️  Absorptive Capacity (", length(abs_cols), " metrics):\n")
      cat("    ", paste(gsub("absorptive_|_.*", "", abs_cols), collapse = ", "), "\n\n")
    }
    if (length(rest_cols) > 0) {
      cat("🔄 Restorative Capacity (", length(rest_cols), " metrics):\n")
      cat("    ", paste(gsub("restorative_|_.*", "", rest_cols), collapse = ", "), "\n\n")
    }
    if (length(adapt_cols) > 0) {
      cat("🔀 Adaptive Capacity (", length(adapt_cols), " metrics):\n")
      cat("    ", paste(gsub("adaptive_|_.*", "", adapt_cols), collapse = ", "), "\n\n")
    }
  }
  
  # Count non-NA values
  if (length(resilience_cols) > 0) {
    non_na_counts <- sapply(resilience_cols, function(col) sum(!is.na(x[[col]])))
    coverage <- round(mean(non_na_counts) / nrow(x) * 100, 1)
    cat("📈 Data Coverage: ", coverage, "% (avg. non-NA values)\n\n")
  }
  
  cat("💡 Next Steps:\n")
  cat("  • Use as_clean_dataframe(result) for a clean data.frame\n")
  cat("  • Use summary(result) for detailed statistics\n")
  cat("  • Access metrics: result$absorptive_vsi_[column_name]\n")
  
  invisible(x)
}

#' Summary method for resilience_data objects
#' @param object resilience_data object
#' @param ... additional arguments
#' @export
summary.resilience_data <- function(object, ...) {
  cat("📊 Resilience Analysis Summary\n")
  cat("=============================\n\n")
  
  # Get analysis info
  analysis_info <- attr(object, "resilience_analysis")
  if (!is.null(analysis_info)) {
    cat("🔧 Configuration:\n")
    cat("  Window width:", analysis_info$window_width, "\n")
    cat("  Scaling method:", analysis_info$scaling_method, "\n")
    cat("  Time series:", paste(analysis_info$ts_cols, collapse = ", "), "\n\n")
  }
  
  # Find resilience columns
  resilience_cols <- grep("^(absorptive|restorative|adaptive)_", names(object), value = TRUE)
  
  if (length(resilience_cols) == 0) {
    cat("No resilience metrics found.\n")
    return(invisible(object))
  }
  
  # Calculate summary statistics for each metric
  cat("📈 Metric Statistics:\n")
  cat("===================\n\n")
  
  for (col in resilience_cols) {
    values <- object[[col]]
    non_na_count <- sum(!is.na(values))
    total_count <- length(values)
    coverage <- round(non_na_count / total_count * 100, 1)
    
    cat("🔹", col, "\n")
    cat("   Coverage: ", coverage, "% (", non_na_count, "/", total_count, ")\n")
    
    if (non_na_count > 0) {
      valid_values <- values[!is.na(values)]
      cat("   Range: [", round(min(valid_values), 3), ", ", round(max(valid_values), 3), "]\n")
      cat("   Mean: ", round(mean(valid_values), 3), " ± ", round(sd(valid_values), 3), "\n")
      cat("   Median: ", round(median(valid_values), 3), "\n")
    } else {
      cat("   No valid values\n")
    }
    cat("\n")
  }
  
  # Overall summary
  all_coverages <- sapply(resilience_cols, function(col) {
    sum(!is.na(object[[col]])) / length(object[[col]]) * 100
  })
  
  cat("🎯 Overall Summary:\n")
  cat("   Average coverage: ", round(mean(all_coverages), 1), "%\n")
  cat("   Best metric: ", names(which.max(all_coverages)), " (", round(max(all_coverages), 1), "%)\n")
  cat("   Worst metric: ", names(which.min(all_coverages)), " (", round(min(all_coverages), 1), "%)\n")
  
  invisible(object)
}

#' Convert resilience results to clean data frame
#'
#' @description
#' Converts resilience analysis results to a clean data frame without 
#' special attributes or classes for easier manipulation.
#'
#' @param x resilience_data object from calculate_resilience_metrics()
#' @param include_original Whether to include original time series columns (default: TRUE)
#' @param include_scaling Whether to include scaling attributes as columns (default: FALSE)
#'
#' @return Clean data frame with resilience metrics
#' @export
as_clean_dataframe <- function(x, include_original = TRUE, include_scaling = FALSE) {
  
  if (!inherits(x, "resilience_data")) {
    stop("Input must be a resilience_data object")
  }
  
  # Get analysis info
  analysis_info <- attr(x, "resilience_analysis")
  
  # Create base data frame by manually extracting columns
  col_names <- names(x)
  result <- data.frame(row.names = rownames(x))
  
  # Copy each column manually to avoid any method dispatch issues
  for (col in col_names) {
    result[[col]] <- x[[col]]
  }
  
  # Optionally remove original time series columns
  if (!include_original && !is.null(analysis_info$ts_cols)) {
    orig_cols <- analysis_info$ts_cols
    orig_cols <- orig_cols[orig_cols %in% names(result)]
    if (length(orig_cols) > 0) {
      result <- result[, !names(result) %in% orig_cols, drop = FALSE]
    }
  }
  
  # Optionally add scaling info as metadata columns
  if (include_scaling && !is.null(analysis_info)) {
    result$analysis_window_width <- analysis_info$window_width
    result$analysis_scaling_method <- analysis_info$scaling_method
    result$analysis_timestamp <- analysis_info$analysis_timestamp
  }
  
  return(result)
}
