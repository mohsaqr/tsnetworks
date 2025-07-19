#' Unified Scaling Functions for Time Series Data
#'
#' @description
#' Provides standardized scaling methods that can be used across all functions
#' in the tsnetworks package. These scaling methods help normalize data for
#' better analysis results and comparability across different time series.
#'
#' @param x Numeric vector or data frame column to scale
#' @param method Scaling method: "none", "center", "standardize", "minmax", "iqr", "robust", "quantile"
#' @param center_value Custom center value (used with "center" method)
#' @param scale_value Custom scale value (used with custom scaling)
#' @param quantile_range Quantile range for "quantile" method (default: c(0.05, 0.95))
#' @param robust_center Whether to use median instead of mean for centering
#' @param na_action How to handle NA values: "omit", "interpolate", "mean", "median"
#'
#' @return List containing scaled data and scaling parameters
#'
#' @examples
#' # Basic scaling examples
#' data <- c(1, 5, 10, 15, 20, 100)
#' 
#' # Standardize (z-score)
#' scaled_std <- apply_scaling(data, method = "standardize")
#' print(scaled_std$data)
#' 
#' # Min-max scaling to [0,1]
#' scaled_minmax <- apply_scaling(data, method = "minmax")
#' print(scaled_minmax$data)
#' 
#' # Robust scaling using IQR
#' scaled_robust <- apply_scaling(data, method = "iqr")
#' print(scaled_robust$data)
#'
#' @importFrom stats median quantile sd IQR na.omit
#' @export
apply_scaling <- function(x, 
                          method = c("none", "center", "standardize", "minmax", "iqr", "robust", "quantile"),
                          center_value = NULL,
                          scale_value = NULL,
                          quantile_range = c(0.05, 0.95),
                          robust_center = FALSE,
                          na_action = c("omit", "interpolate", "mean", "median")) {
  
  method <- match.arg(method)
  na_action <- match.arg(na_action)
  
  # Handle input validation
  if (!is.numeric(x)) {
    stop("Input 'x' must be numeric")
  }
  
  original_x <- x
  original_length <- length(x)
  
  # Handle NA values
  if (any(is.na(x))) {
    if (na_action == "omit") {
      x <- x[!is.na(x)]
      if (length(x) == 0) {
        stop("All values are NA after omitting")
      }
    } else if (na_action == "interpolate") {
      if (requireNamespace("zoo", quietly = TRUE)) {
        x <- zoo::na.approx(x, na.rm = FALSE)
        if (any(is.na(x))) {
          x <- zoo::na.locf(x, na.rm = FALSE)
        }
        if (any(is.na(x))) {
          x <- zoo::na.locf(x, na.rm = FALSE, fromLast = TRUE)
        }
      } else {
        # Fallback to mean imputation
        x[is.na(x)] <- mean(x, na.rm = TRUE)
      }
    } else if (na_action == "mean") {
      x[is.na(x)] <- mean(x, na.rm = TRUE)
    } else if (na_action == "median") {
      x[is.na(x)] <- median(x, na.rm = TRUE)
    }
  }
  
  # Store original statistics
  original_stats <- list(
    mean = mean(original_x, na.rm = TRUE),
    median = median(original_x, na.rm = TRUE),
    sd = sd(original_x, na.rm = TRUE),
    min = min(original_x, na.rm = TRUE),
    max = max(original_x, na.rm = TRUE),
    iqr = IQR(original_x, na.rm = TRUE),
    q25 = quantile(original_x, 0.25, na.rm = TRUE),
    q75 = quantile(original_x, 0.75, na.rm = TRUE)
  )
  
  # Apply scaling method
  scaled_x <- x
  scaling_params <- list(method = method)
  
  if (method == "none") {
    # No scaling
    scaling_params$center <- 0
    scaling_params$scale <- 1
    
  } else if (method == "center") {
    # Center around mean or custom value
    if (is.null(center_value)) {
      center_val <- if (robust_center) median(x, na.rm = TRUE) else mean(x, na.rm = TRUE)
    } else {
      center_val <- center_value
    }
    scaled_x <- x - center_val
    scaling_params$center <- center_val
    scaling_params$scale <- 1
    
  } else if (method == "standardize") {
    # Z-score standardization
    center_val <- if (robust_center) median(x, na.rm = TRUE) else mean(x, na.rm = TRUE)
    scale_val <- if (robust_center) IQR(x, na.rm = TRUE) else sd(x, na.rm = TRUE)
    
    if (scale_val == 0) {
      warning("Scale value is zero. Using centering only.")
      scaled_x <- x - center_val
      scaling_params$scale <- 1
    } else {
      scaled_x <- (x - center_val) / scale_val
      scaling_params$scale <- scale_val
    }
    scaling_params$center <- center_val
    
  } else if (method == "minmax") {
    # Min-max scaling to [0, 1]
    min_val <- min(x, na.rm = TRUE)
    max_val <- max(x, na.rm = TRUE)
    
    if (min_val == max_val) {
      warning("Min and max values are equal. Using centering only.")
      scaled_x <- x - min_val
    } else {
      scaled_x <- (x - min_val) / (max_val - min_val)
    }
    scaling_params$center <- min_val
    scaling_params$scale <- max_val - min_val
    scaling_params$min <- min_val
    scaling_params$max <- max_val
    
  } else if (method == "iqr") {
    # IQR-based scaling (robust to outliers)
    center_val <- median(x, na.rm = TRUE)
    scale_val <- IQR(x, na.rm = TRUE)
    
    if (scale_val == 0) {
      warning("IQR is zero. Using centering only.")
      scaled_x <- x - center_val
      scaling_params$scale <- 1
    } else {
      scaled_x <- (x - center_val) / scale_val
      scaling_params$scale <- scale_val
    }
    scaling_params$center <- center_val
    
  } else if (method == "robust") {
    # Robust scaling using median and MAD
    center_val <- median(x, na.rm = TRUE)
    scale_val <- mad(x, na.rm = TRUE)
    
    if (scale_val == 0) {
      warning("MAD is zero. Using IQR scaling.")
      scale_val <- IQR(x, na.rm = TRUE)
      if (scale_val == 0) {
        warning("IQR is also zero. Using centering only.")
        scaled_x <- x - center_val
        scaling_params$scale <- 1
      } else {
        scaled_x <- (x - center_val) / scale_val
        scaling_params$scale <- scale_val
      }
    } else {
      scaled_x <- (x - center_val) / scale_val
      scaling_params$scale <- scale_val
    }
    scaling_params$center <- center_val
    
  } else if (method == "quantile") {
    # Quantile-based scaling
    q_low <- quantile(x, quantile_range[1], na.rm = TRUE)
    q_high <- quantile(x, quantile_range[2], na.rm = TRUE)
    
    if (q_low == q_high) {
      warning("Quantile range is zero. Using centering only.")
      scaled_x <- x - q_low
      scaling_params$scale <- 1
    } else {
      scaled_x <- (x - q_low) / (q_high - q_low)
    }
    scaling_params$center <- q_low
    scaling_params$scale <- q_high - q_low
    scaling_params$q_low <- q_low
    scaling_params$q_high <- q_high
    scaling_params$quantile_range <- quantile_range
  }
  
  # Store additional parameters
  scaling_params$original_stats <- original_stats
  scaling_params$na_action <- na_action
  scaling_params$robust_center <- robust_center
  
  # Return results
  result <- list(
    data = scaled_x,
    original = original_x,
    scaling_params = scaling_params,
    method = method
  )
  
  class(result) <- c("scaled_data", "list")
  return(result)
}

#' Reverse Scaling Transformation
#'
#' @description
#' Reverses a scaling transformation to return data to its original scale.
#' This is useful for interpreting results or converting predictions back
#' to the original scale.
#'
#' @param scaled_data Scaled data (numeric vector or scaled_data object)
#' @param scaling_params Scaling parameters (from apply_scaling result)
#'
#' @return Numeric vector with original scale restored
#'
#' @examples
#' # Scale and then reverse
#' data <- c(1, 5, 10, 15, 20, 100)
#' scaled_result <- apply_scaling(data, method = "standardize")
#' reversed_data <- reverse_scaling(scaled_result$data, scaled_result$scaling_params)
#' all.equal(data, reversed_data)  # Should be TRUE
#'
#' @export
reverse_scaling <- function(scaled_data, scaling_params) {
  
  if (inherits(scaled_data, "scaled_data")) {
    # Extract data and params from scaled_data object
    data <- scaled_data$data
    params <- scaled_data$scaling_params
  } else {
    # Use provided data and params
    data <- scaled_data
    params <- scaling_params
  }
  
  method <- params$method
  
  if (method == "none") {
    return(data)
    
  } else if (method == "center") {
    return(data + params$center)
    
  } else if (method %in% c("standardize", "iqr", "robust")) {
    return(data * params$scale + params$center)
    
  } else if (method == "minmax") {
    return(data * params$scale + params$center)
    
  } else if (method == "quantile") {
    return(data * params$scale + params$center)
    
  } else {
    stop("Unknown scaling method: ", method)
  }
}

#' Apply Scaling to Data Frame Columns
#'
#' @description
#' Applies scaling to one or more columns in a data frame while preserving
#' other columns. This is useful for preprocessing multiple time series
#' or features simultaneously.
#'
#' @param data Data frame containing the data
#' @param cols Column names or indices to scale
#' @param method Scaling method to apply
#' @param ... Additional arguments passed to apply_scaling()
#'
#' @return Data frame with scaled columns and scaling information as attributes
#'
#' @examples
#' # Create sample data frame
#' df <- data.frame(
#'   time = 1:100,
#'   ts1 = cumsum(rnorm(100)),
#'   ts2 = cumsum(rnorm(100, sd = 2)),
#'   category = sample(letters[1:3], 100, replace = TRUE)
#' )
#' 
#' # Scale multiple time series columns
#' scaled_df <- scale_dataframe(df, cols = c("ts1", "ts2"), method = "standardize")
#' head(scaled_df)
#'
#' @export
scale_dataframe <- function(data, cols, method = "standardize", ...) {
  
  if (!is.data.frame(data)) {
    stop("Input 'data' must be a data frame")
  }
  
  # Validate columns
  if (is.character(cols)) {
    missing_cols <- cols[!cols %in% names(data)]
    if (length(missing_cols) > 0) {
      stop("Columns not found: ", paste(missing_cols, collapse = ", "))
    }
  } else if (is.numeric(cols)) {
    if (any(cols > ncol(data) | cols < 1)) {
      stop("Column indices out of range")
    }
    cols <- names(data)[cols]
  }
  
  # Apply scaling to each column
  scaled_data <- data
  scaling_info <- list()
  
  for (col in cols) {
    scaling_result <- apply_scaling(data[[col]], method = method, ...)
    scaled_data[[col]] <- scaling_result$data
    scaling_info[[col]] <- scaling_result$scaling_params
  }
  
  # Store scaling information as attributes
  attr(scaled_data, "scaling_info") <- scaling_info
  attr(scaled_data, "scaled_columns") <- cols
  attr(scaled_data, "scaling_method") <- method
  
  class(scaled_data) <- c("scaled_dataframe", class(data))
  return(scaled_data)
}

#' Reverse Scaling for Data Frame
#'
#' @description
#' Reverses scaling transformations applied to a data frame.
#'
#' @param scaled_df Scaled data frame (from scale_dataframe)
#' @param cols Columns to reverse (default: all scaled columns)
#'
#' @return Data frame with original scales restored
#'
#' @examples
#' # Scale and reverse a data frame
#' df <- data.frame(x = 1:10, y = (1:10)^2)
#' scaled_df <- scale_dataframe(df, cols = c("x", "y"), method = "minmax")
#' reversed_df <- reverse_scale_dataframe(scaled_df)
#' all.equal(df, reversed_df[, c("x", "y")])  # Should be TRUE
#'
#' @export
reverse_scale_dataframe <- function(scaled_df, cols = NULL) {
  
  if (!inherits(scaled_df, "scaled_dataframe")) {
    stop("Input must be a scaled_dataframe object")
  }
  
  scaling_info <- attr(scaled_df, "scaling_info")
  scaled_columns <- attr(scaled_df, "scaled_columns")
  
  if (is.null(cols)) {
    cols <- scaled_columns
  }
  
  # Validate columns
  missing_cols <- cols[!cols %in% names(scaling_info)]
  if (length(missing_cols) > 0) {
    stop("No scaling information found for columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Reverse scaling
  result_df <- scaled_df
  
  for (col in cols) {
    result_df[[col]] <- reverse_scaling(scaled_df[[col]], scaling_info[[col]])
  }
  
  # Remove scaling attributes
  attr(result_df, "scaling_info") <- NULL
  attr(result_df, "scaled_columns") <- NULL
  attr(result_df, "scaling_method") <- NULL
  class(result_df) <- class(result_df)[class(result_df) != "scaled_dataframe"]
  
  return(result_df)
}

#' Get Scaling Summary
#'
#' @description
#' Provides a summary of scaling transformations applied to data.
#'
#' @param scaled_object Scaled data object (from apply_scaling or scale_dataframe)
#'
#' @return Summary of scaling information
#'
#' @examples
#' data <- rnorm(100, mean = 50, sd = 10)
#' scaled_result <- apply_scaling(data, method = "standardize")
#' get_scaling_summary(scaled_result)
#'
#' @export
get_scaling_summary <- function(scaled_object) {
  
  if (inherits(scaled_object, "scaled_data")) {
    # Single column scaling
    params <- scaled_object$scaling_params
    
    cat("Scaling Summary\n")
    cat("===============\n")
    cat("Method:", params$method, "\n")
    cat("Center:", round(params$center, 4), "\n")
    cat("Scale:", round(params$scale, 4), "\n")
    
    if (!is.null(params$original_stats)) {
      cat("\nOriginal Data Statistics:\n")
      cat("Mean:", round(params$original_stats$mean, 4), "\n")
      cat("SD:", round(params$original_stats$sd, 4), "\n")
      cat("Range: [", round(params$original_stats$min, 4), ", ", 
          round(params$original_stats$max, 4), "]\n", sep = "")
    }
    
  } else if (inherits(scaled_object, "scaled_dataframe")) {
    # Data frame scaling
    scaling_info <- attr(scaled_object, "scaling_info")
    scaled_cols <- attr(scaled_object, "scaled_columns")
    method <- attr(scaled_object, "scaling_method")
    
    cat("Data Frame Scaling Summary\n")
    cat("==========================\n")
    cat("Method:", method, "\n")
    cat("Scaled columns:", paste(scaled_cols, collapse = ", "), "\n\n")
    
    for (col in scaled_cols) {
      cat("Column:", col, "\n")
      cat("  Center:", round(scaling_info[[col]]$center, 4), "\n")
      cat("  Scale:", round(scaling_info[[col]]$scale, 4), "\n")
      if (!is.null(scaling_info[[col]]$original_stats)) {
        cat("  Original range: [", 
            round(scaling_info[[col]]$original_stats$min, 4), ", ",
            round(scaling_info[[col]]$original_stats$max, 4), "]\n", sep = "")
      }
      cat("\n")
    }
    
  } else {
    stop("Input must be a scaled_data or scaled_dataframe object")
  }
}

#' Print method for scaled_data objects
#' @param x A scaled_data object
#' @param ... Additional arguments (unused)
#' @export
print.scaled_data <- function(x, ...) {
  cat("Scaled Data Object\n")
  cat("==================\n")
  cat("Method:", x$scaling_params$method, "\n")
  cat("Length:", length(x$data), "\n")
  cat("Range: [", round(min(x$data, na.rm = TRUE), 4), ", ", 
      round(max(x$data, na.rm = TRUE), 4), "]\n", sep = "")
  
  if (length(x$data) <= 10) {
    cat("Data:", round(x$data, 4), "\n")
  } else {
    cat("Data (first 10):", round(x$data[1:10], 4), "...\n")
  }
}

#' Print method for scaled_dataframe objects
#' @param x A scaled_dataframe object
#' @param ... Additional arguments (unused)
#' @export
print.scaled_dataframe <- function(x, ...) {
  scaled_cols <- attr(x, "scaled_columns")
  method <- attr(x, "scaling_method")
  
  cat("Scaled Data Frame\n")
  cat("=================\n")
  cat("Dimensions:", nrow(x), "x", ncol(x), "\n")
  cat("Scaling method:", method, "\n")
  cat("Scaled columns:", paste(scaled_cols, collapse = ", "), "\n\n")
  
  # Print the data frame
  NextMethod("print")
}

#' Auto-detect Best Scaling Method
#'
#' @description
#' Automatically selects the most appropriate scaling method based on
#' data characteristics such as distribution, outliers, and range.
#'
#' @param x Numeric vector to analyze
#' @param prefer_robust Whether to prefer robust methods for outlier-heavy data
#'
#' @return Recommended scaling method
#'
#' @examples
#' # Normal-like data
#' normal_data <- rnorm(100)
#' auto_detect_scaling(normal_data)
#' 
#' # Data with outliers
#' outlier_data <- c(rnorm(95), 100, -100)
#' auto_detect_scaling(outlier_data)
#'
#' @export
auto_detect_scaling <- function(x, prefer_robust = TRUE) {
  
  if (!is.numeric(x) || length(x) < 3) {
    return("none")
  }
  
  # Remove NAs for analysis
  x_clean <- x[!is.na(x)]
  
  if (length(x_clean) == 0) {
    return("none")
  }
  
  # Calculate statistics
  mean_x <- mean(x_clean)
  median_x <- median(x_clean)
  sd_x <- sd(x_clean)
  iqr_x <- IQR(x_clean)
  range_x <- diff(range(x_clean))
  
  # Check for constant data
  if (sd_x == 0 || range_x == 0) {
    return("none")
  }
  
  # Check for outliers using IQR method
  q1 <- quantile(x_clean, 0.25)
  q3 <- quantile(x_clean, 0.75)
  outlier_threshold <- 1.5 * iqr_x
  outliers <- sum(x_clean < (q1 - outlier_threshold) | x_clean > (q3 + outlier_threshold))
  outlier_prop <- outliers / length(x_clean)
  
  # Check skewness (simple measure)
  skewness <- abs(mean_x - median_x) / sd_x
  
  # Decision logic
  if (outlier_prop > 0.1 && prefer_robust) {
    if (skewness > 0.5) {
      return("quantile")  # Heavy outliers and skewed
    } else {
      return("iqr")  # Heavy outliers but not too skewed
    }
  } else if (outlier_prop > 0.05 && prefer_robust) {
    return("robust")  # Some outliers
  } else if (range_x > 1000 || min(x_clean) >= 0) {
    return("minmax")  # Large range or all positive
  } else {
    return("standardize")  # Normal case
  }
}

#' Batch Apply Scaling to Multiple Vectors
#'
#' @description
#' Applies the same scaling method to multiple vectors, useful for
#' preprocessing multiple time series with consistent scaling.
#'
#' @param data_list List of numeric vectors to scale
#' @param method Scaling method to apply to all vectors
#' @param ... Additional arguments passed to apply_scaling()
#'
#' @return List of scaled vectors with scaling information
#'
#' @examples
#' # Create multiple time series
#' ts_list <- list(
#'   ts1 = cumsum(rnorm(100)),
#'   ts2 = cumsum(rnorm(100, sd = 2)),
#'   ts3 = cumsum(rnorm(100, mean = 1))
#' )
#' 
#' # Scale all with the same method
#' scaled_list <- batch_scale(ts_list, method = "standardize")
#' lapply(scaled_list, function(x) range(x$data))
#'
#' @export
batch_scale <- function(data_list, method = "standardize", ...) {
  
  if (!is.list(data_list)) {
    stop("Input must be a list of numeric vectors")
  }
  
  # Apply scaling to each vector
  scaled_list <- lapply(data_list, function(x) {
    apply_scaling(x, method = method, ...)
  })
  
  # Add batch information
  attr(scaled_list, "batch_method") <- method
  attr(scaled_list, "batch_size") <- length(data_list)
  
  class(scaled_list) <- c("batch_scaled", "list")
  return(scaled_list)
}
