# Time Series Distance Calculator with Enhanced Features
# Updated to use analyze_ts_distances with matrix output and native distance functions

#' Enhanced Time Series Distance Analysis with Matrix Output
#' 
#' Unified interface for time series distance calculation with network-ready output
#' 
#' @param ts_data Time series data (vector or data frame)
#' @param method Distance calculation method
#' @param window_size Window size for sliding window analysis
#' @param step_size Step size for windows (defaults to window_size)
#' @param ts_col_name Column name if ts_data is data frame
#' @param id_col_name ID column name for partitioned analysis
#' @param pairwise If TRUE, all pairwise distances; if FALSE, consecutive only
#' @param output_format Output format ("compact" or "expanded")
#' @param ... Additional arguments passed to distance functions
#' 
#' @return List with adjacency_matrix and additional information
#' @export
#' @importFrom stats setNames
#' @importFrom utils capture.output

analyze_ts_distances <- function(ts_data,
                                method,
                                window_size = NULL, step_size = NULL,
                                id_col_name = NULL, ts_col_name = NULL,
                                output_format = "matrix", pairwise = TRUE, ...) {
  
  # Validate inputs
  if (missing(method)) {
    stop("Method must be specified. Use get_all_methods() to see available options.")
  }
  
  if (!output_format %in% c("matrix", "edges")) {
    stop("output_format must be 'matrix' or 'edges'")
  }
  
  # Determine analysis type and extract data
  if (is.data.frame(ts_data) && !is.null(id_col_name)) {
    # ID-based partitioning
    analysis_type <- "id_partition"
    time_series <- extract_ts_from_dataframe(ts_data, ts_col_name, id_col_name)
    result <- perform_id_partition_analysis(time_series, method, output_format, pairwise, ...)
    
  } else if (!is.null(window_size)) {
    # Windowed analysis
    analysis_type <- "windowed"
    time_series <- extract_ts_vector(ts_data, ts_col_name)
    result <- perform_windowed_analysis(time_series, method, window_size, step_size, output_format, pairwise, ...)
    
  } else {
    stop("Specify either 'window_size' for windowed analysis or 'id_col_name' for ID-based partitioning")
  }
  
  # Add metadata
  result$metadata <- list(
    analysis_type = analysis_type,
    method_used = method,
    output_format = output_format,
    pairwise = pairwise,
    n_segments_or_windows = if (analysis_type == "windowed") result$n_windows else result$n_segments,
    timestamp = Sys.time()
  )
  
  return(result)
}

#' Extract time series vector from various input formats
#' @param ts_data Time series data, can be a vector or a data frame.
#' @param ts_col_name The name of the column containing the time series data if ts_data is a data frame.
#' @return A numeric vector of the time series.
extract_ts_vector <- function(ts_data, ts_col_name = NULL) {
  if (is.vector(ts_data) && !is.list(ts_data)) {
    return(ts_data)
  } else if (is.data.frame(ts_data)) {
    if (is.null(ts_col_name)) {
      numeric_cols <- sapply(ts_data, is.numeric)
      if (sum(numeric_cols) == 1) {
        return(ts_data[[which(numeric_cols)]])
      } else {
        stop("Multiple numeric columns found. Please specify ts_col_name.")
      }
    } else {
      if (!ts_col_name %in% names(ts_data)) {
        stop(paste("Column", ts_col_name, "not found in data frame"))
      }
      return(ts_data[[ts_col_name]])
    }
  } else {
    stop("ts_data must be a vector or data frame")
  }
}

#' Extract time series and ID information from data frame
#' @param ts_data The data frame containing the time series and ID data.
#' @param ts_col_name The name of the column containing the time series data.
#' @param id_col_name The name of the column containing the ID data.
#' @return A list containing the time series vector and a list of segments with their IDs, start and end indices.
extract_ts_from_dataframe <- function(ts_data, ts_col_name, id_col_name) {
  validate_dataframe_input(ts_data, ts_col_name, id_col_name)
  
  # Auto-detect time series column if not specified
  if (is.null(ts_col_name)) {
    numeric_cols <- sapply(ts_data, is.numeric)
    non_id_numeric_cols <- numeric_cols & names(ts_data) != id_col_name
    
    if (sum(non_id_numeric_cols) == 1) {
      ts_col_name <- names(ts_data)[non_id_numeric_cols][1]
    } else {
      stop("Cannot auto-detect time series column. Please specify ts_col_name.")
    }
  }
  
  time_series <- ts_data[[ts_col_name]]
  ids <- ts_data[[id_col_name]]
  
  # Create segment information
  unique_ids <- unique(ids)
  segments <- list(
    ids = character(),
    starts = integer(),
    ends = integer()
  )
  
  for (id in unique_ids) {
    id_indices <- which(ids == id)
    start_idx <- min(id_indices)
    end_idx <- max(id_indices)
    
    segments$ids <- c(segments$ids, as.character(id))
    segments$starts <- c(segments$starts, start_idx)
    segments$ends <- c(segments$ends, end_idx)
  }
  
  return(list(
    series = time_series,
    segments = segments
  ))
}

#' Perform windowed analysis
#' @param time_series The time series data.
#' @param method The distance calculation method.
#' @param window_size The size of the sliding windows.
#' @param step_size The step size for the windows.
#' @param output_format The format of the output, either "matrix" or "edges".
#' @param pairwise If TRUE, all pairwise distances are calculated; if FALSE, only consecutive windows are compared.
#' @param ... Additional arguments passed to the distance functions.
#' @return A list containing the results of the windowed analysis.
perform_windowed_analysis <- function(time_series, method, window_size, step_size, output_format, pairwise, ...) {
  time_series <- validate_time_series(time_series, "time_series")
  window_params <- validate_window_params(window_size, step_size, length(time_series))
  
  current_step_size <- if (is.null(step_size)) window_size else step_size
  
  # Generate windows
  windows_result <- .generate_windows(time_series, window_size, current_step_size, drop_last_if_shorter = TRUE)
  windows <- windows_result$windows
  starts <- windows_result$start_indices
  lengths <- windows_result$actual_lengths
  
  if (length(windows) < 2) {
    stop("Need at least 2 windows for analysis. Try smaller window_size or step_size.")
  }
  
  # Calculate distances
  n_windows <- length(windows)
  if (pairwise) {
    # Calculate all pairwise distances
    distance_matrix <- matrix(0, n_windows, n_windows)
    
    for (i in 1:n_windows) {
      for (j in 1:n_windows) {
        if (i != j) {
          distance_matrix[i, j] <- calculate_unified_distance(windows[[i]], windows[[j]], method, ...)
        }
      }
    }
    
    # Create window names
    window_names <- paste0("W", 1:n_windows, "_", starts, ":", starts + lengths - 1)
    rownames(distance_matrix) <- window_names
    colnames(distance_matrix) <- window_names
    
    # Create network output
    result <- list(
      adjacency_matrix = distance_matrix,
      window_info = data.frame(
        window_id = 1:n_windows,
        start_index = starts,
        end_index = starts + lengths - 1,
        window_size = lengths,
        stringsAsFactors = FALSE
      ),
      method = method,
      type = "windowed"
    )
    
  } else {
    # Calculate only consecutive distances
    distances <- numeric(n_windows - 1)
    
    for (i in 1:(n_windows - 1)) {
      distances[i] <- calculate_unified_distance(windows[[i]], windows[[i + 1]], method, ...)
    }
    
    # Create edge list
    result <- list(
      edges = data.frame(
        from_window = 1:(n_windows - 1),
        to_window = 2:n_windows,
        distance = distances,
        from_start = starts[1:(n_windows - 1)],
        from_end = starts[1:(n_windows - 1)] + lengths[1:(n_windows - 1)] - 1,
        to_start = starts[2:n_windows],
        to_end = starts[2:n_windows] + lengths[2:n_windows] - 1,
        method = method,
        stringsAsFactors = FALSE
      ),
      nodes = data.frame(
        window_id = 1:n_windows,
        start_index = starts,
        end_index = starts + lengths - 1,
        window_size = lengths,
        stringsAsFactors = FALSE
      ),
      method = method,
      type = "windowed"
    )
  }
  
  result$n_windows <- n_windows
  result$original_series <- time_series
  
  return(result)
}

#' Perform ID-based partition analysis
#' @param ts_info A list containing the time series and segment information.
#' @param method The distance calculation method.
#' @param output_format The format of the output, either "matrix" or "edges".
#' @param pairwise If TRUE, all pairwise distances are calculated; if FALSE, only consecutive segments are compared.
#' @param ... Additional arguments passed to the distance functions.
#' @return A list containing the results of the ID-based partition analysis.
perform_id_partition_analysis <- function(ts_info, method, output_format, pairwise, ...) {
  time_series <- validate_time_series(ts_info$series, "time_series")
  segments <- ts_info$segments
  
  if (length(segments$ids) < 2) {
    stop("Need at least 2 ID segments for analysis")
  }
  
  # Extract all segments
  segment_data <- list()
  for (i in 1:length(segments$ids)) {
    seg_start <- segments$starts[i]
    seg_end <- segments$ends[i]
    segment_data[[i]] <- time_series[seg_start:seg_end]
  }
  
  n_segments <- length(segments$ids)
  
  if (pairwise) {
    # Calculate all pairwise distances
    distance_matrix <- matrix(0, n_segments, n_segments)
    
    for (i in 1:n_segments) {
      for (j in 1:n_segments) {
        if (i != j) {
          distance_matrix[i, j] <- calculate_unified_distance(segment_data[[i]], segment_data[[j]], method, ...)
        }
      }
    }
    
    # Add names
    rownames(distance_matrix) <- segments$ids
    colnames(distance_matrix) <- segments$ids
    
    result <- list(
      adjacency_matrix = distance_matrix,
      segment_info = data.frame(
        segment_id = segments$ids,
        start_index = segments$starts,
        end_index = segments$ends,
        segment_size = segments$ends - segments$starts + 1,
        stringsAsFactors = FALSE
      ),
      method = method,
      type = "id_partition"
    )
    
  } else {
    # Calculate only consecutive distances
    distances <- numeric(n_segments - 1)
    
    for (i in 1:(n_segments - 1)) {
      distances[i] <- calculate_unified_distance(segment_data[[i]], segment_data[[i + 1]], method, ...)
    }
    
    # Create edge list
    result <- list(
      edges = data.frame(
        from_segment = segments$ids[1:(n_segments - 1)],
        to_segment = segments$ids[2:n_segments],
        distance = distances,
        from_start = segments$starts[1:(n_segments - 1)],
        from_end = segments$ends[1:(n_segments - 1)],
        to_start = segments$starts[2:n_segments],
        to_end = segments$ends[2:n_segments],
        method = method,
        stringsAsFactors = FALSE
      ),
      nodes = data.frame(
        segment_id = segments$ids,
        start_index = segments$starts,
        end_index = segments$ends,
        segment_size = segments$ends - segments$starts + 1,
        stringsAsFactors = FALSE
      ),
      method = method,
      type = "id_partition"
    )
  }
  
  result$n_segments <- n_segments
  result$original_series <- time_series
  
  return(result)
}

#' Unified distance calculation
#' @param x First vector
#' @param y Second vector
#' @param method Distance method
#' @param ... Additional parameters
#' @return Distance value
#' @importFrom proxy dist
#' @importFrom stats median
#' @importFrom utils capture.output
calculate_unified_distance <- function(x, y, method, ...) {
  # Check if it's a native method first
  if (is_native_method(method)) {
    return(calculate_native_distance(x, y, method, ...))
  }
  
  # Check for legacy dist_* methods
  if (method %in% c("dist_cor", "dist_ccf", "dist_dtw", "dist_nmi", "dist_voi", "dist_mic", "dist_es", "dist_vr")) {
    result <- do.call(method, list(series1 = x, series2 = y, ...))
    return(extract_distance_value(result, method))
  }
  
  # Try proxy package as fallback
  if (requireNamespace("proxy", quietly = TRUE)) {
    tryCatch({
      distance_matrix <- proxy::dist(x = list(x, y), method = method, ...)
      return(as.matrix(distance_matrix)[1, 2])
    }, error = function(e) {
      stop(paste("Unknown distance method:", method, "\nUse get_all_methods() to see available options."))
    })
  }
  
  stop(paste("Unknown distance method:", method, "\nUse get_all_methods() to see available options."))
}

#' Get all available distance methods
#' @return Vector of all available method names
#' @importFrom proxy pr_DB
get_all_methods <- function() {
  native_methods <- get_native_methods()
  legacy_methods <- c("dist_cor", "dist_ccf", "dist_dtw", "dist_nmi", "dist_voi", "dist_mic", "dist_es", "dist_vr")
  
  all_methods <- c(native_methods, legacy_methods)
  
  # Add proxy methods if available
  if (requireNamespace("proxy", quietly = TRUE)) {
    tryCatch({
      proxy_methods <- proxy::pr_DB$get_entry_names()
      all_methods <- c(all_methods, proxy_methods)
    }, error = function(e) {
      # Ignore proxy errors
    })
  }
  
  return(sort(unique(all_methods)))
}

#' Extract distance value from legacy method output
#' @param result The result from a legacy distance method.
#' @param method_name The name of the method used to generate the result.
#' @return The extracted distance value.
extract_distance_value <- function(result, method_name) {
  if (is.null(result) || length(result) == 0) return(NA)
  
  if (is.list(result)) {
    # Handle standard distance methods
    distance_val <- switch(method_name,
      "dist_cor" = result$correlation,
      "dist_ccf" = result$ccf_value,
      "dist_dtw" = result$distance,
      "dist_nmi" = result$NMI,
      "dist_voi" = result$VI,
      "dist_mic" = result$MIC,
      "dist_es" = result$C,
      "dist_vr" = result$van_rossum_distance,
      result$distance  # Default
    )
    
    if (is.null(distance_val)) return(NA)
    if (length(distance_val) > 1) return(distance_val[1])
    return(distance_val)
  }
  
  return(NA)
}

#' Calculate Time Series Distance (Deprecated)
#'
#' This function is deprecated. Please use `analyze_ts_distances()` instead,
#' which provides a unified interface for time series distance calculations
#' with enhanced features and network-ready output.
#'
#' @param ... Arguments passed to `analyze_ts_distances()`.
#'
#' @return The result of `analyze_ts_distances()`.
#' @seealso \code{\link{analyze_ts_distances}}
#' @keywords internal
#' @export
calculate_ts_distance <- function(...) {
  warning("calculate_ts_distance() is deprecated. Use analyze_ts_distances() with unified method parameter.")
  
  args <- list(...)
  
  # Convert old parameters
  if (!is.null(args$series1_data)) args$ts_data <- args$series1_data
  if (!is.null(args$series_col_name)) args$ts_col_name <- args$series_col_name
  if (!is.null(args$proxy_method)) args$method <- args$proxy_method
  
  # Remove old parameters
  args$series1_data <- NULL
  args$series2_data <- NULL
  args$series_col_name <- NULL
  args$proxy_method <- NULL
  args$comparison_mode <- NULL
  
  do.call(analyze_ts_distances, args)
}