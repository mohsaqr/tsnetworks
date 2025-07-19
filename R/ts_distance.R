#' Time Series Distance Calculation
#' 
#' Calculate distances between time series windows using various methods
#' 
#' @param ts_data Time series data (vector)
#' @param method Distance method ("euclidean", "manhattan", "cosine", "correlation", "dtw", etc.)
#' @param window_size Size of sliding windows (if NULL, uses full series)
#' @param step_size Step size for windows (defaults to window_size)
#' @param pairwise If TRUE, calculates all pairwise distances; if FALSE, consecutive windows only
#' @param symmetric If TRUE, makes matrix symmetric (both directions); if FALSE, keeps directional
#' @param ... Additional arguments passed to distance calculation
#' 
#' @return Distance matrix (always square matrix, regardless of pairwise setting)
#' @export
ts_distance <- function(ts_data, 
                       method = "euclidean", 
                       window_size = NULL, 
                       step_size = NULL,
                       pairwise = TRUE,
                       symmetric = TRUE,
                       ...) {
  
  # Validate inputs
  if (is.null(ts_data) || length(ts_data) < 2) {
    stop("ts_data must be a vector with at least 2 elements")
  }
  
  # If no windowing, calculate distance between series
  if (is.null(window_size)) {
    stop("For full series comparison, provide two series or use window_size")
  }
  
  # Use the updated analyze_ts_distances function
  result <- analyze_ts_distances(
    ts_data = ts_data,
    method = method,
    window_size = window_size,
    step_size = step_size,
    pairwise = pairwise,
    ...
  )
  
  # Handle different result structures from analyze_ts_distances
  if (pairwise) {
    # For pairwise=TRUE, result should have adjacency_matrix
    if (!is.null(result$adjacency_matrix)) {
      dist_matrix <- result$adjacency_matrix
      
      # Ensure symmetry if requested (should already be symmetric for pairwise=TRUE)
      if (symmetric && !isSymmetric(dist_matrix)) {
        # Make symmetric by taking average of upper and lower triangular parts
        dist_matrix <- (dist_matrix + t(dist_matrix)) / 2
      }
      
      # Add attributes
      attr(dist_matrix, "pairwise") <- pairwise
      attr(dist_matrix, "symmetric") <- symmetric
      attr(dist_matrix, "method") <- method
      attr(dist_matrix, "window_size") <- window_size
      attr(dist_matrix, "step_size") <- step_size
      
      return(dist_matrix)
    } else {
      warning("No adjacency matrix found in pairwise result.")
      return(NULL)
    }
    
  } else {
    # For pairwise=FALSE, result should have edges
    if (is.null(result$edges) || nrow(result$edges) == 0) {
      warning("No distances calculated. Check your data and parameters.")
      return(NULL)
    }
    
    n_windows <- result$n_windows
    
    # Create square distance matrix filled with infinity (indicating no direct connection)
    dist_matrix <- matrix(Inf, nrow = n_windows, ncol = n_windows)
    diag(dist_matrix) <- 0  # Distance from window to itself is 0
    
    # Fill the matrix with calculated distances (consecutive only)
    for (i in 1:nrow(result$edges)) {
      from_idx <- result$edges$from_window[i]
      to_idx <- result$edges$to_window[i]
      distance <- result$edges$distance[i]
      
      dist_matrix[from_idx, to_idx] <- distance
      
      # Make symmetric if requested (most common for network analysis)
      if (symmetric) {
        dist_matrix[to_idx, from_idx] <- distance
      }
    }
    
    # Add row/column names if available
    if (!is.null(result$nodes) && nrow(result$nodes) > 0) {
      rownames(dist_matrix) <- colnames(dist_matrix) <- paste0("W", result$nodes$window_id)
    }
    
    # Add attributes
    attr(dist_matrix, "pairwise") <- pairwise
    attr(dist_matrix, "symmetric") <- symmetric
    attr(dist_matrix, "method") <- method
    attr(dist_matrix, "window_size") <- window_size
    attr(dist_matrix, "step_size") <- step_size
    
    return(dist_matrix)
  }
}

#' Get Available Distance Methods
#' 
#' Returns list of all available distance calculation methods
#' 
#' @return Character vector of method names
#' @export
ts_distance_methods <- function() {
  return(c(
    "euclidean", "manhattan", "cosine", "correlation", 
    "chebyshev", "canberra", "minkowski",
    "dtw", "cross_correlation", "nmi", "voi", "mic"
  ))
}

#' Calculate Distance Between Two Time Series
#'
#' A simple function to compute the distance between two individual time series (vectors).
#' It supports various distance metrics, leveraging both native implementations and
#' methods available through the `proxy` package.
#'
#' @param ts1 A numeric vector representing the first time series.
#' @param ts2 A numeric vector representing the second time series. Must be of the same length as `ts1`.
#' @param method Character string. The distance calculation method to use (e.g., "euclidean",
#'   "manhattan", "cosine", "dtw"). See `ts_distance_methods()` for available options.
#' @param ... Additional arguments passed to the specific distance function
#'   (e.g., `p` for Minkowski distance).
#'
#' @return A single numeric value representing the calculated distance between `ts1` and `ts2`.
#' @importFrom proxy dist
#' @export
ts_distance_pair <- function(ts1, ts2, method = "euclidean", ...) {
  
  if (length(ts1) != length(ts2)) {
    stop("Time series must have the same length")
  }
  
  # Use native methods or proxy
  if (method %in% c("euclidean", "manhattan", "cosine", "correlation", "chebyshev", "canberra", "minkowski")) {
    result <- calculate_native_distance(ts1, ts2, method, ...)
  } else {
    # Use proxy for other methods
    if (!requireNamespace("proxy", quietly = TRUE)) {
      stop("Package 'proxy' needed for method: ", method)
    }
    dist_obj <- proxy::dist(x = list(ts1, ts2), method = method, ...)
    result <- as.matrix(dist_obj)[1, 2]
  }
  
  return(result)
}
