#' @title Discretization Method Registry
#' @description Registry for time series discretization methods
#' @importFrom stats kmeans quantile sd qnorm dist hclust cutree
#' @importFrom mclust Mclust
#' @importFrom dtwclust tsclust
#' @importFrom proxy dist
#' @importFrom zoo rollapply
#' @importFrom utils head tail
#' @noRd
discretization_methods_registry <- new.env(parent = emptyenv())

#' Register Discretization Method
#' @param name Method name
#' @param fn Function implementing the method
#' @noRd
register_method <- function(name, fn) {
    if (!is.function(fn)) stop("Must provide a function")
    discretization_methods_registry[[name]] <- fn
    invisible()
}

#' Get Registered Method
#' @description Retrieves a registered time series discretization method from the internal registry.
#' @param name Character string. The name of the discretization method to retrieve.
#' @return A function implementing the specified discretization method.
#' @keywords internal
#' @noRd
get_method <- function(name) {
    if (!exists(name, discretization_methods_registry)) {
        stop("Method '", name, "' not found. Available: ",
             paste(ls(discretization_methods_registry), collapse=", "))
    }
    discretization_methods_registry[[name]]
}

#' Discretize using K-means clustering
#' @param series Numeric vector
#' @param num_states Number of states
#' @param state_info Optional state information
#' @param suppress_warnings Whether to suppress warnings
#' @return List containing states and centers
#' @noRd
.discretize_kmeans <- function(series, num_states, state_info, suppress_warnings = FALSE, ...) {
    if (length(unique(na.omit(series))) < num_states) {
        if (!suppress_warnings) 
            warning("Fewer unique values than states requested")
        return(NULL)
    }
    
    # Perform kmeans clustering
    kmeans_result <- try({
        stats::kmeans(na.omit(series), 
                     centers = num_states,
                     nstart = 25)
    }, silent = TRUE)
    
    if (inherits(kmeans_result, "try-error")) {
        if (!suppress_warnings) 
            warning("K-means clustering failed")
        return(NULL)
    }
    
    # Assign states to all points including NAs
    states <- rep(NA, length(series))
    non_na <- !is.na(series)
    states[non_na] <- kmeans_result$cluster
    
    list(
        states = states,
        centers = kmeans_result$centers
    )
}

#' Discretize using quantiles
#' @param series Numeric vector
#' @param num_states Number of states
#' @param state_info Optional state information
#' @param suppress_warnings Whether to suppress warnings
#' @return List containing states and breakpoints
#' @noRd
.discretize_quantile <- function(series, num_states, state_info, suppress_warnings = FALSE, ...) {
    if (length(unique(na.omit(series))) < num_states) {
        if (!suppress_warnings) 
            warning("Fewer unique values than states requested")
        return(NULL)
    }
    
    # Calculate breakpoints
    probs <- seq(0, 1, length.out = num_states + 1)
    breaks <- stats::quantile(series, probs = probs, na.rm = TRUE)
    
    # Handle edge cases
    if (any(duplicated(breaks))) {
        if (!suppress_warnings)
            warning("Duplicate breakpoints found in quantile discretization")
        return(NULL)
    }
    
    # Assign states
    states <- cut(series, 
                 breaks = breaks, 
                 labels = FALSE, 
                 include.lowest = TRUE)
    
    list(
        states = states,
        breaks = breaks
    )
}

#' Discretize using change points
#' @param series Numeric vector
#' @param num_states Number of states
#' @param state_info Optional state information
#' @param suppress_warnings Whether to suppress warnings
#' @return List containing states and change points
#' @noRd
.discretize_change_points <- function(series, num_states, state_info, suppress_warnings = FALSE, ...) {
    if (length(series) < num_states * 2) {
        if (!suppress_warnings)
            warning("Series too short for requested number of states")
        return(NULL)
    }
    
    # Calculate differences
    diffs <- diff(series)
    abs_diffs <- abs(diffs)
    
    # Find major change points
    ordered_points <- order(abs_diffs, decreasing = TRUE)
    change_points <- sort(ordered_points[1:(num_states-1)])
    
    # Create states based on segments
    states <- rep(1, length(series))
    for (i in seq_along(change_points)) {
        states[(change_points[i] + 1):length(series)] <- i + 1
    }
    
    list(
        states = states,
        change_points = change_points
    )
}

#' Discretize using equal width bins
#' @param series Numeric vector
#' @param num_states Number of states
#' @param state_info Optional state information
#' @param suppress_warnings Whether to suppress warnings
#' @return List containing states and breaks
#' @noRd
.discretize_equal_width <- function(series, num_states, state_info, suppress_warnings = FALSE, ...) {
    if (length(unique(na.omit(series))) < num_states) {
        if (!suppress_warnings)
            warning("Fewer unique values than states requested")
        return(NULL)
    }
    
    # Calculate breaks
    range <- range(series, na.rm = TRUE)
    width <- diff(range) / num_states
    breaks <- seq(range[1], range[2], length.out = num_states + 1)
    
    # Assign states
    states <- cut(series, 
                 breaks = breaks, 
                 labels = FALSE, 
                 include.lowest = TRUE)
    
    list(
        states = states,
        breaks = breaks
    )
}

#' Discretize using symbolic representation
#' @param series Numeric vector
#' @param num_states Number of states
#' @param state_info Optional state information
#' @param suppress_warnings Whether to suppress warnings
#' @return List containing states and breakpoints
#' @importFrom stats qnorm
#' @noRd
.discretize_symbolic <- function(series, num_states, state_info, suppress_warnings = FALSE, ...) {
    # Z-normalize the series
    mean_val <- mean(series, na.rm = TRUE)
    sd_val <- sd(series, na.rm = TRUE)
    
    if (sd_val == 0) {
        if (!suppress_warnings)
            warning("Zero standard deviation in series")
        return(NULL)
    }
    
    norm_series <- (series - mean_val) / sd_val
    
    # Use Gaussian breakpoints
    breakpoints <- qnorm(seq(0, 1, length.out = num_states + 1))
    
    # Assign states
    states <- cut(norm_series, 
                 breaks = breakpoints, 
                 labels = FALSE, 
                 include.lowest = TRUE)
    
    list(
        states = states,
        breakpoints = breakpoints,
        mean = mean_val,
        sd = sd_val
    )
}

#' Discretize using Fixed Magnitude Breaks
#' @param series Numeric vector
#' @param num_states Number of states
#' @param state_info Optional state information
#' @param suppress_warnings Whether to suppress warnings
#' @return List containing states and breakpoints
#' @noRd
.discretize_magnitude <- function(series, num_states, state_info, suppress_warnings = FALSE, ...) {
    args <- list(...)
    magnitude_scale <- args$magnitude_scale %||% NULL
    
    if (is.null(magnitude_scale)) {
        # Auto-generate magnitude scale based on data range
        data_range <- diff(range(series, na.rm = TRUE))
        magnitude_scale <- exp(seq(log(data_range/num_states), log(data_range), length.out = num_states))
    }
    
    # Create breaks
    breaks <- sort(unique(c(-Inf, magnitude_scale, Inf)))
    
    # Assign states
    states <- cut(abs(series), 
                 breaks = breaks, 
                 labels = FALSE, 
                 include.lowest = TRUE)
    
    list(
        states = states,
        breaks = breaks,
        magnitude_scale = magnitude_scale
    )
}

#' Discretize using Adaptive Magnitude
#' @param series Numeric vector
#' @param num_states Number of states
#' @param state_info Optional state information
#' @param suppress_warnings Whether to suppress warnings
#' @return List containing states and breakpoints
#' @noRd
.discretize_adaptive_magnitude <- function(series, num_states, state_info, suppress_warnings = FALSE, ...) {
    args <- list(...)
    window_size <- args$window_size %||% 10
    
    if (!requireNamespace("zoo", quietly = TRUE)) {
        if (!suppress_warnings) 
            warning("Package 'zoo' required for 'adaptive_magnitude'")
        return(NULL)
    }
    
    # Calculate rolling statistics
    roll_mean <- zoo::rollapply(series, 
                               width = window_size, 
                               FUN = mean, 
                               fill = NA, 
                               align = "right")
    roll_sd <- zoo::rollapply(series, 
                             width = window_size, 
                             FUN = sd, 
                             fill = NA, 
                             align = "right")
    
    # Calculate z-scores
    z_scores <- (series - roll_mean) / roll_sd
    
    # Create breaks based on z-scores
    breaks <- stats::quantile(z_scores, 
                            probs = seq(0, 1, length.out = num_states + 1),
                            na.rm = TRUE)
    
    # Assign states
    states <- cut(z_scores, 
                 breaks = breaks, 
                 labels = FALSE, 
                 include.lowest = TRUE)
    
    list(
        states = states,
        breaks = breaks,
        z_scores = z_scores
    )
}

#' Discretize using Entropy-based Binning
#' @param series Numeric vector
#' @param num_states Number of states
#' @param state_info Optional state information
#' @param suppress_warnings Whether to suppress warnings
#' @return List containing states and breakpoints
#' @noRd
.discretize_entropy <- function(series, num_states, state_info, suppress_warnings = FALSE, ...) {
    if (length(unique(na.omit(series))) < num_states) {
        if (!suppress_warnings)
            warning("Fewer unique values than states requested")
        return(NULL)
    }
    
    # Calculate entropy for different bin widths
    range_data <- range(series, na.rm = TRUE)
    min_width <- diff(range_data) / (num_states * 2)
    max_width <- diff(range_data) / num_states
    
    # Add small padding to range to ensure all points are included
    padding <- diff(range_data) * 0.001
    range_data <- range_data + c(-padding, padding)
    
    widths <- seq(min_width, max_width, length.out = 20)
    entropies <- sapply(widths, function(w) {
        breaks <- seq(range_data[1], range_data[2], by = w)
        # Ensure breaks cover the data
        if (breaks[1] > min(series, na.rm=TRUE)) breaks[1] <- min(series, na.rm=TRUE)
        if (breaks[length(breaks)] < max(series, na.rm=TRUE)) breaks[length(breaks)] <- max(series, na.rm=TRUE)
        hist_result <- hist(series, breaks = breaks, plot = FALSE)
        probs <- hist_result$counts / sum(hist_result$counts)
        probs <- probs[probs > 0]
        -sum(probs * log(probs))
    })
    
    # Use width that maximizes entropy
    optimal_width <- widths[which.max(entropies)]
    breaks <- seq(range_data[1], range_data[2], by = optimal_width)
    # Ensure breaks cover the data
    if (breaks[1] > min(series, na.rm=TRUE)) breaks[1] <- min(series, na.rm=TRUE)
    if (breaks[length(breaks)] < max(series, na.rm=TRUE)) breaks[length(breaks)] <- max(series, na.rm=TRUE)
    
    # Ensure we have the right number of states
    if (length(breaks) - 1 != num_states) {
        breaks <- seq(min(series, na.rm=TRUE), max(series, na.rm=TRUE), length.out = num_states + 1)
    }
    
    # Discretize using the optimal breaks
    states <- cut(series, breaks = breaks, labels = FALSE, include.lowest = TRUE)
    
    list(
        states = states,
        breaks = breaks
    )
}

#' Discretize using Gaussian mixture model
#' @param series Numeric vector
#' @param num_states Number of states
#' @param state_info Optional state information
#' @param suppress_warnings Whether to suppress warnings
#' @return List containing states and breaks
#' @noRd
.discretize_mixture <- function(series, num_states, state_info, suppress_warnings = FALSE, ...) {
    if (!requireNamespace("mclust", quietly = TRUE)) {
        stop("Package 'mclust' is required for mixture-based discretization")
    }
    
    if (length(unique(na.omit(series))) < num_states) {
        if (!suppress_warnings)
            warning("Fewer unique values than states requested")
        return(NULL)
    }
    
    # Try different model types to ensure convergence
    model <- NULL
    model_types <- c("V", "E")
    
    for (model_type in model_types) {
        tryCatch({
            model <- mclust::Mclust(series, G = num_states, modelNames = model_type, verbose = FALSE)
            if (!is.null(model) && length(unique(model$classification)) == num_states) {
                break
            }
        }, error = function(e) {
            if (!suppress_warnings)
                warning(paste("Model", model_type, "failed:", e$message))
        })
    }
    
    if (is.null(model) || length(unique(model$classification)) != num_states) {
        # Fallback to kmeans if mixture model fails
        if (!suppress_warnings)
            warning("Mixture model failed to find requested number of states, falling back to kmeans")
        kmeans_result <- kmeans(series, centers = num_states)
        states <- kmeans_result$cluster
        means <- kmeans_result$centers[order(kmeans_result$centers)]
        breaks <- c(-Inf, (means[-1] + means[-length(means)]) / 2, Inf)
        return(list(
            states = match(kmeans_result$cluster, order(kmeans_result$centers)),
            breaks = breaks
        ))
    }
    
    # Get state assignments and sort them by mean
    states <- model$classification
    means <- model$parameters$mean
    state_order <- order(means)
    state_map <- match(states, state_order)
    
    # Calculate breaks as midpoints between consecutive means
    sorted_means <- sort(means)
    breaks <- c(-Inf, (sorted_means[-1] + sorted_means[-length(sorted_means)]) / 2, Inf)
    
    list(
        states = state_map,
        breaks = breaks,
        model = model  # Include model for diagnostics
    )
}

#' Discretize using Hierarchical Clustering
#' @param series Numeric vector
#' @param num_states Number of states
#' @param state_info Optional state information
#' @param suppress_warnings Whether to suppress warnings
#' @return List containing states and cluster information
#' @noRd
.discretize_hierarchical <- function(series, num_states, state_info, suppress_warnings = FALSE, ...) {
    if (length(unique(na.omit(series))) < num_states) {
        if (!suppress_warnings)
            warning("Fewer unique values than states requested")
        return(NULL)
    }
    
    # Perform hierarchical clustering
    dist_matrix <- stats::dist(series)
    hclust_result <- stats::hclust(dist_matrix, method = "complete")
    
    # Cut tree to get clusters
    clusters <- stats::cutree(hclust_result, k = num_states)
    
    # Calculate cluster centers
    centers <- tapply(series, clusters, mean)
    
    list(
        states = clusters,
        centers = centers,
        hclust = hclust_result
    )
}

#' Discretize using Percentile-based Magnitude
#' @param series Numeric vector
#' @param num_states Number of states
#' @param state_info Optional state information
#' @param suppress_warnings Whether to suppress warnings
#' @return List containing states and breakpoints
#' @noRd
.discretize_percentile_magnitude <- function(series, num_states, state_info, suppress_warnings = FALSE, ...) {
    args <- list(...)
    magnitude_scale <- args$magnitude_scale %||% NULL
    
    if (is.null(magnitude_scale)) {
        # Create percentile breaks
        probs <- seq(0, 1, length.out = num_states + 1)
    } else {
        probs <- magnitude_scale
        if (min(probs) < 0 || max(probs) > 1) {
            if (!suppress_warnings)
                warning("Percentiles should be 0-1. Clamping.")
            probs <- pmax(0, pmin(1, probs))
        }
        probs <- sort(unique(c(0, probs, 1)))
    }
    
    # Calculate breaks using percentiles
    breaks <- stats::quantile(abs(series), probs = probs, na.rm = TRUE)
    
    # Assign states
    states <- cut(abs(series), 
                 breaks = breaks, 
                 labels = FALSE, 
                 include.lowest = TRUE)
    
    list(
        states = states,
        breaks = breaks,
        percentiles = probs
    )
}

#' Discretize using DTW clustering
#' @param series Numeric vector
#' @param num_states Number of states
#' @param state_info Optional state information
#' @param suppress_warnings Whether to suppress warnings
#' @return List containing states and breaks
#' @noRd
.discretize_dtw <- function(series, num_states, state_info = list(), suppress_warnings = FALSE, ...) {
    if (!requireNamespace("dtwclust", quietly = TRUE)) {
        stop("Package 'dtwclust' is required for DTW clustering.")
    }
    if (is.vector(series)) series <- matrix(series, nrow = 1)
    cl <- dtwclust::tsclust(series, type = "partitional", k = num_states, distance = "dtw", centroid = "pam")
    states <- cl@cluster
    list(states = states, breaks = NA)
}


#' Discretize a time series using windowed hierarchical clustering with proxy distances
#' @param series Numeric vector (the time series)
#' @param num_states Number of clusters
#' @param state_info List with window_size, step, distance_metric, cluster_method
#' @param suppress_warnings Logical
#' @param ... Additional arguments
#' @return List with 'states' (cluster assignment for each point) and 'breaks' (NA)
.discretize_proxy_windowed <- function(
  series, num_states, state_info = list(), suppress_warnings = FALSE, ...
) {
  if (!requireNamespace("proxy", quietly = TRUE)) {
    stop("Package 'proxy' is required for proxy_windowed clustering.")
  }
  window_size <- state_info$window_size %||% 10
  step <- state_info$step %||% 1
  distance_method <- state_info$distance_metric %||% "euclidean"
  cluster_method <- state_info$cluster_method %||% "ward.D2"
  
  n <- length(series)
  starts <- seq(1, n - window_size + 1, by = step)
  num_windows <- length(starts)
  if (num_windows < num_states) {
    stop("Number of windows is less than the number of states requested.")
  }
  windows <- sapply(starts, function(i) series[i:(i + window_size - 1)])
  windows <- t(windows)
  
  dist_mat <- proxy::dist(windows, method = distance_method)
  hc <- hclust(dist_mat, method = cluster_method)
  window_clusters <- cutree(hc, k = num_states)
  
  cluster_assignment <- rep(NA, n)
  window_membership <- lapply(1:num_windows, function(i) starts[i]:(starts[i] + window_size - 1))
  point_clusters <- vector("list", n)
  for (i in seq_along(window_membership)) {
    for (idx in window_membership[[i]]) {
      if (idx >= 1 && idx <= n) {
        point_clusters[[idx]] <- c(point_clusters[[idx]], window_clusters[i])
      }
    }
  }
  for (i in seq_along(point_clusters)) {
    if (length(point_clusters[[i]]) > 0) {
      cluster_assignment[i] <- as.integer(names(sort(table(point_clusters[[i]]), decreasing = TRUE)[1]))
    }
  }
  list(states = cluster_assignment, breaks = NA)
}

#' Discretize a pair of time series
#'
#' This internal helper function discretizes two time series simultaneously
#' using a specified method. This is useful for information-theoretic measures
#' that require joint discretization.
#'
#' @param series1 Numeric vector. The first time series.
#' @param series2 Numeric vector. The second time series.
#' @param n_bins Integer. The number of bins to use for discretization.
#' @param method Character string. The discretization method to use (e.g., "quantile", "equal_width").
#' @return A list containing two discretized numeric vectors, `d1` and `d2`.
#' @keywords internal
.discretize_series_pair <- function(series1, series2, n_bins, method = "quantile") {
  # Combine series for joint discretization if method requires it
  combined_series <- c(series1, series2)
  
  # Use a switch to call the appropriate discretization function
  discretized_combined <- switch(method,
    "quantile" = .discretize_quantile(combined_series, num_states = n_bins),
    "equal_width" = .discretize_equal_width(combined_series, num_states = n_bins),
    stop("Unsupported discretization method for series pair: ", method)
  )
  
  if (is.null(discretized_combined)) {
    stop("Discretization failed for series pair.")
  }
  
  # Split back into two series
  d1 <- discretized_combined$states[1:length(series1)]
  d2 <- discretized_combined$states[(length(series1) + 1):(length(series1) + length(series2))]
  
  return(list(d1 = d1, d2 = d2))
}

# Register all methods
register_method("kmeans", .discretize_kmeans)
register_method("quantile", .discretize_quantile)
register_method("change_points", .discretize_change_points)
register_method("symbolic", .discretize_symbolic)
register_method("magnitude", .discretize_magnitude)
register_method("adaptive_magnitude", .discretize_adaptive_magnitude)
register_method("entropy", .discretize_entropy)
register_method("equal_width", .discretize_equal_width)
register_method("hierarchical", .discretize_hierarchical)
register_method("percentile_magnitude", .discretize_percentile_magnitude)
register_method("mixture", .discretize_mixture)
register_method("dtw", .discretize_dtw)
register_method("proxy_windowed", .discretize_proxy_windowed)