#' @title Visibility Graph Analysis for Time Series
#' @description Functions for creating and analyzing visibility graphs from time series data
#' @importFrom igraph graph_from_edgelist
#' @importFrom parallel makeCluster stopCluster parLapply clusterExport detectCores
#' @importFrom stats quantile
#' @name visibility-graph-package
NULL

#' Compute quantile-based states for time series
#'
#' This function maps time series values to discrete states based on quantiles.
#' It's useful for creating categorical representations of continuous time series data.
#'
#' @param x Numeric vector of time series values
#' @param num_quantiles Integer. Number of quantile-based states to create. Must be >= 2.
#' @param quantile_probs Numeric vector of quantile probabilities. If provided, 
#'   \code{num_quantiles} is ignored. Values must be between 0 and 1 and in 
#'   ascending order.
#' @param labels Character vector of state labels. If NULL, automatic labels are generated.
#'   Length must match the number of states created.
#' @param na_label Character string to use for NA values. Default is "NA".
#' @param include_lowest Logical. Should the lowest break be included? Passed to \code{cut}.
#' @param right Logical. Should intervals be closed on the right? Passed to \code{cut}.
#'
#' @return Character vector of the same length as \code{x} containing state labels
#'
#' @details
#' The function creates states by dividing the time series into quantile-based bins.
#' For example, with \code{num_quantiles = 4}, the function creates quartile-based
#' states: Q1 (0-25%), Q2 (25-50%), Q3 (50-75%), and Q4 (75-100%).
#'
#' @examples
#' \dontrun{
#' # Generate sample data
#' set.seed(123)
#' ts_data <- rnorm(100, mean = 50, sd = 15)
#'
#' # Create quartile-based states
#' states_4 <- compute_quantile_states(ts_data, num_quantiles = 4)
#' table(states_4)
#'
#' # Create tercile states with custom labels
#' states_3 <- compute_quantile_states(ts_data, num_quantiles = 3, 
#'                                    labels = c("Low", "Medium", "High"))
#' table(states_3)
#'
#' # Use custom quantile probabilities
#' states_custom <- compute_quantile_states(ts_data, 
#'                                         quantile_probs = c(0, 0.1, 0.9, 1.0),
#'                                         labels = c("Very_Low", "Normal", "Very_High"))
#' table(states_custom)
#' }
#'
#' @export
compute_quantile_states <- function(x, num_quantiles = 4, quantile_probs = NULL,
                                   labels = NULL, na_label = "NA",
                                   include_lowest = TRUE, right = FALSE) {
  
  # Input validation
  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector")
  }
  
  if (length(x) == 0) {
    return(character(0))
  }
  
  # Handle case where all values are NA
  if (all(is.na(x))) {
    return(rep(na_label, length(x)))
  }
  
  # Validate quantile_probs if provided
  if (!is.null(quantile_probs)) {
    if (!is.numeric(quantile_probs)) {
      stop("'quantile_probs' must be a numeric vector")
    }
    if (length(quantile_probs) < 2) {
      stop("'quantile_probs' must have at least 2 values")
    }
    if (any(quantile_probs < 0) || any(quantile_probs > 1)) {
      stop("'quantile_probs' values must be between 0 and 1")
    }
    if (is.unsorted(quantile_probs)) {
      stop("'quantile_probs' must be in ascending order")
    }
    if (quantile_probs[1] != 0 || quantile_probs[length(quantile_probs)] != 1) {
      stop("'quantile_probs' must start with 0 and end with 1")
    }
    
    num_states <- length(quantile_probs) - 1
    breaks <- stats::quantile(x, probs = quantile_probs, na.rm = TRUE, type = 7)
    
  } else {
    # Validate num_quantiles
    if (!is.numeric(num_quantiles) || length(num_quantiles) != 1 || 
        num_quantiles != round(num_quantiles) || num_quantiles < 2) {
      stop("'num_quantiles' must be an integer >= 2")
    }
    
    num_states <- as.integer(num_quantiles)
    probs <- seq(0, 1, length.out = num_states + 1)
    breaks <- stats::quantile(x, probs = probs, na.rm = TRUE, type = 7)
  }
  
  # Handle case where breaks are not unique (e.g., all values are the same)
  if (length(unique(breaks)) < 2) {
    warning("All non-NA values are identical. Returning single state for all non-NA values.")
    result <- rep(ifelse(is.null(labels), "State_1", labels[1]), length(x))
    result[is.na(x)] <- na_label
    return(result)
  }
  
  # Create labels if not provided
  if (is.null(labels)) {
    labels <- paste0("Q", seq_len(num_states))
  } else {
    if (!is.character(labels)) {
      stop("'labels' must be a character vector")
    }
    if (length(labels) != num_states) {
      stop("Length of 'labels' must equal the number of states created")
    }
  }
  
  # Validate other parameters
  if (!is.character(na_label) || length(na_label) != 1) {
    stop("'na_label' must be a single character string")
  }
  
  if (!is.logical(include_lowest) || length(include_lowest) != 1) {
    stop("'include_lowest' must be a single logical value")
  }
  
  if (!is.logical(right) || length(right) != 1) {
    stop("'right' must be a single logical value")
  }
  
  # Create states using cut
  states <- cut(x, breaks = breaks, labels = labels, 
                include.lowest = include_lowest, right = right)
  
  # Convert to character and handle NA values
  result <- as.character(states)
  result[is.na(x)] <- na_label
  
  return(result)
}

#' Check natural visibility between two points
#'
#' Determines if two points in a time series have natural visibility according
#' to the natural visibility graph algorithm.
#'
#' @param x Numeric vector representing the time series
#' @param i Integer index of the first point (1-based)
#' @param j Integer index of the second point (1-based)
#' @param penetrable Non-negative integer. Number of points allowed to penetrate 
#'   the visibility line for limited penetrable visibility graphs. Default is 0.
#'
#' @return Logical value indicating whether points i and j have natural visibility
#'
#' @details
#' Natural visibility between two points exists if you can draw a straight line
#' between them that doesn't intersect with any intermediate points (or intersects
#' with at most \code{penetrable} points for LPVG variants).
#'
#' @keywords internal
check_natural_visibility <- function(x, i, j, penetrable = 0) {
  
  # Input validation
  if (!is.numeric(x) || length(x) < 2) {
    stop("'x' must be a numeric vector with at least 2 elements")
  }
  
  if (!is.numeric(i) || !is.numeric(j) || length(i) != 1 || length(j) != 1 ||
      i != round(i) || j != round(j) || i < 1 || j < 1 || 
      i > length(x) || j > length(x)) {
    stop("'i' and 'j' must be valid integer indices within the range of 'x'")
  }
  
  if (!is.numeric(penetrable) || length(penetrable) != 1 || 
      penetrable != round(penetrable) || penetrable < 0) {
    stop("'penetrable' must be a non-negative integer")
  }
  
  # Convert to integers
  i <- as.integer(i)
  j <- as.integer(j)
  penetrable <- as.integer(penetrable)
  
  # Check for NA values at endpoints
  if (is.na(x[i]) || is.na(x[j])) {
    return(FALSE)
  }
  
  # If points are the same, no visibility
  if (i == j) {
    return(FALSE)
  }
  
  # If points are adjacent, they are always visible
  if (abs(i - j) == 1) {
    return(TRUE)
  }
  
  # Determine the range between the two points
  start_idx <- min(i, j)
  end_idx <- max(i, j)
  
  # For limited penetrable visibility graphs, count penetrations
  penetrations <- 0L
  
  # Check visibility between points i and j
  # For all points between i and j, check if they block visibility
  for (k in (start_idx + 1L):(end_idx - 1L)) {
    # Skip NA values in the middle
    if (is.na(x[k])) {
      next
    }
    
    # Calculate the y-value of the line connecting i and j at position k
    # Using the standard formula for natural visibility
    line_height <- x[end_idx] + ((x[start_idx] - x[end_idx]) * 
                                (end_idx - k) / (end_idx - start_idx))
    
    # If any point is at or above the line, it's a penetration
    if (x[k] >= line_height) {
      penetrations <- penetrations + 1L
      # If we've exceeded the allowed penetrations, visibility is blocked
      if (penetrations > penetrable) {
        return(FALSE)
      }
    }
  }
  
  # If penetrations are within allowed limit, return TRUE
  return(TRUE)
}

#' Check horizontal visibility between two points
#'
#' Determines if two points in a time series have horizontal visibility according
#' to the horizontal visibility graph algorithm.
#'
#' @param x Numeric vector representing the time series
#' @param i Integer index of the first point (1-based)
#' @param j Integer index of the second point (1-based)
#' @param penetrable Non-negative integer. Number of points allowed to penetrate 
#'   the visibility line for limited penetrable visibility graphs. Default is 0.
#'
#' @return Logical value indicating whether points i and j have horizontal visibility
#'
#' @details
#' Horizontal visibility between two points exists if all intermediate points
#' are strictly lower than both endpoints (or at most \code{penetrable} points
#' violate this condition for LPVG variants).
#'
#' @keywords internal
check_horizontal_visibility <- function(x, i, j, penetrable = 0) {
  
  # Input validation (same as natural visibility)
  if (!is.numeric(x) || length(x) < 2) {
    stop("'x' must be a numeric vector with at least 2 elements")
  }
  
  if (!is.numeric(i) || !is.numeric(j) || length(i) != 1 || length(j) != 1 ||
      i != round(i) || j != round(j) || i < 1 || j < 1 || 
      i > length(x) || j > length(x)) {
    stop("'i' and 'j' must be valid integer indices within the range of 'x'")
  }
  
  if (!is.numeric(penetrable) || length(penetrable) != 1 || 
      penetrable != round(penetrable) || penetrable < 0) {
    stop("'penetrable' must be a non-negative integer")
  }
  
  # Convert to integers
  i <- as.integer(i)
  j <- as.integer(j)
  penetrable <- as.integer(penetrable)
  
  # Check for NA values at endpoints
  if (is.na(x[i]) || is.na(x[j])) {
    return(FALSE)
  }
  
  # If points are the same, no visibility
  if (i == j) {
    return(FALSE)
  }
  
  # If points are adjacent, they are always visible
  if (abs(i - j) == 1) {
    return(TRUE)
  }
  
  # Determine the range between the two points
  start_idx <- min(i, j)
  end_idx <- max(i, j)
  
  # For limited penetrable visibility graphs, count penetrations
  penetrations <- 0L
  
  # Check if any point between i and j is higher than or equal to
  # either of the endpoints
  for (k in (start_idx + 1L):(end_idx - 1L)) {
    # Skip NA values in the middle
    if (is.na(x[k])) {
      next
    }
    
    if (x[k] >= x[start_idx] || x[k] >= x[end_idx]) {
      penetrations <- penetrations + 1L
      # If we've exceeded the allowed penetrations, visibility is blocked
      if (penetrations > penetrable) {
        return(FALSE)
      }
    }
  }
  
  # If penetrations are within allowed limit, return TRUE
  return(TRUE)
}

#' Map time series values to discrete states
#'
#' This function converts continuous time series values to discrete states
#' using various mapping strategies.
#'
#' @param x Numeric vector of time series values
#' @param state_map State mapping specification. Can be:
#'   \itemize{
#'     \item A named vector mapping ranges to states (e.g., c("0-10"="low", "11-50"="medium"))
#'     \item A function that takes a numeric value and returns a state
#'     \item A numeric vector specifying break points for equal-width bins
#'     \item A single number specifying the number of equal-width bins to create
#'   }
#'
#' @return Character vector of state values corresponding to time series values
#'
#' @details
#' The function supports several mapping strategies:
#' \itemize{
#'   \item \strong{Named vector}: Use range strings like "0-10", "<5", ">100", "=50"
#'   \item \strong{Function}: Custom mapping function
#'   \item \strong{Numeric breaks}: Explicit break points for binning
#'   \item \strong{Number of bins}: Creates equal-width bins
#' }
#'
#' @keywords internal
map_to_state <- function(x, state_map) {
  
  # Input validation
  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector")
  }
  
  if (length(x) == 0) {
    return(character(0))
  }
  
  if (is.null(state_map)) {
    return(as.character(x))
  }
  
  # If state_map is a function, apply it directly
  if (is.function(state_map)) {
    tryCatch({
      return(sapply(x, state_map, USE.NAMES = FALSE))
    }, error = function(e) {
      stop("Error applying state_map function: ", e$message)
    })
  }
  
  # If state_map is a named vector, use it as a lookup table
  if (!is.null(names(state_map))) {
    # Create a mapping function
    map_fn <- function(value) {
      if (is.na(value)) {
        return(NA_character_)
      }
      
      for (i in seq_along(state_map)) {
        range_str <- names(state_map)[i]
        
        tryCatch({
          # Parse range from string (e.g., "0-10" or "<5" or ">100")
          if (grepl("-", range_str)) {
            # Range format: "min-max"
            limits <- as.numeric(strsplit(range_str, "-", fixed = TRUE)[[1]])
            if (length(limits) == 2 && !any(is.na(limits))) {
              if (value >= limits[1] && value <= limits[2]) {
                return(as.character(state_map[i]))
              }
            }
          } else if (grepl("^<", range_str)) {
            # Less than format: "<max"
            max_val <- as.numeric(substring(range_str, 2))
            if (!is.na(max_val) && value < max_val) {
              return(as.character(state_map[i]))
            }
          } else if (grepl("^>", range_str)) {
            # Greater than format: ">min"
            min_val <- as.numeric(substring(range_str, 2))
            if (!is.na(min_val) && value > min_val) {
              return(as.character(state_map[i]))
            }
          } else if (grepl("^=", range_str)) {
            # Exact match format: "=value"
            exact_val <- as.numeric(substring(range_str, 2))
            if (!is.na(exact_val) && value == exact_val) {
              return(as.character(state_map[i]))
            }
          } else {
            # Try to parse as a single number for exact match
            exact_val <- suppressWarnings(as.numeric(range_str))
            if (!is.na(exact_val) && value == exact_val) {
              return(as.character(state_map[i]))
            }
          }
        }, error = function(e) {
          warning("Error parsing range '", range_str, "': ", e$message)
        })
      }
      
      # If no match found, return NA
      return(NA_character_)
    }
    
    return(sapply(x, map_fn, USE.NAMES = FALSE))
  }
  
  # If state_map is numeric breaks, create equal-width bins
  if (is.numeric(state_map) && length(state_map) > 1) {
    if (any(is.na(state_map))) {
      stop("state_map breaks cannot contain NA values")
    }
    if (is.unsorted(state_map)) {
      stop("state_map breaks must be in ascending order")
    }
    
    # Create bins and labels
    bins <- cut(x, breaks = state_map, include.lowest = TRUE)
    return(as.character(bins))
  }
  
  # If state_map is a single number, create that many equal-width bins
  if (is.numeric(state_map) && length(state_map) == 1) {
    if (is.na(state_map) || state_map != round(state_map) || state_map < 1) {
      stop("state_map must be a positive integer when specifying number of bins")
    }
    
    num_bins <- as.integer(state_map)
    if (num_bins == 1) {
      return(rep("Bin_1", length(x)))
    }
    
    if (all(is.na(x))) {
      return(rep(NA_character_, length(x)))
    }
    
    x_range <- range(x, na.rm = TRUE)
    if (x_range[1] == x_range[2]) {
      # All values are the same
      result <- rep("Bin_1", length(x))
      result[is.na(x)] <- NA_character_
      return(result)
    }
    
    breaks <- seq(x_range[1], x_range[2], length.out = num_bins + 1)
    bins <- cut(x, breaks = breaks, include.lowest = TRUE)
    return(as.character(bins))
  }
  
  # Default: return the original values as strings
  return(as.character(x))
}

#' Calculate edge weight based on decay
#'
#' Computes the weight of an edge between two time points, optionally applying
#' temporal decay based on their distance.
#'
#' @param i Integer index of the first point
#' @param j Integer index of the second point
#' @param decay_factor Non-negative numeric value controlling temporal decay.
#'   If 0, all edges have weight 1. Higher values cause faster decay.
#'
#' @return Numeric edge weight between 0 and 1
#'
#' @details
#' The weight is calculated as exp(-decay_factor * distance), where distance
#' is the absolute difference between indices i and j.
#'
#' @keywords internal
calculate_edge_weight <- function(i, j, decay_factor) {
  
  # Input validation
  if (!is.numeric(i) || !is.numeric(j) || length(i) != 1 || length(j) != 1) {
    stop("'i' and 'j' must be single numeric values")
  }
  
  if (!is.numeric(decay_factor) || length(decay_factor) != 1 || decay_factor < 0) {
    stop("'decay_factor' must be a single non-negative numeric value")
  }
  
  if (is.na(i) || is.na(j) || is.na(decay_factor)) {
    return(NA_real_)
  }
  
  if (decay_factor == 0) {
    return(1.0)  # No decay, all edges have weight 1
  }
  
  # Calculate distance between points
  distance <- abs(j - i)
  
  # Apply exponential decay: weight = exp(-decay_factor * distance)
  weight <- exp(-decay_factor * distance)
  
  return(weight)
}

#' Generate edges for a range of points in the time series
#'
#' This function generates visibility edges for a subset of time series points,
#' which is useful for parallel processing.
#'
#' @param x Numeric vector representing the time series
#' @param range_indices Integer vector of indices to process
#' @param method Character string specifying visibility method: "nvg" or "hvg"
#' @param directed Logical indicating whether to create directed edges
#' @param limit Integer specifying maximum temporal distance for visibility
#' @param penetrable Integer specifying number of allowed penetrations
#' @param decay_factor Numeric decay factor for edge weights
#'
#' @return Matrix with columns "from", "to", "weight" representing edges
#'
#' @keywords internal
generate_visibility_edges <- function(x, range_indices, method, directed, limit, 
                                     penetrable, decay_factor) {
  
  # Input validation
  if (!is.numeric(x) || length(x) < 2) {
    stop("'x' must be a numeric vector with at least 2 elements")
  }
  
  if (!is.numeric(range_indices) || any(range_indices < 1) || 
      any(range_indices > length(x))) {
    stop("'range_indices' must contain valid indices for 'x'")
  }
  
  if (!method %in% c("nvg", "hvg")) {
    stop("'method' must be either 'nvg' or 'hvg'")
  }
  
  if (!is.logical(directed) || length(directed) != 1) {
    stop("'directed' must be a single logical value")
  }
  
  # Initialize edge matrix
  edges <- matrix(nrow = 0, ncol = 3)
  colnames(edges) <- c("from", "to", "weight")
  
  n <- length(x)
  
  # Select the appropriate visibility check function
  visibility_fn <- if (method == "nvg") check_natural_visibility else check_horizontal_visibility
  
  for (i in range_indices) {
    # Skip if the current point has NA value
    if (is.na(x[i])) {
      next
    }
    
    # Define the range of j to check, respecting the limit parameter
    j_start <- if (directed) i + 1L else 1L
    j_end <- min(n, if (is.null(limit)) n else i + limit)
    
    # Ensure j_start is valid
    if (j_start > n) {
      next
    }
    
    for (j in j_start:j_end) {
      if (i != j && !is.na(x[j])) {
        if (visibility_fn(x, i, j, penetrable)) {
          # Calculate edge weight based on decay
          weight <- calculate_edge_weight(i, j, decay_factor)
          
          # Add edge with weight
          edges <- rbind(edges, c(i, j, weight))
        }
      }
    }
  }
  
  return(edges)
}

#' Create a State-Based Visibility Graph
#'
#' Converts a time-point based visibility graph into a state-based graph by
#' aggregating connections between time points that belong to the same discrete states.
#' This is an internal helper function used by `visibility_graph()`.
#'
#' @param ts Numeric vector of time series values.
#' @param states Character vector of state labels for each time point. Must be
#'   of the same length as `ts`.
#' @param method Character string. The visibility graph construction method:
#'   "nvg" or "hvg".
#' @param directed Logical. If `TRUE`, the state-based graph will be directed.
#' @param limit Integer or `NULL`. Maximum temporal distance for visibility connections.
#' @param penetrable Integer. The number of intermediate points allowed to penetrate
#'   the visibility line.
#' @param decay_factor Numeric. The temporal decay factor for edge weights.
#' @param weight_agg Character string. The method to aggregate edge weights between
#'   state pairs. Options: "sum" (default), "mean", "max", "count".
#'
#' @return A square numeric matrix representing the state-based adjacency matrix,
#'   where rows and columns correspond to unique states. The values in the matrix
#'   represent the aggregated weights of connections between states.
#' @keywords internal
create_state_graph <- function(ts, states, method, directed, limit, penetrable, 
                              decay_factor, weight_agg = "sum") {
  
  # Input validation
  if (!is.numeric(ts)) {
    stop("'ts' must be a numeric vector")
  }
  
  if (!is.character(states) && !is.factor(states)) {
    stop("'states' must be a character vector or factor")
  }
  
  if (length(ts) != length(states)) {
    stop("'ts' and 'states' must have the same length")
  }
  
  if (!weight_agg %in% c("sum", "mean", "max", "count")) {
    stop("'weight_agg' must be one of: 'sum', 'mean', 'max', 'count'")
  }
  
  # Convert states to character if factor
  states <- as.character(states)
  
  # Get unique states, excluding NA
  unique_states <- unique(states[!is.na(states)])
  
  if (length(unique_states) == 0) {
    warning("No valid states found. All states are NA.")
    return(matrix(0, nrow = 0, ncol = 0))
  }
  
  n_states <- length(unique_states)
  
  # Create an empty adjacency matrix for states
  state_matrix <- matrix(0, nrow = n_states, ncol = n_states)
  rownames(state_matrix) <- colnames(state_matrix) <- unique_states
  
  # For mean aggregation, we need to keep track of counts
  if (weight_agg == "mean") {
    count_matrix <- matrix(0, nrow = n_states, ncol = n_states)
    rownames(count_matrix) <- colnames(count_matrix) <- unique_states
  }
  
  # Generate visibility edges between time points
  time_edges <- generate_visibility_edges(ts, seq_along(ts), method, directed, 
                                         limit, penetrable, decay_factor)
  
  # Convert time point edges to state edges
  if (nrow(time_edges) > 0) {
    for (i in seq_len(nrow(time_edges))) {
      from_idx <- time_edges[i, 1]
      to_idx <- time_edges[i, 2]
      weight <- time_edges[i, 3]
      
      # Get states for these time points
      from_state <- states[from_idx]
      to_state <- states[to_idx]
      
      # Skip if either state is NA
      if (is.na(from_state) || is.na(to_state)) {
        next
      }
      
      # Add weight to the corresponding state-to-state edge
      if (weight_agg == "sum") {
        state_matrix[from_state, to_state] <- state_matrix[from_state, to_state] + weight
      } else if (weight_agg == "max") {
        state_matrix[from_state, to_state] <- max(state_matrix[from_state, to_state], weight)
      } else if (weight_agg == "mean") {
        state_matrix[from_state, to_state] <- state_matrix[from_state, to_state] + weight
        count_matrix[from_state, to_state] <- count_matrix[from_state, to_state] + 1
      } else if (weight_agg == "count") {
        state_matrix[from_state, to_state] <- state_matrix[from_state, to_state] + 1
      }
      
      # For undirected graphs, mirror the edge
      if (!directed) {
        if (weight_agg == "sum") {
          state_matrix[to_state, from_state] <- state_matrix[to_state, from_state] + weight
        } else if (weight_agg == "max") {
          state_matrix[to_state, from_state] <- max(state_matrix[to_state, from_state], weight)
        } else if (weight_agg == "mean") {
          state_matrix[to_state, from_state] <- state_matrix[to_state, from_state] + weight
          count_matrix[to_state, from_state] <- count_matrix[to_state, from_state] + 1
        } else if (weight_agg == "count") {
          state_matrix[to_state, from_state] <- state_matrix[to_state, from_state] + 1
        }
      }
    }
  }
  
  # For mean aggregation, divide by counts
  if (weight_agg == "mean") {
    # Avoid division by zero
    count_matrix[count_matrix == 0] <- 1
    state_matrix <- state_matrix / count_matrix
  }
  
  return(state_matrix)
}

#' Compute visibility graph from time series
#'
#' Creates visibility graphs from time series data, where nodes represent time points
#' or states, and edges represent visibility between points according to the selected method.
#' The function supports both Natural Visibility Graphs (NVG) and Horizontal Visibility 
#' Graphs (HVG), with options for limited penetrable visibility, temporal decay, and 
#' state-based aggregation.
#'
#' @param data A data frame or numeric vector containing the time series data.
#'   If a data frame, the \code{ts_value} parameter must specify which column
#'   contains the time series values.
#' @param ts_value Character string or NULL. If \code{data} is a data frame, this
#'   specifies the column name containing time series values. If \code{data} is 
#'   a numeric vector, this parameter is ignored.
#' @param state_value Character string, vector, or NULL. If \code{data} is a data frame,
#'   this should be a column name containing state labels. If \code{data} is a vector,
#'   this should be a vector of the same length containing state labels. If provided,
#'   the output will include state-based graph representations.
#' @param method Character string specifying the visibility graph construction method:
#'   \itemize{
#'     \item \code{"nvg"} (default): Natural Visibility Graph
#'     \item \code{"hvg"}: Horizontal Visibility Graph
#'   }
#' @param directed Logical value indicating whether to create a directed graph.
#'   If \code{TRUE}, edges have direction showing temporal order. Default is \code{FALSE}.
#' @param limit Positive integer or NULL. Maximum temporal distance (in time steps)
#'   for visibility connections. If NULL (default), no limit is applied.
#' @param penetrable Non-negative integer. Number of points allowed to penetrate
#'   the visibility line for Limited Penetrable Visibility Graphs (LPVG).
#'   Default is 0 (standard visibility graph).
#' @param decay_factor Non-negative numeric value controlling temporal decay of
#'   edge weights. Default is 0 (no decay). Higher values cause faster decay
#'   with increasing temporal distance.
#' @param state_map State mapping specification used when \code{state_value} is NULL.
#'   Can be:
#'   \itemize{
#'     \item NULL (default): No state mapping
#'     \item Named vector: Maps value ranges to states (e.g., c("0-25"="low", ">50"="high"))
#'     \item Function: Takes numeric values and returns state labels
#'     \item Numeric vector: Break points for equal-width bins
#'     \item Single number: Number of equal-width bins to create
#'   }
#' @param weight_agg Character string specifying how to aggregate edge weights
#'   when creating state-based graphs:
#'   \itemize{
#'     \item \code{"sum"} (default): Sum all weights between state pairs
#'     \item \code{"mean"}: Average of weights between state pairs
#'     \item \code{"max"}: Maximum weight between state pairs
#'     \item \code{"count"}: Count of connections between state pairs
#'   }
#' @param output_format Character string specifying the output format:
#'   \code{"matrix"} (default) or \code{"edgelist"}. This affects the primary
#'   output slots but all representations are always included.
#' @param default_view Character string specifying which graph view to use
#'   for the primary output slots:
#'   \itemize{
#'     \item \code{"time"}: Time-point based graph (default when no states)
#'     \item \code{"state"}: State-based graph (default when states provided)
#'     \item \code{"value"}: Value-based graph (groups by unique values)
#'   }
#' @param num_cores Positive integer specifying number of CPU cores to use
#'   for parallel processing. Default is 1 (no parallelization). Parallelization
#'   is only used for time series with more than 100 points.
#'
#' @return A list containing multiple graph representations:
#'   \describe{
#'     \item{\code{matrix}}{Primary adjacency matrix (format determined by \code{default_view})}
#'     \item{\code{edge_list}}{Primary edge list (format determined by \code{default_view})}
#'     \item{\code{time_matrix}}{Time-point based adjacency matrix}
#'     \item{\code{time_edge_list}}{Time-point based edge list}
#'     \item{\code{values_matrix}}{Value-based adjacency matrix (groups by unique values)}
#'     \item{\code{values_edge_list}}{Value-based edge list}
#'     \item{\code{state_matrix}}{State-based adjacency matrix (if states provided)}
#'     \item{\code{state_edge_list}}{State-based edge list (if states provided)}
#'     \item{\code{original_values}}{Original time series values}
#'     \item{\code{states}}{State labels (if provided)}
#'     \item{\code{method}}{Visibility method used}
#'     \item{\code{directed}}{Whether the graph is directed}
#'     \item{\code{parameters}}{List of all parameters used}
#'   }
#'
#' @details
#' \subsection{Visibility Graph Methods}{
#' \strong{Natural Visibility Graph (NVG):} Two points can "see" each other if
#' a straight line connecting them doesn't intersect with any intermediate points.
#' This preserves more geometric information from the original time series.
#' 
#' \strong{Horizontal Visibility Graph (HVG):} Two points can "see" each other if
#' all intermediate points are strictly lower than both endpoints. This creates
#' sparser graphs that capture different temporal patterns.
#' }
#' 
#' \subsection{Advanced Features}{
#' \strong{Limited Penetrable Visibility (LPVG):} Allows a specified number of
#' points to penetrate the visibility line, creating more connected graphs.
#' 
#' \strong{Temporal Decay:} Edge weights decrease exponentially with temporal
#' distance, emphasizing local temporal relationships.
#' 
#' \strong{State-based Analysis:} Converts continuous time series to discrete
#' states for categorical analysis of temporal patterns.
#' }
#'
#' @section Performance Considerations:
#' \itemize{
#'   \item Time complexity is O(n²) for the basic algorithm
#'   \item Memory usage scales with the number of edges created
#'   \item Parallel processing is automatically used for n > 100 when num_cores > 1
#'   \item State-based graphs can significantly reduce complexity for large time series
#' }
#'
#' @references
#' \itemize{
#'   \item Lacasa, L., Luque, B., Ballesteros, F., Luque, J., & Nuño, J. C. (2008).
#'     From time series to complex networks: The visibility graph. PNAS, 105(13), 4972-4975.
#'   \item Luque, B., Lacasa, L., Ballesteros, F., & Luque, J. (2009).
#'     Horizontal visibility graphs: Exact results for random time series. Physical Review E, 80(4), 046103.
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with a simple time series
#' set.seed(123)
#' ts_data <- rnorm(50, mean = 50, sd = 15)
#' 
#' # Create natural visibility graph
#' nvg <- visibility_graph(ts_data, method = "nvg")
#' print(dim(nvg$matrix))  # Time-based adjacency matrix
#' 
#' # Create horizontal visibility graph with temporal limit
#' hvg <- visibility_graph(ts_data, method = "hvg", limit = 10)
#' 
#' # Using data frame input
#' df <- data.frame(
#'   time = 1:20,
#'   value = c(5, 15, 25, 75, 45, 10, 30, 20, 60, 35,
#'             40, 80, 55, 25, 70, 45, 35, 65, 50, 30),
#'   category = rep(c("A", "B"), each = 10)
#' )
#' 
#' # Create state-based visibility graph
#' state_vg <- visibility_graph(df, ts_value = "value", state_value = "category")
#' print(state_vg$state_matrix)  # State-based adjacency matrix
#' 
#' # Using quantile-based states
#' quantile_states <- compute_quantile_states(ts_data, num_quantiles = 4)
#' qvg <- visibility_graph(ts_data, state_value = quantile_states)
#' 
#' # Advanced features
#' lpvg <- visibility_graph(ts_data, penetrable = 2, decay_factor = 0.1)
#' 
#' # Custom state mapping
#' custom_map <- c("0-25" = "low", "25-75" = "medium", "75-100" = "high")
#' custom_vg <- visibility_graph(ts_data, state_map = custom_map)
#' 
#' # Parallel processing for large time series
#' large_ts <- rnorm(1000)
#' large_vg <- visibility_graph(large_ts, num_cores = 4)
#' }
#'
#' @seealso
#' \code{\link{compute_quantile_states}} for creating quantile-based states
#' 
#' @export
visibility_graph <- function(data, ts_value = NULL, state_value = NULL,
                           method = "nvg", directed = FALSE, limit = NULL,
                           penetrable = 0, decay_factor = 0, state_map = NULL,
                           weight_agg = "sum", output_format = "matrix", 
                           default_view = NULL, num_cores = 1) {
  
  # Store function call for debugging
  call_info <- match.call()
  
  # Comprehensive input validation
  tryCatch({
    
    # Validate weight_agg parameter
    if (!is.character(weight_agg) || length(weight_agg) != 1 || 
        !weight_agg %in% c("sum", "mean", "max", "count")) {
      stop("'weight_agg' must be one of: 'sum', 'mean', 'max', or 'count'")
    }
    
    # Extract time series values
    if (is.data.frame(data)) {
      if (is.null(ts_value)) {
        stop("When 'data' is a data frame, 'ts_value' must be specified")
      }
      if (!is.character(ts_value) || length(ts_value) != 1) {
        stop("'ts_value' must be a single character string")
      }
      if (!ts_value %in% names(data)) {
        stop("Column '", ts_value, "' not found in data frame")
      }
      x <- data[[ts_value]]
    } else if (is.numeric(data)) {
      x <- data
    } else {
      stop("'data' must be either a data frame or a numeric vector")
    }
    
    # Validate time series values
    if (!is.numeric(x)) {
      stop("Time series values must be numeric")
    }
    
    if (length(x) == 0) {
      stop("Time series cannot be empty")
    }
    
    # Extract or generate states if needed
    states <- NULL
    if (!is.null(state_value)) {
      if (is.data.frame(data)) {
        if (!is.character(state_value) || length(state_value) != 1) {
          stop("'state_value' must be a single character string when data is a data frame")
        }
        if (!state_value %in% names(data)) {
          stop("Column '", state_value, "' not found in data frame")
        }
        states <- data[[state_value]]
      } else {
        if (length(state_value) != length(x)) {
          stop("Length of 'state_value' must match length of time series")
        }
        states <- state_value
      }
      
      # Convert states to character
      states <- as.character(states)
      
    } else if (!is.null(state_map)) {
      # Generate states from state_map
      states <- map_to_state(x, state_map)
    }
    
    # Determine default view if not specified
    if (is.null(default_view)) {
      default_view <- if (!is.null(states)) "state" else "time"
    } else {
      if (!is.character(default_view) || length(default_view) != 1 ||
          !default_view %in% c("time", "state", "value")) {
        stop("'default_view' must be one of: 'time', 'state', or 'value'")
      }
    }
    
    # Validate other parameters
    if (!is.character(method) || length(method) != 1 || 
        !method %in% c("nvg", "hvg")) {
      stop("'method' must be either 'nvg' or 'hvg'")
    }
    
    if (!is.logical(directed) || length(directed) != 1) {
      stop("'directed' must be a single logical value")
    }
    
    if (!is.null(limit)) {
      if (!is.numeric(limit) || length(limit) != 1 || 
          limit != round(limit) || limit <= 0) {
        stop("'limit' must be a positive integer")
      }
      limit <- as.integer(limit)
    }
    
    if (!is.numeric(penetrable) || length(penetrable) != 1 || 
        penetrable != round(penetrable) || penetrable < 0) {
      stop("'penetrable' must be a non-negative integer")
    }
    penetrable <- as.integer(penetrable)
    
    if (!is.numeric(decay_factor) || length(decay_factor) != 1 || decay_factor < 0) {
      stop("'decay_factor' must be a non-negative number")
    }
    
    if (!is.character(output_format) || length(output_format) != 1 ||
        !output_format %in% c("matrix", "edgelist")) {
      stop("'output_format' must be either 'matrix' or 'edgelist'")
    }
    
    if (!is.numeric(num_cores) || length(num_cores) != 1 || 
        num_cores != round(num_cores) || num_cores <= 0) {
      stop("'num_cores' must be a positive integer")
    }
    num_cores <- as.integer(num_cores)
    
    # Check for parallel package availability when needed
    if (num_cores > 1) {
      if (!requireNamespace("parallel", quietly = TRUE)) {
        warning("Package 'parallel' not available. Using single core processing.")
        num_cores <- 1L
      } else {
        # Limit to available cores
        max_cores <- parallel::detectCores()
        if (is.na(max_cores)) max_cores <- 1L
        num_cores <- min(num_cores, max_cores)
      }
    }
    
  }, error = function(e) {
    stop("Input validation failed: ", e$message, call. = FALSE)
  })
  
  n <- length(x)
  
  # Handle edge cases
  if (n < 2) {
    warning("Time series has fewer than 2 points. Creating minimal graph structure.")
    
    # Create minimal result structure
    result <- list(
      matrix = matrix(0, nrow = max(n, 0), ncol = max(n, 0)),
      edge_list = matrix(nrow = 0, ncol = 3, 
                        dimnames = list(NULL, c("from", "to", "weight"))),
      time_matrix = matrix(0, nrow = n, ncol = n),
      time_edge_list = matrix(nrow = 0, ncol = 3,
                             dimnames = list(NULL, c("from", "to", "weight"))),
      values_matrix = matrix(0, nrow = 0, ncol = 0),
      values_edge_list = matrix(nrow = 0, ncol = 3,
                               dimnames = list(NULL, c("from", "to", "weight"))),
      original_values = x,
      method = method,
      directed = directed,
      parameters = list(
        method = method, directed = directed, limit = limit,
        penetrable = penetrable, decay_factor = decay_factor,
        weight_agg = weight_agg, num_cores = num_cores
      )
    )
    
    if (!is.null(states)) {
      result$state_matrix <- matrix(0, nrow = 0, ncol = 0)
      result$state_edge_list <- matrix(nrow = 0, ncol = 3,
                                      dimnames = list(NULL, c("from", "to", "weight")))
      result$states <- states
    }
    
    return(result)
  }
  
  # Check for all NA values
  if (all(is.na(x))) {
    warning("All time series values are NA. Creating empty graph structure.")
    
    result <- list(
      matrix = matrix(0, nrow = n, ncol = n),
      edge_list = matrix(nrow = 0, ncol = 3,
                        dimnames = list(NULL, c("from", "to", "weight"))),
      time_matrix = matrix(0, nrow = n, ncol = n),
      time_edge_list = matrix(nrow = 0, ncol = 3,
                             dimnames = list(NULL, c("from", "to", "weight"))),
      values_matrix = matrix(0, nrow = 0, ncol = 0),
      values_edge_list = matrix(nrow = 0, ncol = 3,
                               dimnames = list(NULL, c("from", "to", "weight"))),
      original_values = x,
      method = method,
      directed = directed,
      parameters = list(
        method = method, directed = directed, limit = limit,
        penetrable = penetrable, decay_factor = decay_factor,
        weight_agg = weight_agg, num_cores = num_cores
      )
    )
    
    if (!is.null(states)) {
      result$state_matrix <- matrix(0, nrow = 0, ncol = 0)
      result$state_edge_list <- matrix(nrow = 0, ncol = 3,
                                      dimnames = list(NULL, c("from", "to", "weight")))
      result$states <- states
    }
    
    return(result)
  }
  
  # Generate the time-point based visibility graph
  tryCatch({
    
    # Use parallelization for larger time series if requested
    if (num_cores > 1 && n > 100) {
      
      # Set up parallel cluster
      cl <- parallel::makeCluster(num_cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      
      # Split the indices across cores
      indices_per_core <- split(seq_len(n), cut(seq_len(n), num_cores, labels = FALSE))
      
      # Export necessary functions and variables to the cluster
      parallel::clusterExport(cl, c("check_natural_visibility", "check_horizontal_visibility", 
                                   "calculate_edge_weight", "generate_visibility_edges", 
                                   "x", "method", "directed", "limit", "penetrable", "decay_factor"), 
                             envir = environment())
      
      # Process each chunk in parallel
      edge_lists <- parallel::parLapply(cl, indices_per_core, function(indices) {
        generate_visibility_edges(x, indices, method, directed, limit, penetrable, decay_factor)
      })
      
      # Combine all edges
      time_edges <- do.call(rbind, edge_lists)
      
    } else {
      # Single-core processing
      time_edges <- generate_visibility_edges(x, seq_len(n), method, directed, 
                                             limit, penetrable, decay_factor)
    }
    
  }, error = function(e) {
    stop("Error in edge generation: ", e$message, call. = FALSE)
  })
  
  # Convert to time-based matrix
  time_matrix <- matrix(0, nrow = n, ncol = n)
  
  if (!is.null(time_edges) && nrow(time_edges) > 0) {
    # Fill the adjacency matrix with edge weights
    for (i in seq_len(nrow(time_edges))) {
      from <- time_edges[i, 1]
      to <- time_edges[i, 2]
      weight <- time_edges[i, 3]
      
      # Ensure indices are within bounds and not NA
      if (!is.na(from) && !is.na(to) && !is.na(weight) &&
          from > 0 && from <= n && to > 0 && to <= n) {
        time_matrix[from, to] <- weight
        
        # For undirected graphs, mirror the edge
        if (!directed) {
          time_matrix[to, from] <- weight
        }
      }
    }
  }
  
  # Initialize result list
  result <- list(
    time_matrix = time_matrix,
    time_edge_list = if (is.null(time_edges)) {
      matrix(nrow = 0, ncol = 3, dimnames = list(NULL, c("from", "to", "weight")))
    } else {
      time_edges
    },
    original_values = x,
    method = method,
    directed = directed,
    parameters = list(
      method = method, directed = directed, limit = limit,
      penetrable = penetrable, decay_factor = decay_factor,
      weight_agg = weight_agg, num_cores = num_cores
    )
  )
  
  # Add state-based graph if states are provided
  if (!is.null(states)) {
    tryCatch({
      state_matrix <- create_state_graph(x, states, method, directed, limit, 
                                        penetrable, decay_factor, weight_agg)
      
      # Convert state matrix to edge list
      state_edge_list <- matrix(nrow = 0, ncol = 3)
      colnames(state_edge_list) <- c("from", "to", "weight")
      
      if (nrow(state_matrix) > 0 && ncol(state_matrix) > 0 && 
          !is.null(rownames(state_matrix)) && !is.null(colnames(state_matrix))) {
        for (i in seq_len(nrow(state_matrix))) {
          for (j in seq_len(ncol(state_matrix))) {
            if (state_matrix[i, j] > 0) {
              state_edge_list <- rbind(state_edge_list, 
                                     c(rownames(state_matrix)[i], 
                                       colnames(state_matrix)[j], 
                                       state_matrix[i, j]))
            }
          }
        }
      }
      
      result$state_matrix <- state_matrix
      result$state_edge_list <- state_edge_list
      result$states <- states
      
    }, error = function(e) {
      warning("Error creating state-based graph: ", e$message)
      result$state_matrix <- matrix(0, nrow = 0, ncol = 0)
      result$state_edge_list <- matrix(nrow = 0, ncol = 3,
                                      dimnames = list(NULL, c("from", "to", "weight")))
      result$states <- states
    })
  }
  
  # Create value-based graph (treating each unique value as a node)
  tryCatch({
    value_map <- as.character(x)
    unique_values <- unique(value_map[!is.na(value_map)])
    
    if (length(unique_values) > 0) {
      # Create value matrix
      value_matrix <- matrix(0, nrow = length(unique_values), ncol = length(unique_values))
      rownames(value_matrix) <- colnames(value_matrix) <- unique_values
      
      # For mean aggregation, we need to keep track of counts
      if (weight_agg == "mean") {
        count_matrix <- matrix(0, nrow = length(unique_values), ncol = length(unique_values))
        rownames(count_matrix) <- colnames(count_matrix) <- unique_values
      }
      
      # Convert time point edges to value edges
      if (!is.null(time_edges) && nrow(time_edges) > 0) {
        for (i in seq_len(nrow(time_edges))) {
          from_idx <- time_edges[i, 1]
          to_idx <- time_edges[i, 2]
          weight <- time_edges[i, 3]
          
          # Check bounds
          if (from_idx >= 1 && from_idx <= length(value_map) &&
              to_idx >= 1 && to_idx <= length(value_map)) {
            
            # Get values for these time points
            from_value <- value_map[from_idx]
            to_value <- value_map[to_idx]
            
            # Skip if either value is NA
            if (!is.na(from_value) && !is.na(to_value)) {
              
              # Add weight to the corresponding value-to-value edge
              if (weight_agg == "sum") {
                value_matrix[from_value, to_value] <- value_matrix[from_value, to_value] + weight
              } else if (weight_agg == "max") {
                value_matrix[from_value, to_value] <- max(value_matrix[from_value, to_value], weight)
              } else if (weight_agg == "mean") {
                value_matrix[from_value, to_value] <- value_matrix[from_value, to_value] + weight
                count_matrix[from_value, to_value] <- count_matrix[from_value, to_value] + 1
              } else if (weight_agg == "count") {
                value_matrix[from_value, to_value] <- value_matrix[from_value, to_value] + 1
              }
              
              # For undirected graphs, mirror the edge
              if (!directed) {
                if (weight_agg == "sum") {
                  value_matrix[to_value, from_value] <- value_matrix[to_value, from_value] + weight
                } else if (weight_agg == "max") {
                  value_matrix[to_value, from_value] <- max(value_matrix[to_value, from_value], weight)
                } else if (weight_agg == "mean") {
                  value_matrix[to_value, from_value] <- value_matrix[to_value, from_value] + weight
                  count_matrix[to_value, from_value] <- count_matrix[to_value, from_value] + 1
                } else if (weight_agg == "count") {
                  value_matrix[to_value, from_value] <- value_matrix[to_value, from_value] + 1
                }
              }
            }
          }
        }
      }
      
      # For mean aggregation, divide by counts
      if (weight_agg == "mean") {
        # Avoid division by zero
        count_matrix[count_matrix == 0] <- 1
        value_matrix <- value_matrix / count_matrix
      }
      
    } else {
      value_matrix <- matrix(0, nrow = 0, ncol = 0)
    }
    
    # Convert value matrix to edge list
    value_edge_list <- matrix(nrow = 0, ncol = 3)
    colnames(value_edge_list) <- c("from", "to", "weight")
    
    if (nrow(value_matrix) > 0 && !is.null(rownames(value_matrix)) && !is.null(colnames(value_matrix))) {
      for (i in seq_len(nrow(value_matrix))) {
        for (j in seq_len(ncol(value_matrix))) {
          if (value_matrix[i, j] > 0) {
            value_edge_list <- rbind(value_edge_list, 
                                   c(rownames(value_matrix)[i], 
                                     colnames(value_matrix)[j], 
                                     value_matrix[i, j]))
          }
        }
      }
    }
    
    result$values_matrix <- value_matrix
    result$values_edge_list <- value_edge_list
    
  }, error = function(e) {
    warning("Error creating value-based graph: ", e$message)
    result$values_matrix <- matrix(0, nrow = 0, ncol = 0)
    result$values_edge_list <- matrix(nrow = 0, ncol = 3,
                                     dimnames = list(NULL, c("from", "to", "weight")))
  })
  
  # Set default view in the main output slots
  tryCatch({
    if (default_view == "time") {
      result$matrix <- result$time_matrix
      result$edge_list <- result$time_edge_list
    } else if (default_view == "state" && !is.null(states) && "state_matrix" %in% names(result)) {
      result$matrix <- result$state_matrix
      result$edge_list <- result$state_edge_list
    } else if (default_view == "value") {
      result$matrix <- result$values_matrix
      result$edge_list <- result$values_edge_list
    } else {
      # Fallback to time if requested view isn't available
      result$matrix <- result$time_matrix
      result$edge_list <- result$time_edge_list
    }
  }, error = function(e) {
    warning("Error setting default view: ", e$message, ". Using time view.")
    result$matrix <- result$time_matrix
    result$edge_list <- result$time_edge_list
  })
  
  return(result)
}
