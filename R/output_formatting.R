# Output Formatting Functions
# Clean, network-ready output formatting for time series distance analysis

#' Create network-ready distance matrix from time series analysis
#' @param analysis_type Type of analysis: "windowed" or "id_partition"
#' @param segments List of segments with metadata
#' @param distances Vector of calculated distances
#' @param method_name Method used for calculation
#' @param output_format Format of the output, either "matrix" or "edges".
#' @return Clean distance matrix or dataframe

create_network_output <- function(analysis_type, segments, distances, method_name, output_format = "matrix") {
  
  if (analysis_type == "windowed") {
    return(create_windowed_network_output(segments, distances, method_name, output_format))
  } else if (analysis_type == "id_partition") {
    return(create_id_network_output(segments, distances, method_name, output_format))
  } else {
    stop("Unknown analysis type")
  }
}

#' Create windowed analysis network output
#' @param segments List of segments with metadata
#' @param distances Vector of calculated distances
#' @param method_name Method used for calculation
#' @param output_format Format of the output, either "matrix" or "edges".
create_windowed_network_output <- function(segments, distances, method_name, output_format) {
  n_windows <- length(segments$windows)
  
  if (output_format == "matrix") {
    # Create adjacency matrix for network
    adj_matrix <- matrix(0, nrow = n_windows, ncol = n_windows)
    
    # Fill consecutive window distances
    for (i in 1:(n_windows - 1)) {
      adj_matrix[i, i + 1] <- distances[i]
      adj_matrix[i + 1, i] <- distances[i]  # Symmetric
    }
    
    # Add row/column names
    window_names <- paste0("W", 1:n_windows, "_", segments$starts, ":", segments$ends)
    rownames(adj_matrix) <- window_names
    colnames(adj_matrix) <- window_names
    
    return(list(
      adjacency_matrix = adj_matrix,
      window_info = data.frame(
        window_id = 1:n_windows,
        start_index = segments$starts,
        end_index = segments$ends,
        window_size = segments$lengths,
        stringsAsFactors = FALSE
      ),
      method = method_name,
      type = "windowed"
    ))
    
  } else if (output_format == "edges") {
    # Create edge list for network
    edges_df <- data.frame(
      from_window = 1:(n_windows - 1),
      to_window = 2:n_windows,
      distance = distances,
      from_start = segments$starts[1:(n_windows - 1)],
      from_end = segments$ends[1:(n_windows - 1)],
      to_start = segments$starts[2:n_windows],
      to_end = segments$ends[2:n_windows],
      method = method_name,
      stringsAsFactors = FALSE
    )
    
    return(list(
      edges = edges_df,
      nodes = data.frame(
        window_id = 1:n_windows,
        start_index = segments$starts,
        end_index = segments$ends,
        window_size = segments$lengths,
        stringsAsFactors = FALSE
      ),
      method = method_name,
      type = "windowed"
    ))
  }
}

#' Create ID-based analysis network output
#' @param segments List of segments with metadata
#' @param distances Vector of calculated distances
#' @param method_name Method used for calculation
#' @param output_format Format of the output, either "matrix" or "edges".
create_id_network_output <- function(segments, distances, method_name, output_format) {
  n_segments <- length(segments$ids)
  
  if (output_format == "matrix") {
    # Create adjacency matrix
    adj_matrix <- matrix(0, nrow = n_segments, ncol = n_segments)
    
    # Fill consecutive segment distances
    for (i in 1:(n_segments - 1)) {
      adj_matrix[i, i + 1] <- distances[i]
      adj_matrix[i + 1, i] <- distances[i]  # Symmetric
    }
    
    # Add names
    rownames(adj_matrix) <- segments$ids
    colnames(adj_matrix) <- segments$ids
    
    return(list(
      adjacency_matrix = adj_matrix,
      segment_info = data.frame(
        segment_id = segments$ids,
        start_index = segments$starts,
        end_index = segments$ends,
        segment_size = segments$ends - segments$starts + 1,
        stringsAsFactors = FALSE
      ),
      method = method_name,
      type = "id_partition"
    ))
    
  } else if (output_format == "edges") {
    # Create edge list
    edges_df <- data.frame(
      from_segment = segments$ids[1:(n_segments - 1)],
      to_segment = segments$ids[2:n_segments],
      distance = distances,
      from_start = segments$starts[1:(n_segments - 1)],
      from_end = segments$ends[1:(n_segments - 1)],
      to_start = segments$starts[2:n_segments],
      to_end = segments$ends[2:n_segments],
      method = method_name,
      stringsAsFactors = FALSE
    )
    
    return(list(
      edges = edges_df,
      nodes = data.frame(
        segment_id = segments$ids,
        start_index = segments$starts,
        end_index = segments$ends,
        segment_size = segments$ends - segments$starts + 1,
        stringsAsFactors = FALSE
      ),
      method = method_name,
      type = "id_partition"
    ))
  }
}

#' Create point-wise mapping for expanded format
#' @param original_series The original time series
#' @param network_result Network result from analysis
#' @return Data frame with one row per point
create_point_mapping <- function(original_series, network_result) {
  n_points <- length(original_series)
  
  if (!is.null(network_result$type) && network_result$type == "windowed") {
    # Handle windowed analysis
    if (!is.null(network_result$nodes)) {
      nodes <- network_result$nodes
    } else {
      nodes <- network_result$window_info
    }
    
    point_mapping <- data.frame(
      point_index = 1:n_points,
      ts_value = original_series,
      window_id = NA,
      position_in_window = NA,
      stringsAsFactors = FALSE
    )
    
    # Add window membership
    for (i in 1:nrow(nodes)) {
      start_idx <- nodes$start_index[i]
      end_idx <- nodes$end_index[i]
      window_id <- if("window_id" %in% names(nodes)) nodes$window_id[i] else i
      
      # Find points in this window
      in_window <- point_mapping$point_index >= start_idx & point_mapping$point_index <= end_idx
      
      if (any(in_window)) {
        existing_windows <- point_mapping$window_id[in_window]
        point_mapping$window_id[in_window] <- ifelse(is.na(existing_windows), 
                                                   window_id, 
                                                   paste(existing_windows, window_id, sep = ","))
        positions <- point_mapping$point_index[in_window] - start_idx + 1
        existing_positions <- point_mapping$position_in_window[in_window]
        point_mapping$position_in_window[in_window] <- ifelse(is.na(existing_positions),
                                                            positions,
                                                            paste(existing_positions, positions, sep = ","))
      }
    }
    
  } else {
    # Handle ID-based analysis
    if (!is.null(network_result$nodes)) {
      nodes <- network_result$nodes
    } else {
      nodes <- network_result$segment_info
    }
    
    point_mapping <- data.frame(
      point_index = 1:n_points,
      ts_value = original_series,
      segment_id = NA,
      position_in_segment = NA,
      stringsAsFactors = FALSE
    )
    
    # Add segment membership
    for (i in 1:nrow(nodes)) {
      start_idx <- nodes$start_index[i]
      end_idx <- nodes$end_index[i]
      segment_id <- if("segment_id" %in% names(nodes)) nodes$segment_id[i] else paste0("seg_", i)
      
      in_segment <- point_mapping$point_index >= start_idx & point_mapping$point_index <= end_idx
      if (any(in_segment)) {
        point_mapping$segment_id[in_segment] <- segment_id
        point_mapping$position_in_segment[in_segment] <- (point_mapping$point_index[in_segment] - start_idx + 1)
      }
    }
  }
  
  return(point_mapping)
}

#' Create Summary Statistics for Network Analysis
#'
#' This internal helper function generates summary statistics for the results
#' of a time series network analysis, providing an overview of the calculated
#' distances and network structure.
#'
#' @param network_result List. The result object from a time series network analysis
#'   (e.g., from `analyze_ts_distances`), containing `type`, `method`, `edges`, and `nodes`.
#' @return A list containing various summary statistics, such as number of windows/segments,
#'   number of edges, and statistics of the distances.
#' @keywords internal
create_network_summary <- function(network_result) {
  if (network_result$type == "windowed") {
    edges <- network_result$edges
    nodes <- network_result$nodes
    
    summary_stats <- list(
      type = "windowed",
      method = network_result$method,
      n_windows = nrow(nodes),
      n_edges = nrow(edges),
      mean_distance = mean(edges$distance, na.rm = TRUE),
      median_distance = median(edges$distance, na.rm = TRUE),
      min_distance = min(edges$distance, na.rm = TRUE),
      max_distance = max(edges$distance, na.rm = TRUE),
      window_size_range = paste(min(nodes$window_size), "-", max(nodes$window_size))
    )
    
  } else if (network_result$type == "id_partition") {
    edges <- network_result$edges
    nodes <- network_result$nodes
    
    summary_stats <- list(
      type = "id_partition",
      method = network_result$method,
      n_segments = nrow(nodes),
      n_edges = nrow(edges),
      mean_distance = mean(edges$distance, na.rm = TRUE),
      median_distance = median(edges$distance, na.rm = TRUE),
      min_distance = min(edges$distance, na.rm = TRUE),
      max_distance = max(edges$distance, na.rm = TRUE),
      segment_size_range = paste(min(nodes$segment_size), "-", max(nodes$segment_size))
    )
  }
  
  return(summary_stats)
}
