#' Build Network from Distance Matrix
#'
#' This is a core function to construct various types of networks from a given
#' distance matrix. It supports different methods for defining connections
#' between nodes based on their distances.
#'
#' @param distance_matrix A square numeric matrix representing the distances
#'   between time series segments or windows. Typically, this would be the output
#'   of `analyze_ts_distances()` or `ts_distance()`.
#' @param method Character string. The network construction method to use.
#'   Options include:
#'   \itemize{
#'     \item \code{"full"} (default): Returns a raw similarity matrix without
#'       applying any specific network processing. This is useful for direct
#'       visualization with tools like `qgraph` where edge filtering is handled
#'       separately.
#'     \item \code{"knn"}: K-nearest neighbors graph. Each node is connected
#'       to its `k` closest neighbors.
#'     \item \code{"threshold"}: Connects nodes whose distance is below a
#'       specified `threshold`.
#'     \item \code{"percentile"}: Connects nodes whose distance falls within
#'       a specified `percentile` range (e.g., the smallest 10% of distances).
#'     \item \code{"gaussian"}: Applies a Gaussian kernel to transform distances
#'       into similarity weights, creating a fully connected weighted graph.
#'     \item \code{"epsilon"}: Epsilon neighborhood graph. Connects nodes whose
#'       distance is less than or equal to `epsilon`.
#'   }
#' @param k Integer. The number of nearest neighbors to consider for the
#'   \code{"knn"} method. Defaults to 5.
#' @param epsilon Numeric. The maximum distance for connections in the
#'   \code{"epsilon"} method. Defaults to 0.1.
#' @param threshold Numeric. The maximum distance for connections in the
#'   \code{"threshold"} method. If `NULL` for this method, it defaults to the
#'   25th percentile of finite, positive distances.
#' @param sigma Numeric. The bandwidth parameter for the Gaussian kernel in the
#'   \code{"gaussian"} method. If `NULL`, it defaults to the median of finite,
#'   positive distances.
#' @param percentile Numeric. The percentile (0 to 1) of smallest distances to
#'   consider for connections in the \code{"percentile"} method. Defaults to 0.1.
#' @param directed Logical. If `TRUE`, the resulting graph will be directed.
#'   Defaults to `FALSE`.
#' @param self_loops Logical. If `TRUE`, self-loops (connections from a node
#'   to itself) are allowed. Defaults to `FALSE`.
#' @param normalize Logical. If `TRUE`, distances are normalized to a [0, 1]
#'   range before applying network construction methods (except for "gaussian"
#'   which handles its own scaling). Defaults to `FALSE`.
#'
#' @return A list containing the constructed network information:
#'   \item{graph}{An `igraph` object representing the network (for non-"full" methods).}
#'   \item{similarity_matrix}{The similarity matrix (for "full" method).}
#'   \item{adjacency_matrix}{The adjacency matrix of the constructed network.}
#'   \item{original_distance_matrix}{The input distance matrix.}
#'   \item{method_used}{The network construction method that was applied.}
#'   \item{parameters}{A list of parameters used for the chosen method.}
#'   \item{network_stats}{A list of basic network statistics (e.g., number of
#'     nodes, edges, density, connectivity).}
#'   \item{nodes}{A data frame with node information, including `window_id` and `label`.
#'     This is useful for compatibility with visualization packages like `qgraph`.}
#'   \item{edges}{A data frame representing the edge list of the network, including
#'     `from`, `to`, and `weight` columns.}
#' @importFrom igraph graph_from_adjacency_matrix V set_graph_attr vcount ecount edge_density is_connected count_components as_data_frame
#' @importFrom stats median quantile

#' @export
build_network <- function(distance_matrix, 
                          method = "full",  # Changed default to "full"
                          k = 5, 
                          epsilon = 0.1, 
                          threshold = 0.5,
                          sigma = NULL,
                          percentile = 0.1,
                          directed = FALSE,
                          self_loops = FALSE,
                          normalize = FALSE) {
  
  # Input validation
  if (!is.matrix(distance_matrix)) {
    stop("Input must be a distance matrix")
  }
  
  # Check if matrix is square
  if (nrow(distance_matrix) != ncol(distance_matrix)) {
    stop("Distance matrix must be square")
  }
  
  n <- nrow(distance_matrix)
  
  # NEW: Handle "full" method - return raw similarities without network processing
  if (method == "full" || method == "raw" || method == "similarities") {
    cat("Method 'full': Returning raw similarity matrix without network processing\n")
    
    # Convert distances to similarities using prepare_for_qgraph
    if (!requireNamespace("tsn", quietly = TRUE)) {
      # Fallback similarity conversion if prepare_for_qgraph not available
      similarity_matrix <- 1 - (distance_matrix / max(distance_matrix[is.finite(distance_matrix)]))
      similarity_matrix[!is.finite(similarity_matrix)] <- 0
      diag(similarity_matrix) <- 0
    } else {
      # Use prepare_for_qgraph for proper similarity conversion
      similarity_matrix <- prepare_for_qgraph(distance_matrix)
    }
    
    # Create minimal result for "full" method
    result <- list(
      method_used = "full",
      similarity_matrix = similarity_matrix,
      adjacency_matrix = similarity_matrix,  # Same as similarity for "full"
      original_distance_matrix = distance_matrix,
      parameters = list(method = "full"),
      network_stats = list(
        nodes = n,
        edges = sum(similarity_matrix > 0),
        density = sum(similarity_matrix > 0) / (n * (n - 1)),
        is_full_similarity = TRUE
      ),
      # For qgraph compatibility
      nodes = data.frame(
        window_id = 1:n,
        label = if (!is.null(rownames(distance_matrix))) rownames(distance_matrix) else paste0("W", 1:n),
        stringsAsFactors = FALSE
      )
    )
    
    cat("[+] Raw similarity matrix ready for qgraph visualization
")
    cat("Use: qgraph(result$similarity_matrix, layout='spring')
")
    
    return(result)
  }
  
  # Continue with existing network methods for non-"full" methods
  
  # Load required libraries for network methods
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for network methods (not needed for 'full' method)")
  }
  
  # Handle infinity values - replace with large finite value for network construction
  if (any(is.infinite(distance_matrix))) {
    finite_dists <- distance_matrix[is.finite(distance_matrix) & distance_matrix > 0]
    if (length(finite_dists) > 0) {
      max_finite <- max(finite_dists)
      distance_matrix[is.infinite(distance_matrix)] <- max_finite * 10
    }
  }
  
  # Normalize distances if requested
  if (normalize) {
    max_dist <- max(distance_matrix[distance_matrix != Inf])
    distance_matrix <- distance_matrix / max_dist
  }
  
  # Initialize adjacency matrix
  adj_matrix <- matrix(0, nrow = n, ncol = n)
  
  # Build network based on specified method
  if (method == "epsilon") {
    # Epsilon neighborhood graph
    adj_matrix[distance_matrix <= epsilon & distance_matrix > 0] <- 1
    
  } else if (method == "knn") {
    # K-nearest neighbors - works well with infinity values
    for (i in 1:n) {
      finite_mask <- is.finite(distance_matrix[i, ]) & distance_matrix[i, ] > 0
      if (sum(finite_mask) > 0) {
        finite_indices <- which(finite_mask)
        finite_distances <- distance_matrix[i, finite_indices]
        k_actual <- min(k, length(finite_distances))
        if (k_actual > 0) {
          nearest_indices <- finite_indices[order(finite_distances)[1:k_actual]]
          adj_matrix[i, nearest_indices] <- 1
        }
      }
    }
    
  } else if (method == "threshold") {
    # Threshold-based connections
    if (is.null(threshold)) {
      finite_dists <- distance_matrix[is.finite(distance_matrix) & distance_matrix > 0]
      threshold <- quantile(finite_dists, 0.25)
    }
    adj_matrix[distance_matrix <= threshold & distance_matrix > 0] <- 1
    
  } else if (method == "percentile") {
    # Connect smallest percentile of distances
    finite_dists <- distance_matrix[is.finite(distance_matrix) & distance_matrix > 0]
    if (length(finite_dists) > 0) {
      threshold_val <- quantile(finite_dists, percentile)
      adj_matrix[distance_matrix <= threshold_val & distance_matrix > 0] <- 1
    }
    
  } else if (method == "gaussian") {
    # Gaussian kernel transformation
    if (is.null(sigma)) {
      finite_dists <- distance_matrix[is.finite(distance_matrix) & distance_matrix > 0]
      sigma <- median(finite_dists)
    }
    adj_matrix <- exp(-distance_matrix^2 / (2 * sigma^2))
    diag(adj_matrix) <- 0
    
  } else {
    stop("Unknown method. Available methods: full (default), knn, threshold, percentile, gaussian, epsilon")
  }
  
  # Remove self-loops if specified
  if (!self_loops) {
    diag(adj_matrix) <- 0
  }
  
  # Create igraph object for network methods
  if (method == "gaussian") {
    # Weighted graph for gaussian (weights are already similarities)
    graph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = ifelse(directed, "directed", "undirected"), weighted = TRUE)
  } else {
    # For other methods, convert binary adjacency to weighted using original distances
    # Convert distances to similarities for edge weights
    weighted_adj <- adj_matrix
    connected_mask <- adj_matrix > 0
    
    if (any(connected_mask)) {
      # Use inverse distance as weight (higher = more similar)
      distance_vals <- distance_matrix[connected_mask]
      # Avoid division by zero
      distance_vals[distance_vals == 0] <- min(distance_vals[distance_vals > 0]) / 1000
      similarity_vals <- 1 / (1 + distance_vals)  # Convert to similarity
      weighted_adj[connected_mask] <- similarity_vals
    }
    
    graph <- igraph::graph_from_adjacency_matrix(weighted_adj, mode = ifelse(directed, "directed", "undirected"), weighted = TRUE)
  }
  
  # Add node names if available
  if (!is.null(rownames(distance_matrix))) {
    igraph::V(graph)$name <- rownames(distance_matrix)
    igraph::V(graph)$label <- rownames(distance_matrix)
  } else {
    igraph::V(graph)$name <- paste0("TS_", 1:n)
    igraph::V(graph)$label <- paste0("TS_", 1:n)
  }
  
  # Add metadata
  graph <- igraph::set_graph_attr(graph, "method", method)
  graph <- igraph::set_graph_attr(graph, "k", k)
  graph <- igraph::set_graph_attr(graph, "epsilon", epsilon)
  graph <- igraph::set_graph_attr(graph, "threshold", threshold)
  graph <- igraph::set_graph_attr(graph, "percentile", percentile)
  
  return(list(
    graph = graph,
    adjacency_matrix = adj_matrix,
    original_distance_matrix = distance_matrix,
    method_used = method,
    parameters = list(k = k, epsilon = epsilon, threshold = threshold, 
                      sigma = sigma, percentile = percentile),
    network_stats = list(
      nodes = igraph::vcount(graph),
      edges = igraph::ecount(graph),
      density = igraph::edge_density(graph),
      is_connected = igraph::is_connected(graph),
      components = igraph::count_components(graph)
    ),
    # Add node and edge information in ts_network format
    nodes = data.frame(
      window_id = 1:n,
      label = if (!is.null(rownames(distance_matrix))) rownames(distance_matrix) else paste0("W", 1:n),
      stringsAsFactors = FALSE
    ),
    edges = if (igraph::ecount(graph) > 0) {
      edges_df <- igraph::as_data_frame(graph, what = "edges")
      edges_df$weight <- if ("weight" %in% names(edges_df)) edges_df$weight else 1
      edges_df
    } else {
      data.frame(from = character(0), to = character(0), weight = numeric(0), stringsAsFactors = FALSE)
    }
  ))
}