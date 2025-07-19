#' Time Series Network Construction
#' 
#' Build networks from distance matrices using various construction methods
#' 
#' @param distance_matrix Symmetric distance matrix
#' @param method Network construction method ("full", "knn", "percentile", "threshold", "gaussian", "epsilon")
#'   - "full": Returns raw similarity matrix (DEFAULT) - no network processing, ready for qgraph
#'   - "knn": K-nearest neighbors network
#'   - "percentile": Percentile threshold network
#'   - "threshold": Distance threshold network
#'   - "gaussian": Gaussian kernel network
#'   - "epsilon": Epsilon neighborhood network
#' @param k Number of nearest neighbors (for knn method)
#' @param percentile Percentile threshold (for percentile method) 
#' @param threshold Distance threshold (for threshold method)
#' @param sigma Gaussian kernel parameter (for gaussian method)
#' @param epsilon Epsilon neighborhood size (for epsilon method)
#' @param normalize Whether to normalize distances first
#' 
#' @return Network object with similarity matrix and optional graph for network methods
#' @export
ts_network <- function(distance_matrix, 
                      method = "full",  # Changed default to "full"
                      k = 5,
                      percentile = 0.2,
                      threshold = NULL,
                      sigma = NULL,
                      epsilon = 0.1,
                      normalize = FALSE) {
  
  # Validate distance matrix
  if (!is.matrix(distance_matrix) || nrow(distance_matrix) != ncol(distance_matrix)) {
    stop("distance_matrix must be a square matrix")
  }
  
  # Use the existing build_network function
  result <- build_network(
    distance_matrix = distance_matrix,
    method = method,
    k = k,
    percentile = percentile,
    threshold = threshold,
    sigma = sigma,
    epsilon = epsilon,
    normalize = normalize
  )
  
  return(result)
}

#' Plot Time Series Network using qgraph
#' 
#' Visualize time series networks using qgraph (publication quality)
#' 
#' @param network_result Network result from ts_network() or build_network()
#' @param layout Layout method ("spring", "circle", "random")
#' @param threshold Threshold for edge display (for "full" method)
#' @param node.color Node color
#' @param node.size Node size
#' @param edge.color Edge color  
#' @param title Plot title
#' @param ... Additional arguments passed to qgraph
#' 
#' @return qgraph plot object (invisible)
#' @export
plot_ts_network <- function(network_result, 
                           layout = "spring", 
                           threshold = 0.1,
                           node.color = "lightblue",
                           node.size = 8,
                           edge.color = "gray50",
                           title = NULL,
                           ...) {
  
  # Load qgraph
  if (!requireNamespace("qgraph", quietly = TRUE)) {
    stop("Package 'qgraph' is required. Install with: install.packages('qgraph')")
  }
  
  # Set default title
  if (is.null(title)) {
    if (network_result$method_used == "full") {
      title <- paste("Time Series Network (Raw Similarities)\nNodes:", 
                     network_result$network_stats$nodes)
    } else {
      title <- paste("Time Series Network -", toupper(network_result$method_used),
                     "\nNodes:", network_result$network_stats$nodes, 
                     "| Edges:", network_result$network_stats$edges)
    }
  }
  
  # Choose matrix to plot
  if (network_result$method_used == "full") {
    # For "full" method, use similarity matrix directly
    plot_matrix <- network_result$similarity_matrix
    cat("Plotting raw similarity matrix with qgraph\n")
  } else {
    # For network methods, use adjacency matrix
    plot_matrix <- network_result$adjacency_matrix
    cat("Plotting network adjacency matrix with qgraph\n")
  }
  
  # Set node labels
  if (!is.null(network_result$nodes)) {
    node_labels <- network_result$nodes$label
  } else if (!is.null(rownames(plot_matrix))) {
    node_labels <- rownames(plot_matrix)
  } else {
    node_labels <- paste0("W", 1:nrow(plot_matrix))
  }
  
  # Create qgraph plot
  plot_result <- qgraph::qgraph(
    plot_matrix,
    layout = layout,
    threshold = threshold,
    color = node.color,
    vsize = node.size,
    edge.color = edge.color,
    title = title,
    labels = node_labels,
    label.cex = 0.8,
    ...
  )
  
  # Print info
  cat("[+] qgraph plot created successfully\n")
  if (network_result$method_used == "full") {
    cat("Showing raw similarities - adjust 'threshold' parameter to control edge display\n")
  }
  
  return(invisible(plot_result))
}

#' Get Available Network Methods
#' 
#' Returns list of all available network construction methods
#' 
#' @return Character vector of method names
#' @export
ts_network_methods <- function() {
  return(c("full", "knn", "percentile", "threshold", "gaussian", "epsilon"))
}

#' Recommend Network Method for Distance Scale
#' 
#' Suggests appropriate network method based on distance distribution
#' 
#' @param distance_matrix Distance matrix
#' @return List with recommended method and parameters
#' @export
ts_network_recommend <- function(distance_matrix) {
  
  # Always recommend "full" as the most flexible default
  recommendation <- list(
    method = "full",
    parameters = list(),
    reason = "Raw similarities provide maximum flexibility for visualization - RECOMMENDED DEFAULT"
  )
  
  # Provide alternative recommendations based on distance scale
  distances <- distance_matrix[upper.tri(distance_matrix) & distance_matrix > 0 & is.finite(distance_matrix)]
  
  if (length(distances) > 0) {
    max_dist <- max(distances)
    median_dist <- median(distances)
    
    alternative <- if (max_dist > 1000) {
      list(
        method = "knn",
        parameters = list(k = 5),
        reason = "Large distances detected - KNN is scale-invariant"
      )
    } else if (max_dist > 100) {
      list(
        method = "percentile", 
        parameters = list(percentile = 0.2),
        reason = "Medium-large distances - percentile method recommended"
      )
    } else {
      list(
        method = "threshold",
        parameters = list(threshold = quantile(distances, 0.25)),
        reason = "Moderate distances - threshold method suitable"
      )
    }
    
    recommendation$alternative <- alternative
  }
  
  return(recommendation)
}

#' Quick Network from Distance Matrix
#'
#' A convenience function for one-step network creation from a distance matrix.
#' It can automatically select a suitable network construction method based on
#' data characteristics or use a default method.
#'
#' @param distance_matrix A numeric distance matrix.
#' @param auto_select Logical. If `TRUE`, the function will automatically select
#'   a network method using `ts_network_recommend()`. If `FALSE` (default),
#'   it will use the "full" method by default.
#'
#' @return A network object, which is a list containing the constructed network
#'   information and the recommendation details.
#' @seealso \code{\link{ts_network_recommend}}, \code{\link{ts_network}}
#' @export
ts_network_auto <- function(distance_matrix, auto_select = FALSE) {
  
  if (auto_select) {
    recommendation <- ts_network_recommend(distance_matrix)
    method <- recommendation$alternative$method  # Use alternative, not the "full" default
    params <- recommendation$alternative$parameters
    
    result <- do.call(ts_network, c(list(distance_matrix = distance_matrix, method = method), params))
    result$recommendation <- recommendation
  } else {
    # Use "full" method by default
    result <- ts_network(distance_matrix, method = "full")
    result$recommendation <- list(
      method = "full", 
      reason = "Default raw similarities - most flexible for qgraph visualization"
    )
  }
  
  return(result)
}
