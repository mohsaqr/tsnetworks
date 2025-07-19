#' Distance to Similarity Conversion and Normalization Utilities
#' 
#' Functions for converting distances to similarities and normalizing for network visualization
#' 

#' Convert Distance Matrix to Similarity Matrix
#' 
#' Transform distance values to similarity values using various methods.
#' Higher similarity values indicate more similar (closer) time series.
#' 
#' @param distance_matrix Square distance matrix
#' @param method Conversion method:
#'   - "normalized_inverse": similarity = 1 - (normalized_distance) [DEFAULT]
#'   - "inverse": similarity = 1 / (1 + distance) 
#'   - "negative_exp": similarity = exp(-distance/sigma)
#'   - "max_minus": similarity = (max_distance - distance) / max_distance
#'   - "gaussian": similarity = exp(-distance^2 / (2*sigma^2))
#'   - "reciprocal": similarity = 1 / distance (for distance > 0)
#' @param sigma Scale parameter for exponential/gaussian methods (auto-detected if NULL)
#' @param remove_self_sim Remove self-similarities (set diagonal to 0)
#' 
#' @return Similarity matrix where higher values = more similar
#' @export
#' @examples
#' # Create example distance matrix
#' distances <- matrix(c(0, 10, 20, 10, 0, 15, 20, 15, 0), nrow = 3)
#' 
#' # Convert to similarities
#' sim_norm_inv <- distance_to_similarity(distances, "normalized_inverse")
#' sim_inverse <- distance_to_similarity(distances, "inverse")
#' sim_exp <- distance_to_similarity(distances, "negative_exp")
#' sim_gaussian <- distance_to_similarity(distances, "gaussian")
#' @importFrom stats median mad quantile sd
#' @importFrom qgraph qgraph
distance_to_similarity <- function(distance_matrix, 
                                  method = "normalized_inverse",
                                  sigma = NULL,
                                  remove_self_sim = TRUE) {
  
  # Validate input
  if (!is.matrix(distance_matrix) || nrow(distance_matrix) != ncol(distance_matrix)) {
    stop("distance_matrix must be a square matrix")
  }
  
  # Get non-zero distances for parameter estimation
  non_zero_distances <- distance_matrix[distance_matrix > 0 & is.finite(distance_matrix)]
  
  if (length(non_zero_distances) == 0) {
    warning("No positive finite distances found")
    return(distance_matrix)
  }
  
  # Auto-detect sigma if needed
  if (is.null(sigma) && method %in% c("negative_exp", "gaussian")) {
    sigma <- median(non_zero_distances)
    if (sigma == 0) sigma <- mean(non_zero_distances)
    if (sigma == 0) sigma <- 1
  }
  
  # Convert based on method
  similarity_matrix <- switch(method,
    "normalized_inverse" = {
      # NEW DEFAULT: 1 - normalized distance
      # First normalize distances to [0,1] range
      max_dist <- max(distance_matrix[is.finite(distance_matrix)])
      if (max_dist == 0) {
        matrix(1, nrow = nrow(distance_matrix), ncol = ncol(distance_matrix))
      } else {
        normalized_dist <- distance_matrix / max_dist
        1 - normalized_dist
      }
    },
    
    "inverse" = {
      # Simple inverse: similarity = 1 / (1 + distance)
      # Works well for any distance scale, always gives values in (0,1]
      1 / (1 + distance_matrix)
    },
    
    "negative_exp" = {
      # Negative exponential: similarity = exp(-distance/sigma)
      # Good for moderate distance ranges
      exp(-distance_matrix / sigma)
    },
    
    "max_minus" = {
      # Maximum minus distance, normalized
      # similarity = (max_distance - distance) / max_distance
      max_dist <- max(distance_matrix[is.finite(distance_matrix)])
      if (max_dist == 0) {
        matrix(1, nrow = nrow(distance_matrix), ncol = ncol(distance_matrix))
      } else {
        (max_dist - distance_matrix) / max_dist
      }
    },
    
    "gaussian" = {
      # Gaussian kernel: similarity = exp(-distance^2 / (2*sigma^2))
      # Good for creating smooth similarity transitions
      exp(-distance_matrix^2 / (2 * sigma^2))
    },
    
    "reciprocal" = {
      # Simple reciprocal: similarity = 1 / distance (only for distance > 0)
      # Can produce very large values for small distances
      result <- matrix(0, nrow = nrow(distance_matrix), ncol = ncol(distance_matrix))
      positive_mask <- distance_matrix > 0
      result[positive_mask] <- 1 / distance_matrix[positive_mask]
      result
    },
    
    stop("Unknown method. Choose from: normalized_inverse, inverse, negative_exp, max_minus, gaussian, reciprocal")
  )
  
  # Handle infinities and NaNs
  similarity_matrix[!is.finite(similarity_matrix)] <- 0
  
  # Remove self-similarities if requested
  if (remove_self_sim) {
    diag(similarity_matrix) <- 0
  }
  
  # Add attributes for reference
  attr(similarity_matrix, "method") <- method
  attr(similarity_matrix, "sigma") <- sigma
  attr(similarity_matrix, "original_distance_range") <- range(non_zero_distances)
  
  return(similarity_matrix)
}

#' Normalize Matrix for Visualization
#' 
#' Normalize matrix values to appropriate range for network visualization
#' 
#' @param matrix Input matrix (distance or similarity)
#' @param method Normalization method:
#'   - "minmax": (x - min) / (max - min) → [0,1]
#'   - "max": x / max(x) → [0,1] 
#'   - "zscore": (x - mean) / sd → standardized
#'   - "robust": (x - median) / mad → robust standardization
#'   - "quantile": scale to quantile range
#' @param target_range Target range for output (e.g., c(0,1))
#' @param exclude_diagonal Exclude diagonal from normalization calculations
#' @param preserve_zeros Whether to keep zero values as zero
#' 
#' @return Normalized matrix
#' @export
#' @examples
#' # Large distance matrix
#' distances <- matrix(runif(9, 1000, 50000), nrow = 3)
#' 
#' # Normalize to 0-1 range
#' norm_minmax <- normalize_matrix(distances, "minmax")
#' norm_max <- normalize_matrix(distances, "max")
normalize_matrix <- function(matrix, 
                            method = "minmax",
                            target_range = c(0, 1),
                            exclude_diagonal = TRUE,
                            preserve_zeros = TRUE) {
  
  if (!is.matrix(matrix)) {
    stop("Input must be a matrix")
  }
  
  # Get values for normalization (excluding diagonal if requested)
  if (exclude_diagonal && nrow(matrix) == ncol(matrix)) {
    # Square matrix - exclude diagonal
    mask <- !diag(nrow(matrix))
    values_for_calc <- matrix[mask]
  } else {
    values_for_calc <- as.vector(matrix)
  }
  
  # Exclude zeros if preserving them
  if (preserve_zeros) {
    non_zero_mask <- values_for_calc != 0
    calc_values <- values_for_calc[non_zero_mask]
  } else {
    calc_values <- values_for_calc
  }
  
  if (length(calc_values) == 0) {
    warning("No values available for normalization")
    return(matrix)
  }
  
  # Calculate normalization parameters
  if (method == "minmax") {
    min_val <- min(calc_values, na.rm = TRUE)
    max_val <- max(calc_values, na.rm = TRUE)
    if (max_val == min_val) {
      normalized <- matrix
    } else {
      normalized <- (matrix - min_val) / (max_val - min_val)
    }
    
  } else if (method == "max") {
    max_val <- max(calc_values, na.rm = TRUE)
    if (max_val == 0) {
      normalized <- matrix
    } else {
      normalized <- matrix / max_val
    }
    
  } else if (method == "zscore") {
    mean_val <- mean(calc_values, na.rm = TRUE)
    sd_val <- sd(calc_values, na.rm = TRUE)
    if (sd_val == 0) {
      normalized <- matrix - mean_val
    } else {
      normalized <- (matrix - mean_val) / sd_val
    }
    
  } else if (method == "robust") {
    median_val <- median(calc_values, na.rm = TRUE)
    mad_val <- mad(calc_values, na.rm = TRUE)
    if (mad_val == 0) {
      normalized <- matrix - median_val
    } else {
      normalized <- (matrix - median_val) / mad_val
    }
    
  } else if (method == "quantile") {
    q05 <- quantile(calc_values, 0.05, na.rm = TRUE)
    q95 <- quantile(calc_values, 0.95, na.rm = TRUE)
    if (q95 == q05) {
      normalized <- matrix
    } else {
      normalized <- (matrix - q05) / (q95 - q05)
    }
    
  } else {
    stop("Unknown method. Choose from: minmax, max, zscore, robust, quantile")
  }
  
  # Preserve zeros if requested
  if (preserve_zeros) {
    zero_mask <- matrix == 0
    normalized[zero_mask] <- 0
  }
  
  # Scale to target range if not zscore/robust
  if (!method %in% c("zscore", "robust") && !is.null(target_range)) {
    current_min <- min(normalized, na.rm = TRUE)
    current_max <- max(normalized, na.rm = TRUE)
    
    if (current_max > current_min) {
      normalized <- (normalized - current_min) / (current_max - current_min)
      normalized <- normalized * (target_range[2] - target_range[1]) + target_range[1]
    }
  }
  
  return(normalized)
}

#' Prepare Matrix for qgraph Visualization
#' 
#' One-stop function to prepare distance matrices for qgraph visualization.
#' Combines distance-to-similarity conversion and normalization.
#' 
#' @param distance_matrix Distance matrix
#' @param conversion_method Method for distance-to-similarity conversion
#' @param normalization Whether to normalize ("none", "minmax", "max")  
#' @param target_range Target range for final values
#' @param threshold_method How to handle thresholding ("none", "percentile", "absolute")
#' @param threshold_value Threshold value (percentile 0-1 or absolute value)
#' 
#' @return Matrix ready for qgraph visualization with metadata
#' @export
#' @examples
#' # Large distance matrix (17,000 to 50,000 range)
#' set.seed(123)
#' large_distances <- matrix(runif(25, 17000, 50000), nrow = 5)
#' large_distances <- (large_distances + t(large_distances)) / 2  # Make symmetric
#' diag(large_distances) <- 0
#' 
#' # Prepare for qgraph
#' viz_matrix <- prepare_for_qgraph(large_distances)
#' 
#' # Use with qgraph
#' if (requireNamespace("qgraph", quietly = TRUE)) {
#'   qgraph::qgraph(viz_matrix, layout = "spring", threshold = 0.7)
#' }
prepare_for_qgraph <- function(distance_matrix,
                              conversion_method = "normalized_inverse",
                              normalization = "minmax",
                              target_range = c(0, 1),
                              threshold_method = "none",
                              threshold_value = 0.7) {
  
  # Step 1: Convert distances to similarities
  if (conversion_method != "none") {
    cat("Converting distances to similarities using method:", conversion_method, "\n")
    similarity_matrix <- distance_to_similarity(distance_matrix, method = conversion_method)
  } else {
    similarity_matrix <- distance_matrix
  }
  
  # Step 2: Normalize if requested
  if (normalization != "none") {
    cat("Normalizing using method:", normalization, "\n")
    similarity_matrix <- normalize_matrix(similarity_matrix, method = normalization, target_range = target_range)
  }
  
  # Step 3: Apply threshold if requested
  if (threshold_method == "percentile") {
    threshold_val <- quantile(similarity_matrix[similarity_matrix > 0], threshold_value, na.rm = TRUE)
    similarity_matrix[similarity_matrix < threshold_val] <- 0
    cat("Applied percentile threshold:", threshold_value, "-> value:", round(threshold_val, 4), "\n")
    
  } else if (threshold_method == "absolute") {
    similarity_matrix[similarity_matrix < threshold_value] <- 0
    cat("Applied absolute threshold:", threshold_value, "\n")
  }
  
  # Add diagnostic information
  non_zero <- similarity_matrix[similarity_matrix > 0]
  cat("Final matrix range:", round(min(non_zero), 4), "to", round(max(non_zero), 4), "\n")
  cat("Number of non-zero connections:", length(non_zero), "\n")
  
  # Add attributes
  attr(similarity_matrix, "conversion_method") <- conversion_method
  attr(similarity_matrix, "normalization") <- normalization
  attr(similarity_matrix, "threshold_method") <- threshold_method
  attr(similarity_matrix, "threshold_value") <- threshold_value
  attr(similarity_matrix, "ready_for_qgraph") <- TRUE
  
  return(similarity_matrix)
}

#' Quick qgraph Plot with Automatic Preparation
#' 
#' Convenience function that prepares distance matrix and plots with qgraph
#' 
#' @param distance_matrix Distance matrix
#' @param layout qgraph layout
#' @param threshold qgraph threshold (applied after conversion/normalization)
#' @param conversion_method Distance-to-similarity conversion method
#' @param title Plot title
#' @param ... Additional arguments for qgraph
#' 
#' @return qgraph plot object (invisible)
#' @export
#' @examples
#' # Quick plot of large distance matrix
#' distances <- matrix(runif(16, 15000, 45000), nrow = 4)
#' distances <- (distances + t(distances)) / 2
#' diag(distances) <- 0
#' 
#' quick_qgraph(distances, layout = "spring", threshold = 0.8)
quick_qgraph <- function(distance_matrix,
                        layout = "spring", 
                        threshold = 0.7,
                        conversion_method = "normalized_inverse",
                        title = "Time Series Network",
                        ...) {
  
  if (!requireNamespace("qgraph", quietly = TRUE)) {
    stop("Package 'qgraph' is required. Install with: install.packages('qgraph')")
  }
  
  # Prepare matrix for visualization
  viz_matrix <- prepare_for_qgraph(
    distance_matrix, 
    conversion_method = conversion_method,
    normalization = "minmax"
  )
  
  # Create plot
  plot_result <- qgraph::qgraph(
    viz_matrix,
    layout = layout,
    threshold = threshold,
    title = title,
    ...
  )
  
  return(invisible(plot_result))
}

#' Show Conversion Examples
#'
#' Demonstrates the effect of different distance-to-similarity conversion methods
#' on an example distance matrix. This function is useful for understanding
#' how various methods transform distance values into similarity scores.
#'
#' @param distance_matrix A numeric distance matrix to use as an example.
#'   If `NULL` (default), a small example matrix will be generated.
#'
#' @return An invisible list of the converted similarity matrices, with attributes
#'   indicating the method and parameters used.
#' @export
#' @examples
#' # Show examples with a default generated distance matrix
#' show_conversion_examples()
#'
#' # Show examples with a custom distance matrix
#' my_distances <- matrix(c(0, 5, 15, 5, 0, 8, 15, 8, 0), nrow = 3, byrow = TRUE)
#' colnames(my_distances) <- rownames(my_distances) <- paste0("A", 1:3)
#' show_conversion_examples(my_distances)
show_conversion_examples <- function(distance_matrix = NULL) {
  
  # Create example if none provided
  if (is.null(distance_matrix)) {
    set.seed(42)
    distance_matrix <- matrix(runif(9, 15000, 45000), nrow = 3)
    distance_matrix <- (distance_matrix + t(distance_matrix)) / 2
    diag(distance_matrix) <- 0
    rownames(distance_matrix) <- colnames(distance_matrix) <- paste0("TS", 1:3)
  }
  
  cat("=== DISTANCE TO SIMILARITY CONVERSION EXAMPLES ===\n\n")
  
  cat("Original distance matrix:\n")
  print(round(distance_matrix, 0))
  cat("Range:", round(min(distance_matrix[distance_matrix > 0]), 0), "to", 
      round(max(distance_matrix), 0), "\n\n")
  
  methods <- c("normalized_inverse", "inverse", "negative_exp", "max_minus", "gaussian")
  results <- list()
  
  for (method in methods) {
    cat("Method:", method, "\n")
    
    sim_matrix <- distance_to_similarity(distance_matrix, method = method)
    results[[method]] <- sim_matrix
    
    print(round(sim_matrix, 4))
    
    non_zero <- sim_matrix[sim_matrix > 0]
    cat("Range:", round(min(non_zero), 4), "to", round(max(non_zero), 4))
    
    if (!is.null(attr(sim_matrix, "sigma"))) {
      cat(" (sigma =", round(attr(sim_matrix, "sigma"), 0), ")")
    }
    cat("\n\n")
  }
  
  return(invisible(results))
}