#' Create TNA Network Plot for Any Classification/Discretization Results
#'
#' @description
#' A unified function to create Temporal Network Analysis (TNA) plots for any
#' classification or discretization results from the tsnetworks package.
#' This function works with results from stna(), detect_regime(), calculate_hurst(),
#' visibility_graph(), or any data frame with state/classification columns.
#'
#' @param data Data frame containing classification/state results
#' @param state_col Column name containing the states/classifications (auto-detected if NULL)
#' @param time_col Column name for time axis (optional, uses row indices if NULL)
#' @param layout Layout algorithm for network: "spring", "circle", "tree", "grid"
#' @param node_size_method Method for node sizing: "frequency", "duration", "uniform"
#' @param edge_weight_method Method for edge weighting: "frequency", "probability", "uniform"
#' @param color_palette Color palette for states: "default", "viridis", "plasma", "rainbow"
#' @param show_labels Whether to show state labels on nodes
#' @param show_edge_labels Whether to show transition probabilities on edges
#' @param min_edge_weight Minimum edge weight to display (filters weak connections)
#' @param title Plot title (auto-generated if NULL)
#' @param ... Additional arguments passed to plotting functions
#'
#' @return Network plot (invisible return of network object)
#'
#' @examples
#' # Example with STNA results
#' data(saqrsteps)
#' stna_result <- stna(saqrsteps, value_column = "Steps", num_states = 4, method = "quantile")
#' plot_tna_network(stna_result)
#'
#' # Example with regime detection results
#' complexity_data <- rolling_measures(saqrsteps, "Steps", window_width = 7)
#' regime_result <- detect_regime(complexity_data, method = "smart")
#' plot_tna_network(regime_result, state_col = "regime_stability")
#'
#' # Example with Hurst analysis results
#' hurst_result <- calculate_hurst(saqrsteps, ts_col = "Steps", method = "dfa")
#' plot_tna_network(hurst_result, state_col = "dfa.state")
#'
#' @importFrom graphics plot text arrows
#' @importFrom grDevices rainbow
#' @export
plot_tna_network <- function(data,
                             state_col = NULL,
                             time_col = NULL,
                             layout = c("spring", "circle", "tree", "grid"),
                             node_size_method = c("frequency", "duration", "uniform"),
                             edge_weight_method = c("frequency", "probability", "uniform"),
                             color_palette = c("default", "viridis", "plasma", "rainbow"),
                             show_labels = TRUE,
                             show_edge_labels = FALSE,
                             min_edge_weight = 0.01,
                             title = NULL,
                             ...) {

  # Match arguments
  layout <- match.arg(layout)
  node_size_method <- match.arg(node_size_method)
  edge_weight_method <- match.arg(edge_weight_method)
  color_palette <- match.arg(color_palette)

  # Auto-detect state column if not provided
  if (is.null(state_col)) {
    # Look for common state column patterns
    state_patterns <- c("\\.state$", "^state", "^regime", "^classification", "^cluster")
    
    for (pattern in state_patterns) {
      matches <- grep(pattern, names(data), value = TRUE)
      if (length(matches) > 0) {
        state_col <- matches[1]
        message("Auto-detected state column: ", state_col)
        break
      }
    }
    
    if (is.null(state_col)) {
      stop("Could not auto-detect state column. Please specify 'state_col' parameter.")
    }
  }

  # Validate state column exists
  if (!state_col %in% names(data)) {
    stop(paste("State column", state_col, "not found in data"))
  }

  # Extract states and clean data
  states <- as.character(data[[state_col]])
  states <- states[!is.na(states)]
  
  if (length(states) == 0) {
    stop("No valid states found in the specified column")
  }

  # Create transition matrix
  transitions <- table(
    From = states[-length(states)],
    To = states[-1]
  )

  # Convert to probability matrix if requested
  if (edge_weight_method == "probability") {
    edge_weights <- prop.table(transitions, margin = 1)
  } else if (edge_weight_method == "frequency") {
    edge_weights <- transitions
  } else {
    edge_weights <- (transitions > 0) * 1
  }

  # Filter weak edges
  edge_weights[edge_weights < min_edge_weight] <- 0

  # Calculate node properties
  unique_states <- rownames(transitions)
  n_states <- length(unique_states)

  # Node sizes based on method
  if (node_size_method == "frequency") {
    state_counts <- table(states)
    node_sizes <- as.numeric(state_counts[unique_states])
    node_sizes <- (node_sizes / max(node_sizes)) * 3 + 1  # Scale to 1-4
  } else if (node_size_method == "duration") {
    # Calculate average duration in each state
    state_runs <- rle(states)
    avg_durations <- tapply(state_runs$lengths, state_runs$values, mean)
    node_sizes <- as.numeric(avg_durations[unique_states])
    node_sizes[is.na(node_sizes)] <- 1
    node_sizes <- (node_sizes / max(node_sizes, na.rm = TRUE)) * 3 + 1
  } else {
    node_sizes <- rep(2, n_states)
  }

  # Node colors based on palette
  if (color_palette == "default") {
    node_colors <- c("#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEAA7", 
                     "#DDA0DD", "#98D8C8", "#F7DC6F", "#BB8FCE", "#85C1E9")[1:n_states]
  } else if (color_palette == "viridis") {
    if (requireNamespace("viridis", quietly = TRUE)) {
      node_colors <- viridis::viridis(n_states)
    } else {
      node_colors <- rainbow(n_states)
    }
  } else if (color_palette == "plasma") {
    if (requireNamespace("viridis", quietly = TRUE)) {
      node_colors <- viridis::plasma(n_states)
    } else {
      node_colors <- rainbow(n_states)
    }
  } else {
    node_colors <- rainbow(n_states)
  }

  # Create layout coordinates
  if (layout == "circle") {
    angles <- seq(0, 2*pi, length.out = n_states + 1)[1:n_states]
    x_coords <- cos(angles)
    y_coords <- sin(angles)
  } else if (layout == "grid") {
    grid_size <- ceiling(sqrt(n_states))
    x_coords <- rep(1:grid_size, length.out = n_states)
    y_coords <- rep(1:grid_size, each = grid_size)[1:n_states]
  } else if (layout == "tree") {
    # Simple tree layout
    levels <- ceiling(log2(n_states + 1))
    x_coords <- rep(1:ceiling(n_states/levels), length.out = n_states)
    y_coords <- rep(1:levels, each = ceiling(n_states/levels))[1:n_states]
  } else {
    # Spring layout (simplified force-directed)
    set.seed(123)  # For reproducibility
    x_coords <- runif(n_states, -1, 1)
    y_coords <- runif(n_states, -1, 1)
    
    # Simple force-directed adjustment
    for (iter in 1:50) {
      for (i in 1:n_states) {
        fx <- fy <- 0
        for (j in 1:n_states) {
          if (i != j) {
            dx <- x_coords[i] - x_coords[j]
            dy <- y_coords[i] - y_coords[j]
            dist <- sqrt(dx^2 + dy^2)
            if (dist > 0) {
              # Repulsion
              fx <- fx + dx / (dist^2) * 0.1
              fy <- fy + dy / (dist^2) * 0.1
              # Attraction for connected nodes
              if (edge_weights[i, j] > 0) {
                fx <- fx - dx * 0.01
                fy <- fy - dy * 0.01
              }
            }
          }
        }
        x_coords[i] <- x_coords[i] + fx * 0.1
        y_coords[i] <- y_coords[i] + fy * 0.1
      }
    }
  }

  # Create the plot
  if (is.null(title)) {
    title <- paste("State Transition Network -", state_col)
  }

  plot(x_coords, y_coords, 
       type = "n",
       xlim = range(x_coords) + c(-0.5, 0.5),
       ylim = range(y_coords) + c(-0.5, 0.5),
       main = title,
       xlab = "", ylab = "",
       axes = FALSE,
       ...)

  # Draw edges
  for (i in 1:n_states) {
    for (j in 1:n_states) {
      if (edge_weights[i, j] > 0) {
        # Calculate edge properties
        edge_width <- sqrt(edge_weights[i, j] / max(edge_weights)) * 5
        edge_alpha <- min(1, edge_weights[i, j] / max(edge_weights) + 0.3)
        
        # Draw arrow
        arrows(x_coords[i], y_coords[i],
               x_coords[j], y_coords[j],
               length = 0.1,
               lwd = edge_width,
               col = adjustcolor("gray30", alpha = edge_alpha))
        
        # Add edge labels if requested
        if (show_edge_labels && edge_weights[i, j] > min_edge_weight * 2) {
          mid_x <- (x_coords[i] + x_coords[j]) / 2
          mid_y <- (y_coords[i] + y_coords[j]) / 2
          text(mid_x, mid_y, 
               sprintf("%.2f", edge_weights[i, j]),
               cex = 0.6, col = "red")
        }
      }
    }
  }

  # Draw nodes
  points(x_coords, y_coords,
         pch = 21,
         cex = node_sizes,
         bg = node_colors,
         col = "black",
         lwd = 2)

  # Add node labels
  if (show_labels) {
    text(x_coords, y_coords,
         labels = unique_states,
         cex = 0.8,
         col = "black",
         font = 2)
  }

  # Add legend
  legend("topright",
         legend = unique_states,
         fill = node_colors,
         cex = 0.8,
         title = "States")

  # Return network information invisibly
  network_info <- list(
    states = unique_states,
    transitions = transitions,
    edge_weights = edge_weights,
    coordinates = data.frame(x = x_coords, y = y_coords, state = unique_states),
    node_sizes = node_sizes,
    node_colors = node_colors
  )
  
  invisible(network_info)
}

#' Create Combined Time Series and Network Plot
#'
#' @description
#' Creates a comprehensive visualization showing both the time series with
#' state classifications and the corresponding state transition network.
#'
#' @param data Data frame with time series and state data
#' @param ts_col Column name for time series values
#' @param state_col Column name for states/classifications
#' @param time_col Column name for time axis (optional)
#' @param network_layout Layout for network plot
#' @param ... Additional arguments passed to plotting functions
#'
#' @return Combined plot (invisible return of both plot objects)
#'
#' @examples
#' data(saqrsteps)
#' result <- stna(saqrsteps, value_column = "Steps", num_states = 4, method = "quantile")
#' plot_combined_tna(result, ts_col = "Steps", state_col = "state")
#'
#' @export
plot_combined_tna <- function(data, ts_col, state_col = NULL, time_col = NULL,
                              network_layout = "spring", ...) {
  
  # Set up layout for two plots
  par(mfrow = c(1, 2))
  
  # Plot 1: Time series with states
  plot_timeseries(data, value_col = ts_col, time_col = time_col, 
                  state_col = state_col, ...)
  
  # Plot 2: Network
  network_info <- plot_tna_network(data, state_col = state_col, 
                                   layout = network_layout, ...)
  
  # Reset layout
  par(mfrow = c(1, 1))
  
  invisible(list(timeseries = data, network = network_info))
}

#' Enhanced Time Series Plotting with State Backgrounds
#'
#' @description
#' An enhanced version of the existing plot_timeseries function with better
#' state visualization and more customization options.
#'
#' @param data Data frame containing the time series
#' @param ts_col Column name for time series values
#' @param time_col Column name for time axis (optional)
#' @param state_col Column name for states/classifications (optional)
#' @param color_palette Color palette for states
#' @param alpha Transparency for state backgrounds
#' @param show_transitions Whether to mark state transitions
#' @param main Plot title
#' @param ... Additional arguments passed to plot()
#'
#' @return Plot (invisible return of state information)
#'
#' @examples
#' data(saqrsteps)
#' result <- stna(saqrsteps, value_column = "Steps", num_states = 4)
#' plot_timeseries_enhanced(result, ts_col = "Steps", state_col = "state")
#'
#' @importFrom graphics plot lines rect abline
#' @export
plot_timeseries_enhanced <- function(data, ts_col, time_col = NULL, state_col = NULL,
                                     color_palette = "default", alpha = 0.3,
                                     show_transitions = TRUE, main = NULL, ...) {
  
  # Determine time axis
  if (!is.null(time_col) && time_col %in% names(data)) {
    time_axis <- data[[time_col]]
    xlab <- time_col
  } else {
    time_axis <- 1:nrow(data)
    xlab <- "Index"
  }
  
  # Set default title
  if (is.null(main)) {
    main <- paste("Time Series:", ts_col)
    if (!is.null(state_col)) {
      main <- paste(main, "with", state_col)
    }
  }
  
  # Create base plot
  plot(time_axis, data[[ts_col]], type = "l", lwd = 2,
       main = main, xlab = xlab, ylab = ts_col, ...)
  
  # Add state coloring if specified
  if (!is.null(state_col) && state_col %in% names(data)) {
    
    # Get state colors
    unique_states <- unique(data[[state_col]])
    unique_states <- unique_states[!is.na(unique_states)]
    
    if (color_palette == "default") {
      colors <- c("#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEAA7", 
                  "#DDA0DD", "#98D8C8", "#F7DC6F", "#BB8FCE", "#85C1E9")
    } else {
      colors <- rainbow(length(unique_states))
    }
    
    names(colors) <- unique_states
    
    # Find state change points
    states <- data[[state_col]]
    state_numeric <- as.numeric(factor(states))
    change_points <- which(diff(state_numeric) != 0)
    change_points <- c(0, change_points, nrow(data))
    
    # Add colored regions
    for (i in 1:(length(change_points) - 1)) {
      start <- change_points[i] + 1
      end <- change_points[i + 1]
      state <- as.character(states[start])
      
      if (!is.na(state) && state %in% names(colors)) {
        rect(time_axis[start], par("usr")[3], time_axis[end], par("usr")[4],
             col = adjustcolor(colors[state], alpha = alpha),
             border = NA)
      }
    }
    
    # Redraw the time series
    lines(time_axis, data[[ts_col]], col = "black", lwd = 2)
    
    # Mark transitions if requested
    if (show_transitions && length(change_points) > 2) {
      transition_points <- change_points[2:(length(change_points)-1)]
      abline(v = time_axis[transition_points], col = "red", lty = 2, alpha = 0.7)
    }
    
    # Add legend
    legend("topright",
           legend = unique_states,
           fill = adjustcolor(colors[unique_states], alpha = alpha),
           cex = 0.8,
           title = "States")
    
    # Return state information
    state_info <- list(
      states = unique_states,
      colors = colors,
      change_points = change_points,
      n_transitions = length(change_points) - 2
    )
    
    invisible(state_info)
  }
}
