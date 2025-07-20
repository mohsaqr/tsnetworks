#' Resilience Plotting Functions for TSnetworks
#'
#' @description
#' Provides resilience-specific plotting functions that are fully compatible
#' with existing TSnetworks plotting infrastructure. These functions extend
#' the existing plotting capabilities to visualize resilience metrics and
#' states while maintaining consistency with TSnetworks aesthetics.
#'
#' @name resilience_plotting
NULL

#' Plot Resilience Capacity Timeline
#'
#' @description
#' Creates resilience capacity plots using TSnetworks plotting infrastructure.
#' This function extends the existing plot_timeseries_enhanced() functionality
#' to visualize resilience capacity metrics over time.
#'
#' @param data Data frame with resilience metrics (from calculate_resilience_metrics())
#' @param capacity_type Which capacity to plot: "all", "absorptive", "restorative", "adaptive", "composite"
#' @param ts_col Time series column name (auto-detected if NULL)
#' @param time_col Time column name (auto-detected if NULL)
#' @param state_col Resilience state column for shading (optional)
#' @param facet_by Column to facet by (optional)
#' @param smooth_line Whether to add smoothed trend line (default: TRUE)
#' @param color_palette Color palette to use (default: "viridis")
#' @param title Plot title (auto-generated if NULL)
#'
#' @return ggplot object compatible with existing TSnetworks plots
#'
#' @examples
#' # Basic resilience capacity timeline
#' data(saqrsteps)
#' resilience_data <- calculate_resilience_metrics(saqrsteps, ts_cols = "Steps")
#' state_data <- classify_resilience_states(resilience_data)
#' 
#' # Plot all capacities
#' plot_resilience_timeline(state_data, capacity_type = "all")
#' 
#' # Plot specific capacity with state shading
#' plot_resilience_timeline(state_data, capacity_type = "absorptive", 
#'                         state_col = "resilience_state")
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_smooth facet_wrap
#' @importFrom ggplot2 labs theme_minimal scale_color_viridis_d scale_fill_manual
#' @importFrom ggplot2 geom_rect alpha
#' @importFrom viridis scale_color_viridis
#' @export
plot_resilience_timeline <- function(data,
                                    capacity_type = c("all", "absorptive", "restorative", "adaptive", "composite"),
                                    ts_col = NULL,
                                    time_col = NULL,
                                    state_col = NULL,
                                    facet_by = NULL,
                                    smooth_line = TRUE,
                                    color_palette = "viridis",
                                    title = NULL) {
  
  capacity_type <- match.arg(capacity_type)
  
  # Validate inputs
  .validate_plotting_inputs(data, capacity_type)
  
  # Auto-detect columns
  if (is.null(ts_col)) {
    ts_col <- .detect_original_ts_column(data)
  }
  
  if (is.null(time_col)) {
    time_col <- .detect_time_column(data)
  }
  
  # Prepare data for plotting
  plot_data <- .prepare_resilience_plot_data(data, capacity_type, ts_col, time_col)
  
  # Create base plot
  p <- ggplot(plot_data, aes(x = time_index, y = value, color = capacity))
  
  # Add state shading if requested
  if (!is.null(state_col) && state_col %in% names(data)) {
    state_data <- .prepare_state_shading_data(data, state_col, time_col)
    p <- p + geom_rect(data = state_data, 
                       aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = state),
                       alpha = 0.2, inherit.aes = FALSE)
    
    # Add state colors
    state_colors <- .get_resilience_state_colors()
    p <- p + scale_fill_manual(values = state_colors, name = "Resilience State")
  }
  
  # Add main plot elements
  p <- p + geom_line(size = 0.8, alpha = 0.8)
  
  if (smooth_line) {
    p <- p + geom_smooth(method = "loess", se = TRUE, alpha = 0.3, size = 0.5)
  }
  
  # Apply color palette
  if (color_palette == "viridis") {
    p <- p + scale_color_viridis_d(name = "Capacity Type")
  } else {
    p <- p + scale_color_brewer(type = "qual", palette = color_palette, name = "Capacity Type")
  }
  
  # Add faceting if requested
  if (!is.null(facet_by) && facet_by %in% names(data)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_by)), scales = "free_y")
  }
  
  # Customize theme (consistent with TSnetworks)
  p <- p + theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      strip.text = element_text(size = 11, face = "bold")
    )
  
  # Add labels
  if (is.null(title)) {
    title <- paste("Resilience", stringr::str_to_title(capacity_type), "Capacity Timeline")
  }
  
  p <- p + labs(
    title = title,
    x = if (!is.null(time_col)) stringr::str_to_title(time_col) else "Time Index",
    y = "Capacity Score"
  )
  
  return(p)
}

#' Plot Resilience State Network
#'
#' @description
#' Creates resilience state transition networks using existing TSnetworks
#' plotting infrastructure. This function leverages plot_tna_network()
#' functionality to visualize resilience state transitions.
#'
#' @param data Data frame with resilience states (from classify_resilience_states())
#' @param state_col Resilience state column name (default: "resilience_state")
#' @param layout Network layout method (default: "spring")
#' @param node_size_by What to size nodes by: "frequency", "duration", "uniform"
#' @param edge_width_by What to size edges by: "transition_prob", "transition_count", "uniform"
#' @param show_labels Whether to show state labels (default: TRUE)
#' @param title Plot title (auto-generated if NULL)
#'
#' @return Network plot using existing TSnetworks infrastructure
#'
#' @examples
#' # Create resilience state network
#' data(saqrsteps)
#' resilience_data <- calculate_resilience_metrics(saqrsteps, ts_cols = "Steps")
#' state_data <- classify_resilience_states(resilience_data)
#' plot_resilience_network(state_data)
#'
#' # Customize network appearance
#' plot_resilience_network(state_data, node_size_by = "duration", 
#'                        edge_width_by = "transition_prob")
#'
#' @export
plot_resilience_network <- function(data,
                                   state_col = "resilience_state",
                                   layout = "spring",
                                   node_size_by = c("frequency", "duration", "uniform"),
                                   edge_width_by = c("transition_prob", "transition_count", "uniform"),
                                   show_labels = TRUE,
                                   title = NULL) {
  
  node_size_by <- match.arg(node_size_by)
  edge_width_by <- match.arg(edge_width_by)
  
  # Validate inputs
  if (!state_col %in% names(data)) {
    stop(paste("State column", state_col, "not found in data"), call. = FALSE)
  }
  
  # Create transition matrix
  transition_matrix <- create_resilience_transition_matrix(data, state_col, normalize = TRUE)
  
  # Get state statistics
  state_summary <- get_resilience_state_summary(data, state_col)
  
  # Prepare network data for plotting
  network_data <- .prepare_network_plot_data(transition_matrix, state_summary, 
                                            node_size_by, edge_width_by)
  
  # Use existing TSnetworks plotting if available, otherwise create custom plot
  if (exists("plot_tna_network", mode = "function")) {
    # Try to use existing TSnetworks network plotting
    tryCatch({
      # Create a data structure compatible with plot_tna_network
      network_compatible_data <- .convert_to_tna_format(data, state_col, transition_matrix)
      plot <- plot_tna_network(network_compatible_data, state_col = state_col)
      
      if (!is.null(title)) {
        plot <- plot + labs(title = title)
      }
      
      return(plot)
    }, error = function(e) {
      # Fall back to custom network plot
      return(.create_custom_network_plot(network_data, layout, show_labels, title))
    })
  } else {
    # Create custom network plot
    return(.create_custom_network_plot(network_data, layout, show_labels, title))
  }
}

#' Plot Resilience Capacity Radar Chart
#'
#' @description
#' Creates radar charts showing resilience capacity profiles over time
#' or for different groups. This provides a multi-dimensional view of
#' resilience performance.
#'
#' @param data Data frame with resilience metrics
#' @param group_by Column to group by (optional, for multiple radar charts)
#' @param time_points Specific time points to plot (optional)
#' @param aggregate_method How to aggregate data: "mean", "median", "latest"
#' @param title Plot title (auto-generated if NULL)
#'
#' @return Radar chart plot
#'
#' @examples
#' # Basic radar chart
#' data(saqrsteps)
#' resilience_data <- calculate_resilience_metrics(saqrsteps, ts_cols = "Steps")
#' state_data <- classify_resilience_states(resilience_data)
#' plot_resilience_radar(state_data)
#'
#' # Radar chart by resilience state
#' plot_resilience_radar(state_data, group_by = "resilience_state")
#'
#' @export
plot_resilience_radar <- function(data,
                                  group_by = NULL,
                                  time_points = NULL,
                                  aggregate_method = c("mean", "median", "latest"),
                                  title = NULL) {
  
  aggregate_method <- match.arg(aggregate_method)
  
  # Prepare radar chart data
  radar_data <- .prepare_radar_chart_data(data, group_by, time_points, aggregate_method)
  
  # Create radar chart
  radar_plot <- .create_radar_chart(radar_data, group_by, title)
  
  return(radar_plot)
}

#' Plot Combined Resilience Dashboard
#'
#' @description
#' Creates a comprehensive resilience visualization dashboard combining
#' multiple plot types. This function extends the existing TSnetworks
#' dashboard approach for resilience analysis.
#'
#' @param data Data frame with resilience analysis
#' @param ts_col Time series column name
#' @param time_col Time column name (optional)
#' @param state_col Resilience state column name
#' @param include_network Whether to include network plot (default: TRUE)
#' @param include_radar Whether to include radar chart (default: TRUE)
#' @param title Dashboard title (auto-generated if NULL)
#'
#' @return Multi-panel plot object
#'
#' @examples
#' # Create comprehensive resilience dashboard
#' data(saqrsteps)
#' resilience_data <- calculate_resilience_metrics(saqrsteps, ts_cols = "Steps")
#' state_data <- classify_resilience_states(resilience_data)
#' plot_resilience_dashboard(state_data, ts_col = "Steps", state_col = "resilience_state")
#'
#' @importFrom gridExtra grid.arrange
#' @export
plot_resilience_dashboard <- function(data,
                                     ts_col,
                                     time_col = NULL,
                                     state_col = "resilience_state",
                                     include_network = TRUE,
                                     include_radar = TRUE,
                                     title = NULL) {
  
  # Create individual plots
  plots <- list()
  
  # 1. Original time series with state shading
  if (exists("plot_timeseries_enhanced", mode = "function")) {
    tryCatch({
      plots$timeseries <- plot_timeseries_enhanced(data, ts_col = ts_col, 
                                                  state_col = state_col,
                                                  title = "Original Time Series")
    }, error = function(e) {
      plots$timeseries <- .create_basic_timeseries_plot(data, ts_col, state_col)
    })
  } else {
    plots$timeseries <- .create_basic_timeseries_plot(data, ts_col, state_col)
  }
  
  # 2. Resilience capacity timeline
  plots$capacity <- plot_resilience_timeline(data, capacity_type = "all", 
                                            ts_col = ts_col, time_col = time_col,
                                            state_col = state_col,
                                            title = "Resilience Capacities")
  
  # 3. Network plot (optional)
  if (include_network) {
    plots$network <- plot_resilience_network(data, state_col = state_col,
                                            title = "State Transition Network")
  }
  
  # 4. Radar chart (optional)
  if (include_radar) {
    plots$radar <- plot_resilience_radar(data, group_by = state_col,
                                        title = "Capacity Profile by State")
  }
  
  # Arrange plots
  if (length(plots) == 2) {
    # Just timeseries and capacity
    combined_plot <- gridExtra::grid.arrange(plots$timeseries, plots$capacity, 
                                           ncol = 1, heights = c(1, 1))
  } else if (length(plots) == 3) {
    # Three plots
    combined_plot <- gridExtra::grid.arrange(plots$timeseries, plots$capacity, 
                                           plots[[3]], ncol = 1, heights = c(1, 1, 1))
  } else if (length(plots) == 4) {
    # All four plots
    combined_plot <- gridExtra::grid.arrange(plots$timeseries, plots$capacity,
                                           plots$network, plots$radar,
                                           ncol = 2, nrow = 2)
  }
  
  # Add overall title if provided
  if (!is.null(title)) {
    combined_plot <- gridExtra::grid.arrange(
      gridExtra::textGrob(title, gp = grid::gpar(fontsize = 16, fontface = "bold")),
      combined_plot,
      heights = c(0.1, 0.9)
    )
  }
  
  return(combined_plot)
}

# Internal helper functions

.validate_plotting_inputs <- function(data, capacity_type) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame", call. = FALSE)
  }
  
  if (nrow(data) == 0) {
    stop("'data' cannot be empty", call. = FALSE)
  }
  
  # Check for resilience metrics
  resilience_cols <- grep("^(absorptive|restorative|adaptive)_", names(data), value = TRUE)
  if (length(resilience_cols) == 0) {
    stop("No resilience metrics found in data. Run calculate_resilience_metrics() first.", call. = FALSE)
  }
}

.detect_original_ts_column <- function(data) {
  # Look for resilience metric columns to infer original ts_cols
  resilience_cols <- grep("^(absorptive|restorative|adaptive)_", names(data), value = TRUE)
  if (length(resilience_cols) > 0) {
    # Extract original column names from resilience metric column names
    ts_cols <- unique(gsub("^(absorptive|restorative|adaptive)_[^_]+_", "", resilience_cols))
    return(ts_cols[1])  # Return first one
  } else {
    # Fallback to first numeric column
    numeric_cols <- names(data)[sapply(data, is.numeric)]
    return(numeric_cols[1])
  }
}

.detect_time_column <- function(data) {
  # Look for common time column names
  time_candidates <- c("time", "date", "timestamp", "index", "row_number")
  
  for (candidate in time_candidates) {
    if (candidate %in% names(data)) {
      return(candidate)
    }
  }
  
  # If no time column found, create index
  return(NULL)
}

.prepare_resilience_plot_data <- function(data, capacity_type, ts_col, time_col) {
  
  # Create time index if no time column
  if (is.null(time_col)) {
    data$time_index <- 1:nrow(data)
    time_col <- "time_index"
  }
  
  # Prepare data based on capacity type
  if (capacity_type == "all") {
    # Get all capacity score columns
    score_cols <- grep("^resilience_state_(absorptive|restorative|adaptive)_score$", 
                      names(data), value = TRUE)
    
    if (length(score_cols) == 0) {
      stop("No capacity score columns found. Run classify_resilience_states() first.", call. = FALSE)
    }
    
    # Reshape to long format
    plot_data <- data.frame()
    for (col in score_cols) {
      capacity_name <- gsub("^resilience_state_(.+)_score$", "\\1", col)
      temp_data <- data.frame(
        time_index = data[[time_col]],
        value = data[[col]],
        capacity = capacity_name
      )
      plot_data <- rbind(plot_data, temp_data)
    }
    
  } else if (capacity_type == "composite") {
    # Use composite score
    composite_col <- "resilience_state_composite_score"
    if (!composite_col %in% names(data)) {
      stop("Composite score column not found. Run classify_resilience_states() first.", call. = FALSE)
    }
    
    plot_data <- data.frame(
      time_index = data[[time_col]],
      value = data[[composite_col]],
      capacity = "composite"
    )
    
  } else {
    # Specific capacity type
    score_col <- paste0("resilience_state_", capacity_type, "_score")
    if (!score_col %in% names(data)) {
      stop(paste("Score column", score_col, "not found. Run classify_resilience_states() first."), call. = FALSE)
    }
    
    plot_data <- data.frame(
      time_index = data[[time_col]],
      value = data[[score_col]],
      capacity = capacity_type
    )
  }
  
  return(plot_data)
}

.prepare_state_shading_data <- function(data, state_col, time_col) {
  
  if (is.null(time_col)) {
    data$time_index <- 1:nrow(data)
    time_col <- "time_index"
  }
  
  states <- data[[state_col]]
  time_vals <- data[[time_col]]
  
  # Create run-length encoding for states
  state_runs <- rle(as.character(states))
  
  # Create shading data
  shading_data <- data.frame()
  current_pos <- 1
  
  for (i in seq_along(state_runs$lengths)) {
    if (!is.na(state_runs$values[i])) {
      start_pos <- current_pos
      end_pos <- current_pos + state_runs$lengths[i] - 1
      
      shading_data <- rbind(shading_data, data.frame(
        start = time_vals[start_pos],
        end = time_vals[end_pos],
        state = state_runs$values[i]
      ))
    }
    current_pos <- current_pos + state_runs$lengths[i]
  }
  
  return(shading_data)
}

.get_resilience_state_colors <- function() {
  state_defs <- .get_resilience_state_definitions()
  colors <- sapply(state_defs, function(x) x$color)
  return(colors)
}

.prepare_network_plot_data <- function(transition_matrix, state_summary, node_size_by, edge_width_by) {
  
  # Prepare node data
  nodes <- data.frame(
    id = rownames(transition_matrix),
    label = rownames(transition_matrix),
    stringsAsFactors = FALSE
  )
  
  # Add node sizes
  if (node_size_by == "frequency") {
    nodes$size <- sapply(nodes$id, function(x) state_summary$state_statistics[[x]]$count)
  } else if (node_size_by == "duration") {
    nodes$size <- sapply(nodes$id, function(x) state_summary$state_statistics[[x]]$duration_stats$mean_duration)
  } else {
    nodes$size <- rep(1, nrow(nodes))
  }
  
  # Add node colors
  state_colors <- .get_resilience_state_colors()
  nodes$color <- state_colors[nodes$id]
  
  # Prepare edge data
  edges <- data.frame()
  for (i in 1:nrow(transition_matrix)) {
    for (j in 1:ncol(transition_matrix)) {
      if (transition_matrix[i, j] > 0) {
        edges <- rbind(edges, data.frame(
          from = rownames(transition_matrix)[i],
          to = colnames(transition_matrix)[j],
          weight = transition_matrix[i, j],
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Add edge widths
  if (edge_width_by == "transition_prob") {
    edges$width <- edges$weight
  } else if (edge_width_by == "transition_count") {
    # Would need count matrix for this - use weight as proxy
    edges$width <- edges$weight
  } else {
    edges$width <- rep(1, nrow(edges))
  }
  
  return(list(nodes = nodes, edges = edges))
}

.convert_to_tna_format <- function(data, state_col, transition_matrix) {
  # Try to create a format compatible with existing plot_tna_network
  # This is a best-effort conversion
  
  result <- list(
    data = data,
    transition_matrix = transition_matrix,
    state_column = state_col
  )
  
  return(result)
}

.create_custom_network_plot <- function(network_data, layout, show_labels, title) {
  
  # Create a simple network plot using ggplot2
  # This is a fallback when TSnetworks plotting is not available
  
  nodes <- network_data$nodes
  edges <- network_data$edges
  
  # Simple circular layout
  n_nodes <- nrow(nodes)
  angles <- seq(0, 2*pi, length.out = n_nodes + 1)[1:n_nodes]
  nodes$x <- cos(angles)
  nodes$y <- sin(angles)
  
  # Add edge coordinates
  edges$x <- nodes$x[match(edges$from, nodes$id)]
  edges$y <- nodes$y[match(edges$from, nodes$id)]
  edges$xend <- nodes$x[match(edges$to, nodes$id)]
  edges$yend <- nodes$y[match(edges$to, nodes$id)]
  
  # Create plot
  p <- ggplot() +
    geom_segment(data = edges, aes(x = x, y = y, xend = xend, yend = yend, 
                                  size = width), alpha = 0.6, color = "gray50") +
    geom_point(data = nodes, aes(x = x, y = y, size = size, color = id), alpha = 0.8) +
    scale_size_continuous(range = c(3, 10), guide = "none") +
    scale_color_manual(values = nodes$color, name = "State") +
    theme_void() +
    theme(legend.position = "bottom")
  
  if (show_labels) {
    p <- p + geom_text(data = nodes, aes(x = x, y = y, label = label), 
                      size = 3, fontface = "bold")
  }
  
  if (!is.null(title)) {
    p <- p + labs(title = title)
  }
  
  return(p)
}

.prepare_radar_chart_data <- function(data, group_by, time_points, aggregate_method) {
  
  # Get capacity score columns
  score_cols <- grep("^resilience_state_(absorptive|restorative|adaptive)_score$", 
                    names(data), value = TRUE)
  
  if (length(score_cols) == 0) {
    stop("No capacity score columns found. Run classify_resilience_states() first.", call. = FALSE)
  }
  
  # Prepare aggregation
  if (!is.null(group_by) && group_by %in% names(data)) {
    # Group by specified column
    groups <- unique(data[[group_by]])
    radar_data <- data.frame()
    
    for (group in groups) {
      group_data <- data[data[[group_by]] == group, ]
      
      group_scores <- sapply(score_cols, function(col) {
        values <- group_data[[col]]
        if (aggregate_method == "mean") {
          mean(values, na.rm = TRUE)
        } else if (aggregate_method == "median") {
          median(values, na.rm = TRUE)
        } else {
          tail(values[!is.na(values)], 1)
        }
      })
      
      radar_data <- rbind(radar_data, data.frame(
        group = group,
        capacity = gsub("^resilience_state_(.+)_score$", "\\1", names(group_scores)),
        value = as.numeric(group_scores)
      ))
    }
    
  } else {
    # Overall aggregation
    overall_scores <- sapply(score_cols, function(col) {
      values <- data[[col]]
      if (aggregate_method == "mean") {
        mean(values, na.rm = TRUE)
      } else if (aggregate_method == "median") {
        median(values, na.rm = TRUE)
      } else {
        tail(values[!is.na(values)], 1)
      }
    })
    
    radar_data <- data.frame(
      group = "Overall",
      capacity = gsub("^resilience_state_(.+)_score$", "\\1", names(overall_scores)),
      value = as.numeric(overall_scores)
    )
  }
  
  return(radar_data)
}

.create_radar_chart <- function(radar_data, group_by, title) {
  
  # Create a simple radar-like chart using ggplot2
  # For a true radar chart, would need additional packages like ggradar
  
  p <- ggplot(radar_data, aes(x = capacity, y = value)) +
    geom_col(aes(fill = capacity), alpha = 0.7) +
    coord_polar() +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 8),
      legend.position = "none"
    ) +
    labs(x = "", y = "Capacity Score")
  
  if (!is.null(group_by) && length(unique(radar_data$group)) > 1) {
    p <- p + facet_wrap(~ group)
  }
  
  if (!is.null(title)) {
    p <- p + labs(title = title)
  }
  
  return(p)
}

.create_basic_timeseries_plot <- function(data, ts_col, state_col) {
  
  # Create basic time series plot as fallback
  data$time_index <- 1:nrow(data)
  
  p <- ggplot(data, aes(x = time_index, y = !!sym(ts_col))) +
    geom_line(color = "steelblue", size = 0.8) +
    theme_minimal() +
    labs(x = "Time Index", y = ts_col, title = "Time Series")
  
  # Add state shading if available
  if (!is.null(state_col) && state_col %in% names(data)) {
    state_data <- .prepare_state_shading_data(data, state_col, "time_index")
    state_colors <- .get_resilience_state_colors()
    
    p <- p + geom_rect(data = state_data, 
                       aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = state),
                       alpha = 0.2, inherit.aes = FALSE) +
      scale_fill_manual(values = state_colors, name = "Resilience State")
  }
  
  return(p)
}
