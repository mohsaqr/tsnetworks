#' Plot Time Series with State Mapping
#'
#' A robust, self-contained function to visualize time series data with state overlays
#'
#' @param data A data.frame containing time series data
#' @param id_col Character name of the ID column (optional, creates single series if NULL)
#' @param value_col Character name of the value column  
#' @param time_col Character name of the time column (optional, uses row order if NULL)
#' @param state_col Character name of the state column
#' @param selected Character vector of specific series IDs to plot
#' @param overlay Character string: "v" (vertical), "h" (horizontal), or NULL (no overlay)
#' @param points Logical, whether to show points colored by state
#' @param ncol Integer, number of columns for facets
#' @param max_series Integer, maximum number of series to plot
#' @param trend Logical, whether to add trend line
#' @param scales Character, facet scales: "fixed", "free", "free_x", "free_y"
#' @return A ggplot object
#' @export

plot_timeseries <- function(data, 
                            id_col = NULL, 
                            value_col = "value", 
                            time_col = NULL, 
                            state_col = "state",
                            selected = NULL, 
                            overlay = "v", 
                            points = FALSE,
                            ncol = NULL, 
                            max_series = 10, 
                            trend = FALSE,
                            scales = "free") {
  
  # Input validation
  if (!is.data.frame(data)) stop("data must be a data.frame")
  if (nrow(data) == 0) stop("data cannot be empty")
  
  # Create ID column if not provided
  if (is.null(id_col)) {
    id_col <- ".series_id"
    data[[id_col]] <- "Series_1"
  }
  
  # Create time column if not provided
  if (is.null(time_col)) {
    time_col <- ".time_index"
    # Ensure data is sorted by id first to maintain proper order
    data <- data[order(data[[id_col]]), ]
    # Create sequential time index within each series
    data[[time_col]] <- unlist(lapply(split(seq_len(nrow(data)), data[[id_col]]), seq_along))
  }
  
  required_cols <- c(value_col, state_col)
  # Only check for id_col if it was originally provided
  if (!id_col %in% c(".series_id") && !is.null(id_col)) {
    required_cols <- c(required_cols, id_col)
  }
  # Only check for time_col if it was originally provided  
  if (!time_col %in% c(".time_index") && !is.null(time_col)) {
    required_cols <- c(required_cols, time_col)
  }
  
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (!is.null(overlay) && !overlay %in% c("h", "v")) {
    stop("overlay must be 'h', 'v', or NULL")
  }
  
  if (!scales %in% c("fixed", "free", "free_x", "free_y")) {
    stop("scales must be one of: 'fixed', 'free', 'free_x', 'free_y'")
  }
  
  # Select series to plot
  all_ids <- unique(data[[id_col]])
  if (is.null(selected)) {
    selected <- all_ids[seq_len(min(length(all_ids), max_series))]
  } else {
    selected <- selected[seq_len(min(length(selected), max_series))]
    selected <- intersect(selected, all_ids)  # Ensure selected IDs exist
    if (length(selected) == 0) {
      warning("No valid series IDs found in selected, using all available series")
      selected <- all_ids[seq_len(min(length(all_ids), max_series))]
    }
  }
  
  # Filter data
  plot_data <- data[data[[id_col]] %in% selected, ]
  if (nrow(plot_data) == 0) {
    stop("No data found for selected series")
  }
  plot_data <- plot_data[order(plot_data[[id_col]], plot_data[[time_col]]), ]
  
  # Create state overlay rectangles
  overlay_data <- NULL
  if (!is.null(overlay) && nrow(plot_data) > 1) {
    tryCatch({
      segment_col <- if (overlay == "v") time_col else value_col
      
      overlay_data <- plot_data %>%
        dplyr::arrange(!!dplyr::sym(id_col), !!dplyr::sym(segment_col)) %>%
        dplyr::group_by(!!dplyr::sym(id_col)) %>%
        dplyr::mutate(
          state_change = !!dplyr::sym(state_col) != dplyr::lag(!!dplyr::sym(state_col), default = dplyr::first(!!dplyr::sym(state_col))),
          group_id = cumsum(state_change),
          segment_lag = dplyr::lag(!!dplyr::sym(segment_col), default = dplyr::first(!!dplyr::sym(segment_col))),
          segment_lead = dplyr::lead(!!dplyr::sym(segment_col), default = dplyr::last(!!dplyr::sym(segment_col)))
        ) %>%
        dplyr::group_by(!!dplyr::sym(id_col), group_id, !!dplyr::sym(state_col)) %>%
        dplyr::summarise(
          segment_min = 0.5 * (min(!!dplyr::sym(segment_col)) + min(segment_lag, na.rm = TRUE)),
          segment_max = 0.5 * (max(!!dplyr::sym(segment_col)) + max(segment_lead, na.rm = TRUE)),
          .groups = "drop"
        )
      
      # Set rectangle coordinates based on overlay type
      if (overlay == "v") {
        overlay_data$xmin <- overlay_data$segment_min
        overlay_data$xmax <- overlay_data$segment_max
        overlay_data$ymin <- -Inf
        overlay_data$ymax <- Inf
      } else {  # overlay == "h"
        overlay_data$xmin <- -Inf
        overlay_data$xmax <- Inf
        overlay_data$ymin <- overlay_data$segment_min
        overlay_data$ymax <- overlay_data$segment_max
      }
    }, error = function(e) {
      warning("Could not create overlay rectangles: ", e$message)
      overlay_data <<- NULL
    })
  }
  
  # Color palette function
  get_state_colors <- function(n_states) {
    if (n_states <= 8) {
      # Accent palette colors
      c("#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0", "#f0027f", "#bf5b17", "#666666")[1:n_states]
    } else {
      # Set3 palette colors (expanded)
      grDevices::rainbow(n_states, alpha = 0.7)
    }
  }
  
  # Create base plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = !!dplyr::sym(time_col), y = !!dplyr::sym(value_col)))
  
  # Add state overlay rectangles
  if (!is.null(overlay_data) && nrow(overlay_data) > 0) {
    p <- p + ggplot2::geom_rect(
      data = overlay_data,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = !!dplyr::sym(state_col)),
      alpha = 0.5,
      inherit.aes = FALSE
    )
  }
  
  # Add main line
  p <- p + ggplot2::geom_line(linewidth = 0.5)
  
  # Add points if requested
  if (points) {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(fill = !!dplyr::sym(state_col)),
      shape = 21,
      size = 1.5,
      show.legend = FALSE
    )
  }
  
  # Add trend line if requested
  if (trend) {
    p <- p + ggplot2::geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      color = "darkgray",
      linetype = "dashed",
      linewidth = 0.8
    )
  }
  
  # Add facets for multiple series
  if (length(selected) > 1) {
    p <- p + ggplot2::facet_wrap(
      dplyr::vars(!!dplyr::sym(id_col)),
      ncol = ncol,
      scales = scales
    )
  }
  
  # Customize colors and theme
  unique_states <- sort(unique(plot_data[[state_col]]))
  state_colors <- get_state_colors(length(unique_states))
  names(state_colors) <- unique_states
  
  # Create appropriate x-axis label
  x_label <- if (is.null(time_col) || time_col == ".time_index") "Index" else "Time"
  series_count <- length(unique(plot_data[[id_col]]))
  title_text <- if (id_col == ".series_id") "Time Series Plot" else paste("Time Series Plot -", series_count, "Series")
  
  p <- p +
    ggplot2::scale_fill_manual(
      values = state_colors,
      name = "State",
      breaks = unique_states
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = x_label,
      y = "Value",
      title = title_text
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  
  return(p)
}

#' Helper function to add time column if missing
#'
#' This internal helper function adds a sequential time column to a data frame
#' if it's missing, grouping by an ID column if specified. This is useful for
#' plotting functions that require a time index.
#'
#' @param data A data.frame.
#' @param id_col Character string. The name of the ID column to group by when
#'   creating the time index. Defaults to "id".
#' @param time_col Character string. The name of the time column to create.
#'   Defaults to "time".
#' @return The data.frame with the time column added if it was missing.
#' @importFrom dplyr group_by mutate ungroup arrange sym lag first lead summarise row_number vars
#' @importFrom ggplot2 ggplot aes geom_rect geom_line geom_point geom_smooth facet_wrap scale_fill_manual theme_minimal labs theme element_blank element_text
#' @importFrom utils globalVariables
#' @importFrom grDevices rainbow
#' @name add_time_column
#' @keywords internal

utils::globalVariables(c("state_change", "group_id", "segment_lag", "segment_lead", "xmin", "xmax", "ymin", "ymax"))

add_time_column <- function(data, id_col = "id", time_col = "time") {
  if (time_col %in% names(data)) {
    return(data)
  }
  
  # Add sequential time column for each series
  data %>%
    dplyr::group_by(!!dplyr::sym(id_col)) %>%
    dplyr::mutate(!!time_col := dplyr::row_number()) %>%
    dplyr::ungroup()
}
