#' Integration Functions for TSnetworks Resilience Analysis
#'
#' @description
#' Provides seamless integration functions that allow resilience analysis
#' to be added to existing TSnetworks workflows without modifying any
#' existing functions. These functions act as bridges between resilience
#' analysis and existing TSnetworks capabilities.
#'
#' @name resilience_integration
NULL

#' Add Resilience Analysis to Existing TSnetworks Results
#'
#' @description
#' Takes results from existing TSnetworks functions and adds comprehensive
#' resilience analysis. This function automatically detects the type of
#' TSnetworks result and applies appropriate resilience analysis while
#' maintaining all existing functionality.
#'
#' @param tsnetworks_result Result from stna(), rolling_measures(), or other TSnetworks functions
#' @param capacity_types Which resilience capacities to analyze
#' @param window_width Rolling window size for resilience metrics
#' @param scaling_method Scaling method for resilience metrics
#' @param classify_states Whether to classify resilience states (default: TRUE)
#' @param state_method State classification method
#' @param ... Additional arguments passed to resilience functions
#'
#' @return Enhanced TSnetworks object with resilience analysis added
#'
#' @examples
#' # Add resilience to existing STNA analysis
#' data(saqrsteps)
#' stna_result <- stna(saqrsteps, "Steps", num_states = 4, method = "quantile")
#' enhanced_result <- add_resilience_to_tsnetworks(stna_result)
#' 
#' # Use enhanced result with existing TSnetworks functions
#' plot_tna_network(enhanced_result$data, state_col = "resilience_state")
#'
#' # Add resilience to rolling measures analysis
#' rolling_result <- rolling_measures(saqrsteps, ts_cols = "Steps", window_width = 7)
#' enhanced_rolling <- add_resilience_to_tsnetworks(rolling_result)
#'
#' @export
add_resilience_to_tsnetworks <- function(tsnetworks_result,
                                        capacity_types = c("absorptive", "restorative", "adaptive"),
                                        window_width = 7L,
                                        scaling_method = "auto",
                                        classify_states = TRUE,
                                        state_method = "threshold",
                                        ...) {
  
  # Detect the type of TSnetworks result
  result_type <- .detect_tsnetworks_result_type(tsnetworks_result)
  
  # Extract data component
  if (result_type == "stna") {
    data <- tsnetworks_result$data
    ts_cols <- .extract_ts_columns_from_stna(tsnetworks_result)
  } else if (result_type == "rolling_measures") {
    data <- tsnetworks_result
    ts_cols <- .extract_ts_columns_from_rolling(tsnetworks_result)
  } else if (result_type == "data_frame") {
    data <- tsnetworks_result
    ts_cols <- .auto_detect_ts_columns(tsnetworks_result)
  } else {
    stop("Unsupported TSnetworks result type. Please provide stna() result, rolling_measures() result, or data frame.", call. = FALSE)
  }
  
  # Apply resilience analysis
  resilience_data <- calculate_resilience_metrics(
    data = data,
    ts_cols = ts_cols,
    window_width = window_width,
    scaling_method = scaling_method,
    capacity_types = capacity_types,
    ...
  )
  
  # Classify states if requested
  if (classify_states) {
    resilience_data <- classify_resilience_states(
      data = resilience_data,
      method = state_method,
      ts_cols = ts_cols,
      ...
    )
  }
  
  # Integrate with original result structure
  if (result_type == "stna") {
    # Enhance STNA result
    enhanced_result <- tsnetworks_result
    enhanced_result$data <- resilience_data
    enhanced_result$resilience_analysis <- list(
      capacity_types = capacity_types,
      window_width = window_width,
      scaling_method = scaling_method,
      state_classification = classify_states,
      state_method = if (classify_states) state_method else NULL
    )
    
    # Add resilience transition matrix if states were classified
    if (classify_states) {
      enhanced_result$resilience_transition_matrix <- create_resilience_transition_matrix(
        resilience_data, state_col = "resilience_state"
      )
      enhanced_result$resilience_summary <- get_resilience_state_summary(
        resilience_data, state_col = "resilience_state", ts_cols = ts_cols
      )
    }
    
    class(enhanced_result) <- c("enhanced_stna", class(enhanced_result))
    
  } else if (result_type == "rolling_measures") {
    # Enhance rolling measures result
    enhanced_result <- resilience_data
    attr(enhanced_result, "original_rolling_measures") <- tsnetworks_result
    attr(enhanced_result, "resilience_enhancement") <- list(
      capacity_types = capacity_types,
      window_width = window_width,
      scaling_method = scaling_method,
      state_classification = classify_states
    )
    
    class(enhanced_result) <- c("enhanced_rolling_measures", class(enhanced_result))
    
  } else {
    # Data frame result
    enhanced_result <- resilience_data
    attr(enhanced_result, "resilience_enhancement") <- list(
      capacity_types = capacity_types,
      window_width = window_width,
      scaling_method = scaling_method,
      state_classification = classify_states
    )
    
    class(enhanced_result) <- c("enhanced_tsnetworks_data", class(enhanced_result))
  }
  
  return(enhanced_result)
}

#' Complete Resilience Analysis Pipeline
#'
#' @description
#' Performs comprehensive resilience analysis that integrates seamlessly
#' with TSnetworks workflow. This function can optionally include existing
#' TSnetworks analyses (STNA, rolling measures, Hurst analysis) alongside
#' resilience analysis.
#'
#' @param data Input data (data frame with time series)
#' @param ts_cols Time series column names
#' @param analysis_types Types of analysis to perform
#' @param resilience_config Configuration for resilience analysis
#' @param stna_config Configuration for STNA analysis (if included)
#' @param rolling_config Configuration for rolling measures (if included)
#' @param hurst_config Configuration for Hurst analysis (if included)
#' @param integrate_results Whether to create integrated visualizations
#'
#' @return Comprehensive analysis object with all requested analyses
#'
#' @examples
#' # Complete analysis pipeline
#' data(saqrsteps)
#' comprehensive_result <- resilience_analysis_pipeline(
#'   data = saqrsteps,
#'   ts_cols = "Steps",
#'   analysis_types = c("resilience", "stna", "rolling"),
#'   integrate_results = TRUE
#' )
#'
#' # Access different analysis components
#' print(comprehensive_result$resilience_summary)
#' print(comprehensive_result$stna_summary)
#' 
#' # Create integrated dashboard
#' plot_integrated_dashboard(comprehensive_result)
#'
#' @export
resilience_analysis_pipeline <- function(data,
                                        ts_cols,
                                        analysis_types = c("resilience", "stna", "rolling"),
                                        resilience_config = list(),
                                        stna_config = list(),
                                        rolling_config = list(),
                                        hurst_config = list(),
                                        integrate_results = TRUE) {
  
  # Validate inputs
  .validate_pipeline_inputs(data, ts_cols, analysis_types)
  
  # Initialize results container
  pipeline_results <- list(
    data = data,
    ts_cols = ts_cols,
    analysis_types = analysis_types,
    timestamp = Sys.time()
  )
  
  # Apply scaling if requested
  scaling_method <- resilience_config$scaling_method %||% "auto"
  if (scaling_method != "none") {
    if (scaling_method == "auto") {
      scaling_method <- auto_detect_scaling(data[[ts_cols[1]]])
    }
    scaled_data <- scale_dataframe(data, cols = ts_cols, method = scaling_method)
    pipeline_results$scaled_data <- scaled_data
    working_data <- scaled_data
  } else {
    working_data <- data
  }
  
  # Perform resilience analysis
  if ("resilience" %in% analysis_types) {
    resilience_args <- list(
      data = working_data,
      ts_cols = ts_cols,
      window_width = resilience_config$window_width %||% 7L,
      scaling_method = "none",  # Already scaled above
      capacity_types = resilience_config$capacity_types %||% c("absorptive", "restorative", "adaptive")
    )
    
    resilience_result <- do.call(calculate_resilience_metrics, resilience_args)
    
    # Classify resilience states
    state_args <- list(
      data = resilience_result,
      method = resilience_config$state_method %||% "threshold",
      ts_cols = ts_cols
    )
    resilience_result <- do.call(classify_resilience_states, state_args)
    
    pipeline_results$resilience_data <- resilience_result
    pipeline_results$resilience_summary <- get_resilience_state_summary(
      resilience_result, ts_cols = ts_cols
    )
    pipeline_results$resilience_transition_matrix <- create_resilience_transition_matrix(
      resilience_result
    )
  }
  
  # Perform STNA analysis if requested
  if ("stna" %in% analysis_types) {
    stna_args <- list(
      time_series = working_data,
      value_column = ts_cols[1],
      num_states = stna_config$num_states %||% 4,
      method = stna_config$method %||% "quantile"
    )
    
    stna_result <- do.call(stna, stna_args)
    pipeline_results$stna_result <- stna_result
    pipeline_results$stna_summary <- stna_result$state_summary
  }
  
  # Perform rolling measures analysis if requested
  if ("rolling" %in% analysis_types) {
    rolling_args <- list(
      data = working_data,
      ts_cols = ts_cols,
      window_width = rolling_config$window_width %||% 7L,
      measures = rolling_config$measures %||% c("complexity", "fluctuation", "distribution")
    )
    
    rolling_result <- do.call(rolling_measures, rolling_args)
    pipeline_results$rolling_result <- rolling_result
  }
  
  # Perform Hurst analysis if requested and available
  if ("hurst" %in% analysis_types && exists("calculate_hurst", mode = "function")) {
    tryCatch({
      hurst_args <- list(
        data = working_data,
        ts_col = ts_cols[1],
        method = hurst_config$method %||% "dfa"
      )
      
      hurst_result <- do.call(calculate_hurst, hurst_args)
      pipeline_results$hurst_result <- hurst_result
    }, error = function(e) {
      warning("Hurst analysis failed: ", e$message, call. = FALSE)
    })
  }
  
  # Create integrated data if multiple analyses were performed
  if (integrate_results && length(analysis_types) > 1) {
    pipeline_results$integrated_data <- .create_integrated_data(pipeline_results)
  }
  
  # Add metadata
  pipeline_results$config <- list(
    resilience_config = resilience_config,
    stna_config = stna_config,
    rolling_config = rolling_config,
    hurst_config = hurst_config
  )
  
  class(pipeline_results) <- c("resilience_pipeline", "list")
  return(pipeline_results)
}

#' Convert Resilience Results to TSnetworks Format
#'
#' @description
#' Ensures resilience analysis results are fully compatible with existing
#' TSnetworks functions. This function formats resilience data to work
#' seamlessly with all existing TSnetworks plotting and analysis functions.
#'
#' @param resilience_data Data frame with resilience analysis
#' @param target_format Target format: "stna", "rolling_measures", "generic"
#' @param preserve_attributes Whether to preserve resilience-specific attributes
#'
#' @return Data frame formatted for existing TSnetworks functions
#'
#' @examples
#' # Format resilience data for existing TSnetworks functions
#' data(saqrsteps)
#' resilience_data <- calculate_resilience_metrics(saqrsteps, ts_cols = "Steps")
#' state_data <- classify_resilience_states(resilience_data)
#' 
#' # Format for use with existing plotting functions
#' formatted_data <- format_for_tsnetworks(state_data, target_format = "generic")
#' 
#' # Now works with existing TSnetworks functions
#' # plot_timeseries_enhanced(formatted_data, ts_col = "Steps", state_col = "resilience_state")
#'
#' @export
format_for_tsnetworks <- function(resilience_data,
                                  target_format = c("generic", "stna", "rolling_measures"),
                                  preserve_attributes = TRUE) {
  
  target_format <- match.arg(target_format)
  
  # Validate input
  if (!is.data.frame(resilience_data)) {
    stop("'resilience_data' must be a data frame", call. = FALSE)
  }
  
  # Create formatted data
  formatted_data <- resilience_data
  
  # Ensure column naming follows TSnetworks conventions
  if (target_format == "stna") {
    # Ensure compatibility with STNA result expectations
    if (!"state" %in% names(formatted_data) && "resilience_state" %in% names(formatted_data)) {
      formatted_data$state <- formatted_data$resilience_state
    }
    
    # Add any missing STNA-expected columns
    if (!"original_value" %in% names(formatted_data)) {
      ts_cols <- .auto_detect_ts_columns(formatted_data)
      if (length(ts_cols) > 0) {
        formatted_data$original_value <- formatted_data[[ts_cols[1]]]
      }
    }
    
  } else if (target_format == "rolling_measures") {
    # Ensure compatibility with rolling measures expectations
    # Add any missing rolling measures columns as needed
    
  } else {
    # Generic format - ensure basic compatibility
    # Make sure state columns are properly named
  }
  
  # Preserve attributes if requested
  if (preserve_attributes) {
    # Copy resilience-specific attributes
    resilience_attrs <- attributes(resilience_data)
    for (attr_name in names(resilience_attrs)) {
      if (attr_name %in% c("resilience_analysis", "resilience_states")) {
        attr(formatted_data, attr_name) <- resilience_attrs[[attr_name]]
      }
    }
  }
  
  # Add compatibility metadata
  attr(formatted_data, "tsnetworks_compatible") <- TRUE
  attr(formatted_data, "formatted_for") <- target_format
  attr(formatted_data, "formatting_timestamp") <- Sys.time()
  
  return(formatted_data)
}

#' Extract Resilience Metrics from Enhanced Objects
#'
#' @description
#' Helper function to extract resilience-specific components from enhanced
#' TSnetworks objects. This allows users to access resilience results
#' separately while maintaining the integrated object structure.
#'
#' @param enhanced_object Object with resilience analysis (from add_resilience_to_tsnetworks())
#' @param component Which component to extract: "metrics", "states", "summary", "all"
#'
#' @return List of resilience metrics and states
#'
#' @examples
#' # Extract resilience components
#' data(saqrsteps)
#' stna_result <- stna(saqrsteps, "Steps", num_states = 4, method = "quantile")
#' enhanced_result <- add_resilience_to_tsnetworks(stna_result)
#' 
#' # Extract specific resilience components
#' resilience_metrics <- extract_resilience_results(enhanced_result, "metrics")
#' resilience_summary <- extract_resilience_results(enhanced_result, "summary")
#'
#' @export
extract_resilience_results <- function(enhanced_object,
                                       component = c("all", "metrics", "states", "summary")) {
  
  component <- match.arg(component)
  
  # Detect object type and extract data
  if (inherits(enhanced_object, "enhanced_stna")) {
    data <- enhanced_object$data
    resilience_info <- enhanced_object$resilience_analysis
    summary_info <- enhanced_object$resilience_summary
  } else if (inherits(enhanced_object, "enhanced_rolling_measures")) {
    data <- enhanced_object
    resilience_info <- attr(enhanced_object, "resilience_enhancement")
    summary_info <- NULL
  } else if (inherits(enhanced_object, "enhanced_tsnetworks_data")) {
    data <- enhanced_object
    resilience_info <- attr(enhanced_object, "resilience_enhancement")
    summary_info <- NULL
  } else {
    stop("Object does not appear to be an enhanced TSnetworks object", call. = FALSE)
  }
  
  # Extract requested components
  if (component == "metrics") {
    # Extract resilience metric columns
    metric_cols <- grep("^(absorptive|restorative|adaptive)_", names(data), value = TRUE)
    return(data[metric_cols])
    
  } else if (component == "states") {
    # Extract resilience state columns
    state_cols <- grep("^resilience_state", names(data), value = TRUE)
    return(data[state_cols])
    
  } else if (component == "summary") {
    # Return summary information
    return(summary_info)
    
  } else {
    # Return all resilience components
    metric_cols <- grep("^(absorptive|restorative|adaptive)_", names(data), value = TRUE)
    state_cols <- grep("^resilience_state", names(data), value = TRUE)
    
    return(list(
      metrics = data[metric_cols],
      states = data[state_cols],
      summary = summary_info,
      analysis_info = resilience_info
    ))
  }
}

#' Plot Integrated TSnetworks and Resilience Dashboard
#'
#' @description
#' Creates comprehensive visualization combining TSnetworks and resilience
#' analysis results. This function automatically detects available analyses
#' and creates an integrated dashboard.
#'
#' @param pipeline_result Result from resilience_analysis_pipeline()
#' @param ts_col Time series column to focus on
#' @param layout Dashboard layout: "grid", "vertical", "horizontal"
#' @param include_network Whether to include network plots
#' @param title Overall dashboard title
#'
#' @return Integrated dashboard plot
#'
#' @examples
#' # Create integrated dashboard
#' data(saqrsteps)
#' pipeline_result <- resilience_analysis_pipeline(
#'   data = saqrsteps,
#'   ts_cols = "Steps",
#'   analysis_types = c("resilience", "stna", "rolling")
#' )
#' plot_integrated_dashboard(pipeline_result, ts_col = "Steps")
#'
#' @export
plot_integrated_dashboard <- function(pipeline_result,
                                     ts_col,
                                     layout = c("grid", "vertical", "horizontal"),
                                     include_network = TRUE,
                                     title = "Integrated TSnetworks and Resilience Analysis") {
  
  layout <- match.arg(layout)
  
  # Validate input
  if (!inherits(pipeline_result, "resilience_pipeline")) {
    stop("Input must be a result from resilience_analysis_pipeline()", call. = FALSE)
  }
  
  # Create individual plots based on available analyses
  plots <- list()
  
  # 1. Original time series
  if ("integrated_data" %in% names(pipeline_result)) {
    data_for_plotting <- pipeline_result$integrated_data
  } else if ("resilience_data" %in% names(pipeline_result)) {
    data_for_plotting <- pipeline_result$resilience_data
  } else {
    data_for_plotting <- pipeline_result$data
  }
  
  plots$timeseries <- .create_basic_timeseries_plot(data_for_plotting, ts_col, "resilience_state")
  plots$timeseries <- plots$timeseries + labs(title = "Original Time Series")
  
  # 2. Resilience capacity timeline
  if ("resilience_data" %in% names(pipeline_result)) {
    plots$resilience <- plot_resilience_timeline(
      pipeline_result$resilience_data,
      capacity_type = "all",
      ts_col = ts_col,
      title = "Resilience Capacities"
    )
  }
  
  # 3. Rolling measures (if available)
  if ("rolling_result" %in% names(pipeline_result)) {
    rolling_data <- pipeline_result$rolling_result
    # Create rolling measures plot
    plots$rolling <- .create_rolling_measures_plot(rolling_data, ts_col)
  }
  
  # 4. Network plots (if requested and available)
  if (include_network) {
    if ("resilience_data" %in% names(pipeline_result)) {
      plots$resilience_network <- plot_resilience_network(
        pipeline_result$resilience_data,
        title = "Resilience State Network"
      )
    }
    
    if ("stna_result" %in% names(pipeline_result)) {
      # Create STNA network plot if plotting function exists
      if (exists("plot_tna_network", mode = "function")) {
        tryCatch({
          plots$stna_network <- plot_tna_network(
            pipeline_result$stna_result$data,
            state_col = "state"
          )
          plots$stna_network <- plots$stna_network + labs(title = "STNA Network")
        }, error = function(e) {
          # Skip if plotting fails
        })
      }
    }
  }
  
  # Arrange plots based on layout
  if (layout == "grid") {
    # Grid layout
    n_plots <- length(plots)
    if (n_plots <= 4) {
      combined_plot <- gridExtra::grid.arrange(grobs = plots, ncol = 2)
    } else {
      combined_plot <- gridExtra::grid.arrange(grobs = plots, ncol = 3)
    }
  } else if (layout == "vertical") {
    # Vertical layout
    combined_plot <- gridExtra::grid.arrange(grobs = plots, ncol = 1)
  } else {
    # Horizontal layout
    combined_plot <- gridExtra::grid.arrange(grobs = plots, nrow = 1)
  }
  
  # Add overall title
  if (!is.null(title)) {
    combined_plot <- gridExtra::grid.arrange(
      gridExtra::textGrob(title, gp = grid::gpar(fontsize = 16, fontface = "bold")),
      combined_plot,
      heights = c(0.05, 0.95)
    )
  }
  
  return(combined_plot)
}

# Internal helper functions

.detect_tsnetworks_result_type <- function(result) {
  if (is.list(result) && "data" %in% names(result) && "transition_matrix" %in% names(result)) {
    return("stna")
  } else if (is.data.frame(result)) {
    # Check if it looks like rolling measures result
    rolling_cols <- grep("_(complexity|fluctuation|distribution|autocorrelation|max|min|variance)$", 
                        names(result), value = TRUE)
    if (length(rolling_cols) > 0) {
      return("rolling_measures")
    } else {
      return("data_frame")
    }
  } else {
    return("unknown")
  }
}

.extract_ts_columns_from_stna <- function(stna_result) {
  # Extract time series columns from STNA result
  if ("parameters" %in% names(stna_result) && "value_column" %in% names(stna_result$parameters)) {
    return(stna_result$parameters$value_column)
  } else {
    # Try to infer from data
    return(.auto_detect_ts_columns(stna_result$data))
  }
}

.extract_ts_columns_from_rolling <- function(rolling_result) {
  # Extract original time series columns from rolling measures result
  rolling_cols <- grep("_(complexity|fluctuation|distribution|autocorrelation|max|min|variance)$", 
                      names(rolling_result), value = TRUE)
  
  if (length(rolling_cols) > 0) {
    # Extract base column names
    ts_cols <- unique(gsub("_(complexity|fluctuation|distribution|autocorrelation|max|min|variance)$", 
                          "", rolling_cols))
    return(ts_cols)
  } else {
    return(.auto_detect_ts_columns(rolling_result))
  }
}

.auto_detect_ts_columns <- function(data) {
  # Auto-detect time series columns
  numeric_cols <- names(data)[sapply(data, is.numeric)]
  
  # Exclude common non-time-series columns
  exclude_patterns <- c("^(absorptive|restorative|adaptive)_", "^resilience_state", 
                       "_score$", "_state$", "complexity$", "fluctuation$", 
                       "distribution$", "autocorrelation$")
  
  for (pattern in exclude_patterns) {
    numeric_cols <- numeric_cols[!grepl(pattern, numeric_cols)]
  }
  
  return(numeric_cols)
}

.validate_pipeline_inputs <- function(data, ts_cols, analysis_types) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame", call. = FALSE)
  }
  
  if (!all(ts_cols %in% names(data))) {
    missing_cols <- ts_cols[!ts_cols %in% names(data)]
    stop("Time series columns not found: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  
  valid_analyses <- c("resilience", "stna", "rolling", "hurst")
  if (!all(analysis_types %in% valid_analyses)) {
    invalid <- analysis_types[!analysis_types %in% valid_analyses]
    stop("Invalid analysis types: ", paste(invalid, collapse = ", "), call. = FALSE)
  }
}

.create_integrated_data <- function(pipeline_results) {
  # Start with base data
  integrated_data <- pipeline_results$data
  
  # Add resilience data if available
  if ("resilience_data" %in% names(pipeline_results)) {
    resilience_cols <- grep("^(absorptive|restorative|adaptive|resilience_state)", 
                           names(pipeline_results$resilience_data), value = TRUE)
    for (col in resilience_cols) {
      integrated_data[[col]] <- pipeline_results$resilience_data[[col]]
    }
  }
  
  # Add STNA data if available
  if ("stna_result" %in% names(pipeline_results)) {
    stna_cols <- c("state", "ranked_state", "original_value", "transformed_value")
    available_stna_cols <- intersect(stna_cols, names(pipeline_results$stna_result$data))
    for (col in available_stna_cols) {
      integrated_data[[paste0("stna_", col)]] <- pipeline_results$stna_result$data[[col]]
    }
  }
  
  # Add rolling measures if available
  if ("rolling_result" %in% names(pipeline_results)) {
    rolling_cols <- grep("_(complexity|fluctuation|distribution|autocorrelation|max|min|variance)$", 
                        names(pipeline_results$rolling_result), value = TRUE)
    for (col in rolling_cols) {
      integrated_data[[col]] <- pipeline_results$rolling_result[[col]]
    }
  }
  
  return(integrated_data)
}

.create_rolling_measures_plot <- function(rolling_data, ts_col) {
  # Create a basic rolling measures plot
  rolling_cols <- grep("_(complexity|fluctuation|distribution)$", names(rolling_data), value = TRUE)
  
  if (length(rolling_cols) == 0) {
    return(NULL)
  }
  
  # Prepare data for plotting
  rolling_data$time_index <- 1:nrow(rolling_data)
  plot_data <- data.frame()
  
  for (col in rolling_cols[1:min(3, length(rolling_cols))]) {  # Limit to first 3 measures
    measure_name <- gsub("^.+_(.+)$", "\\1", col)
    temp_data <- data.frame(
      time_index = rolling_data$time_index,
      value = rolling_data[[col]],
      measure = measure_name
    )
    plot_data <- rbind(plot_data, temp_data)
  }
  
  # Create plot
  p <- ggplot(plot_data, aes(x = time_index, y = value, color = measure)) +
    geom_line(size = 0.8) +
    facet_wrap(~ measure, scales = "free_y") +
    theme_minimal() +
    labs(title = "Rolling Measures", x = "Time Index", y = "Value") +
    theme(legend.position = "none")
  
  return(p)
}

# Utility operator for default values
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Print method for resilience_pipeline objects
#' @param x A resilience_pipeline object
#' @param ... Additional arguments (unused)
#' @export
print.resilience_pipeline <- function(x, ...) {
  cat("Resilience Analysis Pipeline Result\n")
  cat("===================================\n")
  cat("Analysis types:", paste(x$analysis_types, collapse = ", "), "\n")
  cat("Time series columns:", paste(x$ts_cols, collapse = ", "), "\n")
  cat("Analysis timestamp:", format(x$timestamp), "\n")
  
  cat("\nAvailable components:\n")
  components <- names(x)[!names(x) %in% c("data", "ts_cols", "analysis_types", "timestamp", "config")]
  for (comp in components) {
    cat("  -", comp, "\n")
  }
  
  if ("resilience_summary" %in% names(x)) {
    cat("\nResilience states found:", length(x$resilience_summary$state_names), "\n")
    cat("State names:", paste(x$resilience_summary$state_names, collapse = ", "), "\n")
  }
  
  if ("integrated_data" %in% names(x)) {
    cat("\nIntegrated data dimensions:", nrow(x$integrated_data), "x", ncol(x$integrated_data), "\n")
  }
}
