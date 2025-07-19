#' @title Markov Transition Network Generator
#' @description Main interface for generating Markov transition networks
#'
#' @importFrom stats na.omit quantile
#' @importFrom dplyr %>% group_by mutate ungroup
#' @name stna-package
NULL

#' Sequence/Time Series Network Analysis (STNA)
#'
#' This function provides a comprehensive framework for analyzing time series data
#' by converting it into discrete states and constructing Markov transition networks.
#' It supports various discretization methods, global or user-level analysis,
#' and optional state ranking and data transformation.
#'
#' @param time_series A data frame containing the time series data. It must include
#'   a column for the time series values and optionally a column for user/group IDs.
#' @param num_states Integer. The desired number of discrete states to discretize
#'   the time series into.
#' @param method Character string. The discretization method to use. Options include:
#'   "kmeans", "quantile", "equal_width", "entropy", "mixture", "hierarchical",
#'   "change_points", "dtw", and "proxy_windowed". Defaults to "kmeans".
#' @param analysis_level Character string. Specifies whether to perform analysis
#'   at a "global" level (treating the entire `time_series` as one unit) or at a
#'   "user" level (performing analysis independently for each unique user/group
#'   defined by `user_column` and then aggregating transitions). Defaults to "global".
#' @param window_size Integer. Relevant for window-based discretization methods
#'   like "proxy_windowed". Specifies the size of the sliding window. Defaults to 1.
#' @param distance_metric Character string. Relevant for distance-based discretization
#'   methods like "proxy_windowed". Specifies the distance metric to use (e.g., "euclidean").
#'   Defaults to "euclidean".
#' @param user_column Character string. The name of the column in `time_series` that
#'   contains user or group identifiers. Required if `analysis_level` is "user".
#'   Defaults to `NULL`.
#' @param value_column Character string. The name of the column in `time_series` that
#'   contains the numeric time series values to be analyzed.
#' @param rank_states Logical. If `TRUE` (default), the discrete states will be
#'   ranked based on their mean original time series values.
#' @param transform Character string. Specifies a transformation to apply to the
#'   time series data before discretization. Options: "none" (default), "log",
#'   "sqrt", "scale" (z-score standardization), "center" (mean-centering).
#'
#' @return A list containing comprehensive results of the STNA, including:
#'   \item{data}{The original input data frame augmented with `state`, `original_value`,
#'     `transformed_value`, and `ranked_state` (if `rank_states` is TRUE) columns.}
#'   \item{states}{A vector of the discrete state assignments for each data point.}
#'   \item{ranked_states}{A vector of ranked state assignments if `rank_states` is TRUE.}
#'   \item{state_mapping}{A mapping from original state IDs to ranked state IDs if `rank_states` is TRUE.}
#'   \item{transition_matrix}{The Markov transition probability matrix between states.}
#'   \item{global_statistics}{Detailed statistics for each state across the entire time series.}
#'   \item{state_summary}{A summary data frame of state statistics.}
#'   \item{user_results}{A list of results for each user/group if `analysis_level` is "user".}
#'   \item{method}{The discretization method used.}
#'   \item{breaks}{Information about the breaks used for discretization (method-dependent).}
#'   \item{analysis_level}{The analysis level performed ("global" or "user").}
#'   \item{transformation}{The transformation method applied.}
#'   \item{transformation_parameters}{Parameters used for the transformation.}
#'   \item{parameters}{A list of all parameters used in the function call.}
#' @export
stna <- function(time_series,
                num_states,
                method = "kmeans",
                analysis_level = "global",
                window_size = 1,
                distance_metric = "euclidean",
                user_column = NULL,
                value_column,
                rank_states = TRUE,
                transform = "none") {
    
    # Validate and prepare input data
    validated_input <- .validate_input(time_series, value_column, user_column)
    
    # Store original values
    original_values <- validated_input$series
    
    # Transform series if requested
    transformed <- transform_series(validated_input$series, transform)
    validated_input$series <- transformed$transformed
    
    # Prepare discretization arguments
    discretization_fun <- get_method(method)
    discretization_args <- list(
        num_states = num_states,
        state_info = list(
            window_size = window_size,
            distance_metric = distance_metric
        ),
        suppress_warnings = FALSE
    )

    # If window_size is set and user_column is not NULL, process per user
    if (!is.null(window_size) && window_size > 1 && !is.null(user_column)) {
        user_ids <- validated_input$users
        all_states <- rep(NA, length(validated_input$series))
        unique_users <- unique(user_ids)
        for (user in unique_users) {
            idx <- which(user_ids == user)
            user_series <- validated_input$series[idx]
            user_args <- discretization_args
            user_args$series <- user_series
            user_result <- do.call(discretization_fun, user_args)
            all_states[idx] <- user_result$states
        }
        discretization_result <- list(states = all_states, breaks = NA)
    } else {
        discretization_args$series <- validated_input$series
        discretization_result <- do.call(discretization_fun, discretization_args)
    }
    
    # Calculate global statistics
    global_stats <- calculate_state_statistics(
        discretization_result$states,
        validated_input$series,
        validated_input$data,
        value_column
    )
    
    # Initialize storage for user-level results
    user_results <- list()

    # Calculate transition matrices and statistics
    if (analysis_level == "global") {
        transitions <- calculate_transitions(discretization_result$states, num_states)
        transition_matrix <- normalize_transitions(transitions)
        
        # Always build user_results, even in global mode
        if (!is.null(validated_input$users)) {
            user_splits <- split(seq_along(validated_input$series), validated_input$users)
            for (user in names(user_splits)) {
                user_indices <- user_splits[[user]]
                user_states <- discretization_result$states[user_indices]
                user_values <- validated_input$series[user_indices]
                user_data <- validated_input$data[user_indices, ]
                user_stats <- calculate_state_statistics(user_states, user_values, user_data, value_column)
                # Augment user_data with state, value, and transformed_value
                user_data$state <- user_states
                user_data$value <- original_values[user_indices]
                user_data$transformed_value <- transformed$transformed[user_indices]
                user_results[[user]] <- list(
                    states = user_states,
                    statistics = user_stats,
                    data = user_data
                )
            }
        }
    } else if (analysis_level == "user") {
        if (is.null(validated_input$users)) {
            stop("User-level analysis requested but no user_column provided")
        }
        
        # Split data by user
        user_splits <- split(seq_along(validated_input$series), validated_input$users)
        
        # Initialize aggregate transition counts matrix
        aggregate_transitions <- matrix(0, num_states, num_states)
        
        # Calculate per-user results
        for (user in names(user_splits)) {
            user_indices <- user_splits[[user]]
            user_states <- discretization_result$states[user_indices]
            user_values <- validated_input$series[user_indices]
            user_data <- validated_input$data[user_indices, ]
            
            # Ensure at least 2 states for transitions
            if (length(user_states) >= 2) {
                # Calculate user-specific transitions and statistics
                user_transitions <- calculate_transitions(user_states, num_states)
                user_trans_matrix <- normalize_transitions(user_transitions)
                user_stats <- calculate_state_statistics(user_states, user_values, user_data, value_column)
                
                # Add to aggregate transitions
                aggregate_transitions <- aggregate_transitions + user_transitions
                
                # Prepare user-specific output data
                user_output_data <- user_data
                user_output_data$state <- user_states
                user_output_data$original_value <- original_values[user_indices]
                user_output_data$transformed_value <- transformed$transformed[user_indices]
                
                user_results[[user]] <- list(
                    states = user_states,
                    transitions = user_transitions,
                    transition_matrix = user_trans_matrix,
                    statistics = user_stats,
                    data = user_output_data
                )
            }
        }
        
        # Calculate aggregate transition matrix
        if (length(user_results) > 0) {
            transition_matrix <- normalize_transitions(aggregate_transitions)
        } else {
            stop("No valid user transitions found")
        }
    } else {
        stop("Invalid analysis_level. Must be 'global' or 'user'")
    }
    
    # Optional state ranking
    if (rank_states) {
        ranking <- rank_states_by_mean(discretization_result$states, validated_input$series)
        state_mapping <- ranking$mapping
    } else {
        state_mapping <- NULL
    }
    
    # Prepare output data frame
    output_data <- validated_input$data
    output_data$state <- discretization_result$states
    output_data$original_value <- original_values
    output_data$transformed_value <- transformed$transformed
    if (rank_states) {
        output_data$ranked_state <- ranking$ranked_states
    }
    
    # Return comprehensive results
    list(
        # Data and states
        data = output_data,
        states = discretization_result$states,
        ranked_states = if (rank_states) ranking$ranked_states else NULL,
        state_mapping = state_mapping,
        
        # Transition information
        transition_matrix = transition_matrix,
        
        # Statistics and parameters
        global_statistics = global_stats$detailed,
        state_summary = global_stats$summary,
        user_results = user_results,
        
        # Method information
        method = method,
        breaks = discretization_result$breaks,
        analysis_level = analysis_level,
        
        # Transformation information
        transformation = transformed$method,
        transformation_parameters = transformed$parameters,
        
        parameters = list(
            num_states = num_states,
            window_size = window_size,
            distance_metric = distance_metric,
            rank_states = rank_states,
            transform = transform
        )
    )
}
