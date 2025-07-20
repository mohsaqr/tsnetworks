#' Classify Resilience States for TSnetworks
#'
#' @description
#' Creates resilience state classifications that are fully compatible with
#' existing TSnetworks state analysis and plotting functions. This function
#' analyzes resilience capacity metrics to classify time points into distinct
#' resilience states that can be used with all existing TSnetworks visualization
#' and network analysis functions.
#'
#' @param data Data frame with resilience metrics (from calculate_resilience_metrics())
#' @param method Classification method: "threshold", "kmeans", "quantile", "composite"
#' @param state_column_name Name for the main resilience state column (default: "resilience_state")
#' @param capacity_weights Weights for combining capacity scores (absorptive, restorative, adaptive)
#' @param threshold_values Custom threshold values for classification
#' @param n_states Number of states for methods that require it (default: 5)
#' @param ts_cols Time series columns to consider (auto-detected if NULL)
#'
#' @return Data frame with resilience state columns added, compatible with
#'   plot_tna_network(), plot_timeseries_enhanced(), and other TSnetworks functions
#'
#' @examples
#' # Basic resilience state classification
#' data(saqrsteps)
#' resilience_data <- calculate_resilience_metrics(saqrsteps, ts_cols = "Steps")
#' state_data <- classify_resilience_states(resilience_data)
#' 
#' # Use with existing TSnetworks plotting
#' plot_tna_network(state_data, state_col = "resilience_state")
#' plot_timeseries_enhanced(state_data, ts_col = "Steps", state_col = "resilience_state")
#'
#' # Custom classification with different weights
#' custom_states <- classify_resilience_states(
#'   resilience_data,
#'   method = "composite",
#'   capacity_weights = c(absorptive = 0.5, restorative = 0.3, adaptive = 0.2)
#' )
#'
#' @importFrom stats kmeans quantile
#' @export
classify_resilience_states <- function(data,
                                      method = c("threshold", "kmeans", "quantile", "composite"),
                                      state_column_name = "resilience_state",
                                      capacity_weights = c(absorptive = 0.4, restorative = 0.3, adaptive = 0.3),
                                      threshold_values = NULL,
                                      n_states = 5,
                                      ts_cols = NULL) {
  
  method <- match.arg(method)
  
  # Validate inputs
  .validate_state_classification_inputs(data, method, state_column_name, capacity_weights, n_states)
  
  # Auto-detect time series columns if needed
  if (is.null(ts_cols)) {
    # Look for resilience metric columns to infer original ts_cols
    resilience_cols <- grep("^(absorptive|restorative|adaptive)_", names(data), value = TRUE)
    if (length(resilience_cols) > 0) {
      # Extract original column names from resilience metric column names
      ts_cols <- unique(gsub("^(absorptive|restorative|adaptive)_[^_]+_", "", resilience_cols))
    } else {
      # Fallback to numeric columns
      ts_cols <- names(data)[sapply(data, is.numeric)]
    }
  }
  
  # Calculate composite capacity scores
  capacity_scores <- .calculate_composite_capacity_scores(data, ts_cols, capacity_weights)
  
  # Classify states based on method
  if (method == "threshold") {
    state_classifications <- .classify_states_threshold(capacity_scores, threshold_values)
  } else if (method == "kmeans") {
    state_classifications <- .classify_states_kmeans(capacity_scores, n_states)
  } else if (method == "quantile") {
    state_classifications <- .classify_states_quantile(capacity_scores, n_states)
  } else if (method == "composite") {
    state_classifications <- .classify_states_composite(capacity_scores, threshold_values)
  }
  
  # Add state classifications to data
  result_data <- data
  result_data[[state_column_name]] <- state_classifications$main_state
  
  # Add individual capacity states
  result_data[[paste0(state_column_name, "_absorptive")]] <- state_classifications$absorptive_state
  result_data[[paste0(state_column_name, "_restorative")]] <- state_classifications$restorative_state
  result_data[[paste0(state_column_name, "_adaptive")]] <- state_classifications$adaptive_state
  
  # Add capacity scores for reference
  result_data[[paste0(state_column_name, "_absorptive_score")]] <- capacity_scores$absorptive
  result_data[[paste0(state_column_name, "_restorative_score")]] <- capacity_scores$restorative
  result_data[[paste0(state_column_name, "_adaptive_score")]] <- capacity_scores$adaptive
  result_data[[paste0(state_column_name, "_composite_score")]] <- capacity_scores$composite
  
  # Add metadata attributes
  attr(result_data, "resilience_states") <- list(
    method = method,
    state_column_name = state_column_name,
    capacity_weights = capacity_weights,
    n_states = n_states,
    state_definitions = .get_resilience_state_definitions(),
    classification_timestamp = Sys.time()
  )
  
  class(result_data) <- c("resilience_states", class(result_data))
  return(result_data)
}

#' Create Resilience State Transition Matrix
#'
#' @description
#' Creates transition matrices compatible with existing TSnetworks network analysis.
#' The resulting transition matrix can be used with existing network visualization
#' and analysis functions.
#'
#' @param data Data frame with resilience states (from classify_resilience_states())
#' @param state_col Resilience state column name (default: "resilience_state")
#' @param normalize Whether to normalize to probabilities (default: TRUE)
#' @param id_col Grouping column for user-level analysis (optional)
#'
#' @return Transition matrix compatible with existing TSnetworks network functions
#'
#' @examples
#' # Create resilience state transition matrix
#' data(saqrsteps)
#' resilience_data <- calculate_resilience_metrics(saqrsteps, ts_cols = "Steps")
#' state_data <- classify_resilience_states(resilience_data)
#' transition_matrix <- create_resilience_transition_matrix(state_data)
#' print(transition_matrix)
#'
#' # Use with existing network analysis
#' # (This would work with any existing TSnetworks network analysis functions)
#'
#' @export
create_resilience_transition_matrix <- function(data,
                                               state_col = "resilience_state",
                                               normalize = TRUE,
                                               id_col = NULL) {
  
  # Validate inputs
  if (!state_col %in% names(data)) {
    stop(paste("State column", state_col, "not found in data"), call. = FALSE)
  }
  
  states <- data[[state_col]]
  
  if (all(is.na(states))) {
    stop("All state values are NA", call. = FALSE)
  }
  
  # Remove NA values
  valid_indices <- !is.na(states)
  states <- states[valid_indices]
  
  if (length(states) < 2) {
    stop("Need at least 2 valid state observations for transition analysis", call. = FALSE)
  }
  
  # Calculate transitions
  if (!is.null(id_col) && id_col %in% names(data)) {
    # Grouped calculation (same logic as existing TSnetworks functions)
    groups <- split(seq_along(data[[state_col]]), data[[id_col]], drop = TRUE)
    
    # Get unique states
    unique_states <- sort(unique(states))
    n_states <- length(unique_states)
    
    # Initialize transition matrix
    transition_matrix <- matrix(0, nrow = n_states, ncol = n_states,
                               dimnames = list(unique_states, unique_states))
    
    # Calculate transitions for each group
    for (group_indices in groups) {
      group_states <- data[[state_col]][group_indices]
      group_states <- group_states[!is.na(group_states)]
      
      if (length(group_states) >= 2) {
        # Calculate transitions for this group
        for (i in 1:(length(group_states) - 1)) {
          from_state <- as.character(group_states[i])
          to_state <- as.character(group_states[i + 1])
          
          if (from_state %in% unique_states && to_state %in% unique_states) {
            transition_matrix[from_state, to_state] <- transition_matrix[from_state, to_state] + 1
          }
        }
      }
    }
    
  } else {
    # Global calculation
    unique_states <- sort(unique(states))
    n_states <- length(unique_states)
    
    # Initialize transition matrix
    transition_matrix <- matrix(0, nrow = n_states, ncol = n_states,
                               dimnames = list(unique_states, unique_states))
    
    # Calculate transitions
    for (i in 1:(length(states) - 1)) {
      from_state <- as.character(states[i])
      to_state <- as.character(states[i + 1])
      
      transition_matrix[from_state, to_state] <- transition_matrix[from_state, to_state] + 1
    }
  }
  
  # Normalize if requested
  if (normalize) {
    # Convert to probabilities (same logic as existing TSnetworks functions)
    row_sums <- rowSums(transition_matrix)
    for (i in 1:nrow(transition_matrix)) {
      if (row_sums[i] > 0) {
        transition_matrix[i, ] <- transition_matrix[i, ] / row_sums[i]
      }
    }
  }
  
  # Add attributes for compatibility
  attr(transition_matrix, "type") <- "resilience_transitions"
  attr(transition_matrix, "normalized") <- normalize
  attr(transition_matrix, "state_column") <- state_col
  
  return(transition_matrix)
}

#' Get Resilience State Summary Statistics
#'
#' @description
#' Provides comprehensive summary statistics for resilience states,
#' compatible with existing TSnetworks analysis patterns.
#'
#' @param data Data frame with resilience states
#' @param state_col Resilience state column name
#' @param ts_cols Time series columns to analyze
#'
#' @return List with state summary statistics
#'
#' @examples
#' data(saqrsteps)
#' resilience_data <- calculate_resilience_metrics(saqrsteps, ts_cols = "Steps")
#' state_data <- classify_resilience_states(resilience_data)
#' summary_stats <- get_resilience_state_summary(state_data, ts_cols = "Steps")
#' print(summary_stats)
#'
#' @export
get_resilience_state_summary <- function(data, state_col = "resilience_state", ts_cols = NULL) {
  
  if (!state_col %in% names(data)) {
    stop(paste("State column", state_col, "not found in data"), call. = FALSE)
  }
  
  # Auto-detect ts_cols if needed
  if (is.null(ts_cols)) {
    ts_cols <- names(data)[sapply(data, is.numeric)]
    ts_cols <- setdiff(ts_cols, grep("^(absorptive|restorative|adaptive)_", ts_cols, value = TRUE))
  }
  
  states <- data[[state_col]]
  unique_states <- sort(unique(states[!is.na(states)]))
  
  # Calculate state statistics
  state_stats <- list()
  
  for (state in unique_states) {
    state_indices <- which(states == state)
    state_data <- data[state_indices, ]
    
    state_info <- list(
      count = length(state_indices),
      proportion = length(state_indices) / nrow(data),
      duration_stats = .calculate_state_durations(states, state),
      capacity_scores = list()
    )
    
    # Calculate capacity score statistics for this state
    capacity_cols <- grep("^resilience_state_(absorptive|restorative|adaptive)_score$", names(data), value = TRUE)
    for (cap_col in capacity_cols) {
      if (cap_col %in% names(data)) {
        cap_values <- state_data[[cap_col]]
        state_info$capacity_scores[[cap_col]] <- list(
          mean = mean(cap_values, na.rm = TRUE),
          median = median(cap_values, na.rm = TRUE),
          sd = sd(cap_values, na.rm = TRUE),
          min = min(cap_values, na.rm = TRUE),
          max = max(cap_values, na.rm = TRUE)
        )
      }
    }
    
    # Calculate time series statistics for this state
    state_info$ts_stats <- list()
    for (ts_col in ts_cols) {
      if (ts_col %in% names(data)) {
        ts_values <- state_data[[ts_col]]
        state_info$ts_stats[[ts_col]] <- list(
          mean = mean(ts_values, na.rm = TRUE),
          median = median(ts_values, na.rm = TRUE),
          sd = sd(ts_values, na.rm = TRUE),
          min = min(ts_values, na.rm = TRUE),
          max = max(ts_values, na.rm = TRUE)
        )
      }
    }
    
    state_stats[[state]] <- state_info
  }
  
  # Calculate transition statistics
  transition_matrix <- create_resilience_transition_matrix(data, state_col, normalize = TRUE)
  
  # Overall summary
  summary_info <- list(
    n_states = length(unique_states),
    state_names = unique_states,
    total_observations = nrow(data),
    valid_observations = sum(!is.na(states)),
    state_statistics = state_stats,
    transition_matrix = transition_matrix,
    state_definitions = .get_resilience_state_definitions()
  )
  
  class(summary_info) <- c("resilience_state_summary", "list")
  return(summary_info)
}

# Internal helper functions

.validate_state_classification_inputs <- function(data, method, state_column_name, capacity_weights, n_states) {
  
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
  
  if (!is.character(state_column_name) || length(state_column_name) != 1) {
    stop("'state_column_name' must be a single character string", call. = FALSE)
  }
  
  if (!is.numeric(capacity_weights) || length(capacity_weights) != 3) {
    stop("'capacity_weights' must be a numeric vector of length 3", call. = FALSE)
  }
  
  if (abs(sum(capacity_weights) - 1) > 1e-6) {
    warning("capacity_weights do not sum to 1. Normalizing weights.", call. = FALSE)
  }
  
  if (!is.numeric(n_states) || length(n_states) != 1 || n_states < 2) {
    stop("'n_states' must be a single integer >= 2", call. = FALSE)
  }
}

.calculate_composite_capacity_scores <- function(data, ts_cols, capacity_weights) {
  
  # Normalize weights
  capacity_weights <- capacity_weights / sum(capacity_weights)
  
  n_obs <- nrow(data)
  
  # Initialize capacity scores
  absorptive_scores <- rep(NA_real_, n_obs)
  restorative_scores <- rep(NA_real_, n_obs)
  adaptive_scores <- rep(NA_real_, n_obs)
  
  # Calculate capacity scores for each time series column
  for (ts_col in ts_cols) {
    
    # Absorptive capacity metrics
    vsi_col <- paste0("absorptive_vsi_", ts_col)
    arch_col <- paste0("absorptive_arch_", ts_col)
    cv_col <- paste0("absorptive_cv_", ts_col)
    
    if (all(c(vsi_col, arch_col, cv_col) %in% names(data))) {
      # Normalize metrics to [0,1] scale
      vsi_norm <- .normalize_to_01(data[[vsi_col]])
      arch_norm <- 1 - .normalize_to_01(data[[arch_col]])  # Lower ARCH is better
      cv_norm <- 1 - .normalize_to_01(data[[cv_col]])      # Lower CV is better
      
      # Combine absorptive metrics (equal weights)
      absorptive_col <- (vsi_norm + arch_norm + cv_norm) / 3
      
      # Average across time series columns
      if (all(is.na(absorptive_scores))) {
        absorptive_scores <- absorptive_col
      } else {
        absorptive_scores <- (absorptive_scores + absorptive_col) / 2
      }
    }
    
    # Restorative capacity metrics
    recovery_time_col <- paste0("restorative_recovery_time_", ts_col)
    eng_resilience_col <- paste0("restorative_eng_resilience_", ts_col)
    recovery_slope_col <- paste0("restorative_recovery_slope_", ts_col)
    
    if (all(c(recovery_time_col, eng_resilience_col, recovery_slope_col) %in% names(data))) {
      # Normalize metrics to [0,1] scale
      recovery_time_norm <- 1 - .normalize_to_01(data[[recovery_time_col]])  # Lower recovery time is better
      eng_resilience_norm <- .normalize_to_01(data[[eng_resilience_col]])
      recovery_slope_norm <- .normalize_to_01(data[[recovery_slope_col]])
      
      # Combine restorative metrics (equal weights)
      restorative_col <- (recovery_time_norm + eng_resilience_norm + recovery_slope_norm) / 3
      
      # Average across time series columns
      if (all(is.na(restorative_scores))) {
        restorative_scores <- restorative_col
      } else {
        restorative_scores <- (restorative_scores + restorative_col) / 2
      }
    }
    
    # Adaptive capacity metrics
    entropy_col <- paste0("adaptive_sample_entropy_", ts_col)
    dfa_col <- paste0("adaptive_dfa_", ts_col)
    ratio_col <- paste0("adaptive_capacity_ratio_", ts_col)
    
    if (all(c(entropy_col, dfa_col, ratio_col) %in% names(data))) {
      # Normalize metrics to [0,1] scale
      entropy_norm <- .normalize_to_01(data[[entropy_col]])
      dfa_norm <- .normalize_to_01(abs(data[[dfa_col]] - 0.5))  # Closer to 0.5 is better
      dfa_norm <- 1 - dfa_norm  # Invert so higher is better
      ratio_norm <- .normalize_to_01(data[[ratio_col]])
      
      # Combine adaptive metrics (equal weights)
      adaptive_col <- (entropy_norm + dfa_norm + ratio_norm) / 3
      
      # Average across time series columns
      if (all(is.na(adaptive_scores))) {
        adaptive_scores <- adaptive_col
      } else {
        adaptive_scores <- (adaptive_scores + adaptive_col) / 2
      }
    }
  }
  
  # Calculate composite score
  composite_scores <- capacity_weights[1] * absorptive_scores + 
                     capacity_weights[2] * restorative_scores + 
                     capacity_weights[3] * adaptive_scores
  
  return(list(
    absorptive = absorptive_scores,
    restorative = restorative_scores,
    adaptive = adaptive_scores,
    composite = composite_scores
  ))
}

.normalize_to_01 <- function(x) {
  if (all(is.na(x))) return(x)
  
  x_clean <- x[!is.na(x)]
  if (length(x_clean) == 0) return(x)
  
  min_val <- min(x_clean)
  max_val <- max(x_clean)
  
  if (min_val == max_val) {
    return(rep(0.5, length(x)))  # If all values are the same, return middle value
  }
  
  normalized <- (x - min_val) / (max_val - min_val)
  return(normalized)
}

.classify_states_threshold <- function(capacity_scores, threshold_values) {
  
  # Default thresholds if not provided
  if (is.null(threshold_values)) {
    threshold_values <- list(
      high = 0.8,
      medium_high = 0.6,
      medium = 0.4,
      low = 0.2
    )
  }
  
  composite <- capacity_scores$composite
  absorptive <- capacity_scores$absorptive
  restorative <- capacity_scores$restorative
  adaptive <- capacity_scores$adaptive
  
  # Classify main state based on composite score
  main_state <- character(length(composite))
  main_state[composite >= threshold_values$high] <- "resilient_stable"
  main_state[composite >= threshold_values$medium_high & composite < threshold_values$high] <- "resilient_robust"
  main_state[composite >= threshold_values$medium & composite < threshold_values$medium_high] <- "resilient_recovering"
  main_state[composite >= threshold_values$low & composite < threshold_values$medium] <- "resilient_vulnerable"
  main_state[composite < threshold_values$low] <- "resilient_critical"
  main_state[is.na(composite)] <- NA_character_
  
  # Classify individual capacity states
  absorptive_state <- .classify_individual_capacity(absorptive, "absorptive")
  restorative_state <- .classify_individual_capacity(restorative, "restorative")
  adaptive_state <- .classify_individual_capacity(adaptive, "adaptive")
  
  return(list(
    main_state = main_state,
    absorptive_state = absorptive_state,
    restorative_state = restorative_state,
    adaptive_state = adaptive_state
  ))
}

.classify_states_kmeans <- function(capacity_scores, n_states) {
  
  # Prepare data for clustering
  cluster_data <- data.frame(
    absorptive = capacity_scores$absorptive,
    restorative = capacity_scores$restorative,
    adaptive = capacity_scores$adaptive
  )
  
  # Remove rows with any NA values
  complete_rows <- complete.cases(cluster_data)
  
  if (sum(complete_rows) < n_states) {
    warning("Not enough complete observations for k-means clustering. Using threshold method.", call. = FALSE)
    return(.classify_states_threshold(capacity_scores, NULL))
  }
  
  # Perform k-means clustering
  kmeans_result <- kmeans(cluster_data[complete_rows, ], centers = n_states, nstart = 20)
  
  # Map clusters to resilience state names
  cluster_centers <- kmeans_result$centers
  cluster_composite <- rowMeans(cluster_centers)
  cluster_order <- order(cluster_composite, decreasing = TRUE)
  
  state_names <- c("resilient_stable", "resilient_robust", "resilient_recovering", 
                   "resilient_vulnerable", "resilient_critical")[1:n_states]
  
  # Create state assignments
  main_state <- rep(NA_character_, length(capacity_scores$composite))
  main_state[complete_rows] <- state_names[kmeans_result$cluster]
  
  # Classify individual capacities
  absorptive_state <- .classify_individual_capacity(capacity_scores$absorptive, "absorptive")
  restorative_state <- .classify_individual_capacity(capacity_scores$restorative, "restorative")
  adaptive_state <- .classify_individual_capacity(capacity_scores$adaptive, "adaptive")
  
  return(list(
    main_state = main_state,
    absorptive_state = absorptive_state,
    restorative_state = restorative_state,
    adaptive_state = adaptive_state
  ))
}

.classify_states_quantile <- function(capacity_scores, n_states) {
  
  composite <- capacity_scores$composite
  
  # Calculate quantile thresholds
  quantile_probs <- seq(0, 1, length.out = n_states + 1)
  thresholds <- quantile(composite, probs = quantile_probs, na.rm = TRUE)
  
  # Classify states
  main_state <- character(length(composite))
  state_names <- c("resilient_critical", "resilient_vulnerable", "resilient_recovering", 
                   "resilient_robust", "resilient_stable")[1:n_states]
  
  for (i in 1:n_states) {
    if (i == 1) {
      mask <- composite >= thresholds[i] & composite <= thresholds[i + 1]
    } else {
      mask <- composite > thresholds[i] & composite <= thresholds[i + 1]
    }
    main_state[mask & !is.na(mask)] <- state_names[i]
  }
  
  main_state[is.na(composite)] <- NA_character_
  
  # Classify individual capacities
  absorptive_state <- .classify_individual_capacity(capacity_scores$absorptive, "absorptive")
  restorative_state <- .classify_individual_capacity(capacity_scores$restorative, "restorative")
  adaptive_state <- .classify_individual_capacity(capacity_scores$adaptive, "adaptive")
  
  return(list(
    main_state = main_state,
    absorptive_state = absorptive_state,
    restorative_state = restorative_state,
    adaptive_state = adaptive_state
  ))
}

.classify_states_composite <- function(capacity_scores, threshold_values) {
  
  # Use threshold method but with more sophisticated logic
  threshold_result <- .classify_states_threshold(capacity_scores, threshold_values)
  
  # Refine classifications based on individual capacity patterns
  main_state <- threshold_result$main_state
  absorptive <- capacity_scores$absorptive
  restorative <- capacity_scores$restorative
  adaptive <- capacity_scores$adaptive
  
  # Special cases based on capacity combinations
  for (i in seq_along(main_state)) {
    if (is.na(main_state[i])) next
    
    abs_val <- absorptive[i]
    res_val <- restorative[i]
    ada_val <- adaptive[i]
    
    if (!is.na(abs_val) && !is.na(res_val) && !is.na(ada_val)) {
      # High absorptive but low restorative = vulnerable
      if (abs_val > 0.7 && res_val < 0.3) {
        main_state[i] <- "resilient_vulnerable"
      }
      # High restorative but low absorptive = recovering
      else if (res_val > 0.7 && abs_val < 0.3) {
        main_state[i] <- "resilient_recovering"
      }
      # Balanced high performance = stable
      else if (abs_val > 0.6 && res_val > 0.6 && ada_val > 0.6) {
        main_state[i] <- "resilient_stable"
      }
      # All low = critical
      else if (abs_val < 0.3 && res_val < 0.3 && ada_val < 0.3) {
        main_state[i] <- "resilient_critical"
      }
    }
  }
  
  return(list(
    main_state = main_state,
    absorptive_state = threshold_result$absorptive_state,
    restorative_state = threshold_result$restorative_state,
    adaptive_state = threshold_result$adaptive_state
  ))
}

.classify_individual_capacity <- function(capacity_values, capacity_type) {
  
  state <- character(length(capacity_values))
  
  # Use quantile-based classification
  q25 <- quantile(capacity_values, 0.25, na.rm = TRUE)
  q50 <- quantile(capacity_values, 0.50, na.rm = TRUE)
  q75 <- quantile(capacity_values, 0.75, na.rm = TRUE)
  
  state[capacity_values >= q75] <- paste0(capacity_type, "_high")
  state[capacity_values >= q50 & capacity_values < q75] <- paste0(capacity_type, "_medium")
  state[capacity_values >= q25 & capacity_values < q50] <- paste0(capacity_type, "_low")
  state[capacity_values < q25] <- paste0(capacity_type, "_very_low")
  state[is.na(capacity_values)] <- NA_character_
  
  return(state)
}

.calculate_state_durations <- function(states, target_state) {
  
  # Calculate run lengths for the target state
  state_runs <- rle(states == target_state & !is.na(states))
  target_runs <- state_runs$lengths[state_runs$values]
  
  if (length(target_runs) == 0) {
    return(list(
      mean_duration = 0,
      median_duration = 0,
      max_duration = 0,
      min_duration = 0,
      n_episodes = 0
    ))
  }
  
  return(list(
    mean_duration = mean(target_runs),
    median_duration = median(target_runs),
    max_duration = max(target_runs),
    min_duration = min(target_runs),
    n_episodes = length(target_runs)
  ))
}

.get_resilience_state_definitions <- function() {
  
  list(
    "resilient_stable" = list(
      description = "High performance across all resilience capacities",
      characteristics = "High absorptive and restorative capacity, good adaptive capacity",
      color = "#2E8B57",  # Sea Green
      priority = 1
    ),
    
    "resilient_robust" = list(
      description = "Balanced performance across resilience capacities",
      characteristics = "Moderate to high performance in all three capacity dimensions",
      color = "#4169E1",  # Royal Blue
      priority = 2
    ),
    
    "resilient_recovering" = list(
      description = "Strong recovery capabilities, rebuilding resilience",
      characteristics = "Good restorative capacity, improving adaptive capacity",
      color = "#FF8C00",  # Dark Orange
      priority = 3
    ),
    
    "resilient_vulnerable" = list(
      description = "Low absorptive capacity, at-risk state",
      characteristics = "Weak ability to absorb disturbances, may have some recovery capacity",
      color = "#DC143C",  # Crimson
      priority = 4
    ),
    
    "resilient_critical" = list(
      description = "Very low resilience across all capacities",
      characteristics = "Poor performance in absorptive, restorative, and adaptive capacities",
      color = "#8B0000",  # Dark Red
      priority = 5
    )
  )
}

#' Print method for resilience_states objects
#' @param x A resilience_states object
#' @param ... Additional arguments (unused)
#' @export
print.resilience_states <- function(x, ...) {
  cat("Resilience States Data\n")
  cat("======================\n")
  
  # Get resilience states info
  states_info <- attr(x, "resilience_states")
  
  if (!is.null(states_info)) {
    cat("Classification method:", states_info$method, "\n")
    cat("State column name:", states_info$state_column_name, "\n")
    cat("Capacity weights:", paste(names(states_info$capacity_weights), "=", 
        round(states_info$capacity_weights, 3), collapse = ", "), "\n")
    cat("Classification timestamp:", format(states_info$classification_timestamp), "\n")
  }
  
  cat("Dimensions:", nrow(x), "x", ncol(x), "\n")
  
  # Show state distribution
  main_state_col <- if (!is.null(states_info)) states_info$state_column_name else "resilience_state"
  if (main_state_col %in% names(x)) {
    state_counts <- table(x[[main_state_col]], useNA = "ifany")
    cat("\nState distribution:\n")
    print(state_counts)
  }
  
  # Identify resilience columns
  resilience_cols <- grep("^resilience_state_", names(x), value = TRUE)
  if (length(resilience_cols) > 0) {
    cat("\nResilience state columns:", length(resilience_cols), "\n")
    if (length(resilience_cols) <= 8) {
      cat("  ", paste(resilience_cols, collapse = ", "), "\n")
    } else {
      cat("  ", paste(resilience_cols[1:8], collapse = ", "), "...\n")
    }
  }
  
  cat("\nFirst few rows:\n")
  print(head(x))
}

#' Print method for resilience_state_summary objects
#' @param x A resilience_state_summary object
#' @param ... Additional arguments (unused)
#' @export
print.resilience_state_summary <- function(x, ...) {
  cat("Resilience State Summary\n")
  cat("========================\n")
  cat("Number of states:", x$n_states, "\n")
  cat("State names:", paste(x$state_names, collapse = ", "), "\n")
  cat("Total observations:", x$total_observations, "\n")
  cat("Valid observations:", x$valid_observations, "\n")
  
  cat("\nState Statistics:\n")
  for (state_name in x$state_names) {
    state_info <- x$state_statistics[[state_name]]
    cat("\n", state_name, ":\n")
    cat("  Count:", state_info$count, "(", round(state_info$proportion * 100, 1), "%)\n")
    cat("  Episodes:", state_info$duration_stats$n_episodes, "\n")
    cat("  Mean duration:", round(state_info$duration_stats$mean_duration, 2), "\n")
    cat("  Max duration:", state_info$duration_stats$max_duration, "\n")
  }
  
  cat("\nTransition Matrix:\n")
  print(round(x$transition_matrix, 3))
}
