#' @importFrom stats median sd

#' @title Ranking Utilities
#' @description Functions for ranking and ordering states

#' Rank Results by Mean
#' @param states Original state assignments
#' @param series Original time series values
#' @return List with new states and mapping
#' @noRd
.rank_results_by_mean <- function(states, series) {
    if (length(states) != length(series))
        stop("Length mismatch between states and series")
    
    # Calculate mean value per state
    state_means <- tapply(series, states, mean, na.rm = TRUE)
    
    # Create ranking map
    rank_map <- order(state_means)
    names(rank_map) <- sort(unique(states))
    
    # Apply ranking
    new_states <- rank_map[as.character(states)]
    
    list(
        states = new_states,
        rank_map = rank_map,
        state_means = state_means
    )
}

#' Apply Ranking Map
#' @param states State assignments to transform
#' @param rank_map Mapping from old to new states
#' @return Transformed state assignments
#' @noRd
.apply_ranking_map <- function(states, rank_map) {
    if (!all(states %in% names(rank_map)))
        stop("Some states not found in ranking map")
    
    rank_map[as.character(states)]
}

#' Rank states by their mean values
#' @param states Vector of state assignments
#' @param values Original time series values
#' @return List containing ranked states and mapping
#' @noRd
rank_states_by_mean <- function(states, values) {
    # Calculate mean value for each state
    state_means <- tapply(values, states, mean)
    
    # Create rank mapping (from original state to ranked state)
    mapping <- rank(state_means, ties.method = "first")
    names(mapping) <- names(state_means)
    
    # Apply mapping to states
    ranked_states <- mapping[as.character(states)]
    
    # Return both ranked states and mapping
    list(
        ranked_states = ranked_states,
        mapping = mapping,
        state_means = state_means
    )
}

#' Apply ranking map to state assignments
#' @param states Vector of state assignments
#' @param mapping Named vector mapping old states to new states
#' @return Vector of mapped state assignments
#' @noRd
apply_state_ranking <- function(states, mapping) {
    # Convert states to character to use as index
    char_states <- as.character(states)
    
    # Apply mapping
    mapped_states <- mapping[char_states]
    
    return(mapped_states)
}

#' Create state summary with rankings
#'
#' This internal helper function generates a summary data frame for discrete states,
#' including counts, mean, median, standard deviation, and rankings based on mean and median values.
#'
#' @param states A vector of discrete state assignments for each data point.
#' @param values A numeric vector of the original time series values.
#' @return A data frame summarizing statistics for each unique state, including:
#'   \item{state}{The unique state identifier.}
#'   \item{count}{The number of data points assigned to this state.}
#'   \item{mean}{The mean of `values` for data points in this state.}
#'   \item{median}{The median of `values` for data points in this state.}
#'   \item{sd}{The standard deviation of `values` for data points in this state.}
#'   \item{rank_by_mean}{The rank of the state based on its mean value (lower rank for higher mean).}
#'   \item{rank_by_median}{The rank of the state based on its median value (lower rank for higher median).}
#' @keywords internal
create_state_summary <- function(states, values) {
    # Calculate statistics for each state
    state_summary <- data.frame(
        state = sort(unique(states)),
        count = as.vector(table(states)),
        mean = tapply(values, states, mean)[as.character(sort(unique(states)))],
        median = tapply(values, states, median)[as.character(sort(unique(states)))],
        sd = tapply(values, states, sd)[as.character(sort(unique(states)))]
    )
    
    # Add rankings
    state_summary$rank_by_mean <- rank(-state_summary$mean)  # Higher values get lower ranks
    state_summary$rank_by_median <- rank(-state_summary$median)
    
    return(state_summary)
}
