#' @importFrom stats na.omit

#' @title Transition Matrix Utilities
#' @description Functions for calculating and manipulating transition matrices

#' Calculate Transition Counts
#' @param states Vector of state assignments
#' @param num_states Total number of possible states (optional)
#' @return Matrix of transition counts
#' @noRd
.calculate_transition_counts <- function(states, num_states = NULL) {
    states <- na.omit(states)
    if (length(states) < 2) 
        stop("Need at least 2 states to calculate transitions")
    
    if (is.null(num_states)) {
        unique_states <- sort(unique(states))
        n_states <- length(unique_states)
        state_names <- unique_states
    } else {
        n_states <- num_states
        state_names <- 1:num_states
    }
    
    # Initialize transition matrix
    trans_matrix <- matrix(0, 
                          nrow = n_states, 
                          ncol = n_states,
                          dimnames = list(state_names, state_names))
    
    # Count transitions
    for (i in 1:(length(states) - 1)) {
        if (!is.na(states[i]) && !is.na(states[i + 1])) {
            trans_matrix[states[i], states[i + 1]] <- 
                trans_matrix[states[i], states[i + 1]] + 1
        }
    }
    trans_matrix
}

#' Calculate transitions with fixed number of states
#' @param states Vector of state assignments
#' @param num_states Total number of possible states
#' @return Matrix of transition counts
#' @export
calculate_transitions <- function(states, num_states = NULL) {
    .calculate_transition_counts(states, num_states)
}

#' Normalize Transition Counts to Probabilities
#' @param counts Transition count matrix
#' @return Matrix of transition probabilities
#' @noRd
.normalize_transition_counts <- function(counts) {
    # Row-wise normalization
    row_sums <- rowSums(counts)
    probs <- counts / row_sums
    # Handle rows with all zeros
    probs[row_sums == 0, ] <- 0
    probs
}

#' Create Network Data Structure
#' @param trans_matrix Transition probability matrix
#' @return List containing network data
#' @noRd
.create_network_data <- function(trans_matrix) {
    # Extract edges and weights
    edges <- which(trans_matrix > 0, arr.ind = TRUE)
    weights <- trans_matrix[edges]
    
    list(
        nodes = rownames(trans_matrix),
        edges = edges,
        weights = weights,
        transition_matrix = trans_matrix
    )
}

#' Normalize transition counts to probabilities
#' @param transitions Matrix of transition counts
#' @return Matrix of transition probabilities
#' @noRd
normalize_transitions <- function(transitions) {
    # Add small constant to avoid division by zero
    epsilon <- 1e-10
    row_sums <- rowSums(transitions) + epsilon
    
    # Normalize each row to get probabilities
    prob_matrix <- sweep(transitions, 1, row_sums, FUN = "/")
    
    return(prob_matrix)
}

#' Create Network Data Structure from Transition Matrix
#'
#' This internal helper function converts a transition probability matrix into
#' a list structure suitable for network analysis or visualization, including
#' nodes, edges, and weights.
#'
#' @param transition_matrix A square matrix of transition probabilities.
#' @return A list containing:
#'   \item{nodes}{A vector of node (state) identifiers.}
#'   \item{edges}{A matrix where each row represents an edge, with columns for
#'     `from` state index and `to` state index.}
#'   \item{weights}{A vector of edge weights (transition probabilities) corresponding
#'     to the `edges` matrix.}
#'   \item{transition_matrix}{The input transition probability matrix.}
#' @keywords internal
#' @noRd
create_network_data <- function(transition_matrix) {
    n_states <- nrow(transition_matrix)
    
    # Create edge list
    edges <- which(transition_matrix > 0, arr.ind = TRUE)
    weights <- transition_matrix[edges]
    
    # Create network data structure
    list(
        nodes = 1:n_states,
        edges = edges,
        weights = weights,
        transition_matrix = transition_matrix
    )
}
