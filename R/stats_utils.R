#' @title Statistical Utility Functions
#' @description Helper functions for statistical calculations

#' Calculate Skewness
#' @param x Numeric vector
#' @return Skewness value
#' @noRd
.calculate_skewness <- function(x) {
    n <- length(x)
    if (n < 3) return(NA)
    m <- mean(x, na.rm = TRUE)
    s <- sd(x, na.rm = TRUE)
    if (s == 0) return(0)
    sum((x - m)^3) / ((n - 1) * s^3)
}

#' Calculate Kurtosis
#'
#' This internal helper function calculates the excess kurtosis of a numeric vector.
#' Kurtosis is a measure of the "tailedness" of the probability distribution of a real-valued random variable.
#'
#' @param x A numeric vector.
#' @return A numeric value representing the excess kurtosis. Returns `NA` if the vector has fewer than 4 elements or `0` if the standard deviation is zero.
#' @keywords internal
.calculate_kurtosis <- function(x) {
    n <- length(x)
    if (n < 4) return(NA)
    m <- mean(x, na.rm = TRUE)
    s <- sd(x, na.rm = TRUE)
    if (s == 0) return(0)
    (sum((x - m)^4) / ((n - 1) * s^4)) - 3
}

#' Calculate comprehensive statistics for states
#' @param states Vector of state assignments
#' @param values Original time series values
#' @param data Original data frame
#' @param value_column Name of value column
#' @return List containing state statistics and summary
#' @export
calculate_state_statistics <- function(states, values, data, value_column) {
    # Calculate statistics for each state
    state_stats <- lapply(sort(unique(states)), function(s) {
        state_values <- values[states == s]
        state_indices <- which(states == s)
        
        list(
            n = length(state_values),
            frequency = length(state_values) / length(values),
            mean = mean(state_values),
            median = median(state_values),
            sd = sd(state_values),
            min = min(state_values),
            max = max(state_values),
            q25 = quantile(state_values, 0.25),
            q75 = quantile(state_values, 0.75),
            skewness = .calculate_skewness(state_values),
            kurtosis = .calculate_kurtosis(state_values),
            indices = state_indices
        )
    })
    
    # Safe naming of state_stats
    unique_states <- sort(unique(states))
    if (length(unique_states) == 0) {
        warning("No unique states found in calculate_state_statistics")
        return(list(state_stats = state_stats, summary = data.frame()))
    }
    if (length(state_stats) == length(unique_states)) {
        names(state_stats) <- paste0("state_", unique_states)
    } else {
        warning(sprintf("Mismatch: state_stats length (%d) vs unique_states (%d)",
                       length(state_stats), length(unique_states)))
    }
    
    # Create summary data frame
    summary_df <- do.call(rbind, lapply(names(state_stats), function(state_name) {
        stats <- state_stats[[state_name]]
        data.frame(
            state = gsub("state_", "", state_name),
            n = stats$n,
            frequency = stats$frequency,
            mean = stats$mean,
            median = stats$median,
            sd = stats$sd,
            min = stats$min,
            max = stats$max,
            q25 = stats$q25,
            q75 = stats$q75,
            skewness = stats$skewness,
            kurtosis = stats$kurtosis
        )
    }))
    
    list(
        detailed = state_stats,
        summary = summary_df
    )
}

#' Transform time series data
#' @param series Numeric vector of time series values
#' @param method Transformation method ("none", "log", "sqrt", "scale", "center")
#' @return List with transformed series and transformation info
#' @export
transform_series <- function(series, method = "none") {
    result <- list(
        original = series,
        transformed = series,
        method = method,
        parameters = list()
    )
    
    if (method == "log") {
        # Handle negative values
        min_val <- min(series)
        if (min_val <= 0) {
            shift <- abs(min_val) + 1
            result$transformed <- log(series + shift)
            result$parameters$shift <- shift
        } else {
            result$transformed <- log(series)
        }
    } else if (method == "sqrt") {
        # Handle negative values
        min_val <- min(series)
        if (min_val < 0) {
            shift <- abs(min_val)
            result$transformed <- sqrt(series + shift)
            result$parameters$shift <- shift
        } else {
            result$transformed <- sqrt(series)
        }
    } else if (method == "scale") {
        result$transformed <- scale(series)
        result$parameters$center <- attr(result$transformed, "scaled:center")
        result$parameters$scale <- attr(result$transformed, "scaled:scale")
        # Convert from matrix to vector
        result$transformed <- as.vector(result$transformed)
    } else if (method == "center") {
        result$transformed <- scale(series, scale = FALSE)
        result$parameters$center <- attr(result$transformed, "scaled:center")
        # Convert from matrix to vector
        result$transformed <- as.vector(result$transformed)
    }
    
    return(result)
}

#' Reverse transform time series data
#' @param transformed Transformed series
#' @param method Original transformation method
#' @param parameters Parameters used in the original transformation
#' @return Original series
#' @export
reverse_transform <- function(transformed, method, parameters) {
    if (method == "log") {
        if (!is.null(parameters$shift)) {
            exp(transformed) - parameters$shift
        } else {
            exp(transformed)
        }
    } else if (method == "sqrt") {
        if (!is.null(parameters$shift)) {
            transformed^2 - parameters$shift
        } else {
            transformed^2
        }
    } else if (method == "scale") {
        transformed * parameters$scale + parameters$center
    } else if (method == "center") {
        transformed + parameters$center
    } else {
        transformed
    }
}
