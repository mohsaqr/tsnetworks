% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visibility_graph.R
\name{compute_quantile_states}
\alias{compute_quantile_states}
\title{Compute quantile-based states for time series}
\usage{
compute_quantile_states(
  x,
  num_quantiles = 4,
  quantile_probs = NULL,
  labels = NULL,
  na_label = "NA",
  include_lowest = TRUE,
  right = FALSE
)
}
\arguments{
\item{x}{Numeric vector of time series values}

\item{num_quantiles}{Integer. Number of quantile-based states to create. Must be >= 2.}

\item{quantile_probs}{Numeric vector of quantile probabilities. If provided, 
\code{num_quantiles} is ignored. Values must be between 0 and 1 and in 
ascending order.}

\item{labels}{Character vector of state labels. If NULL, automatic labels are generated.
Length must match the number of states created.}

\item{na_label}{Character string to use for NA values. Default is "NA".}

\item{include_lowest}{Logical. Should the lowest break be included? Passed to \code{cut}.}

\item{right}{Logical. Should intervals be closed on the right? Passed to \code{cut}.}
}
\value{
Character vector of the same length as \code{x} containing state labels
}
\description{
This function maps time series values to discrete states based on quantiles.
It's useful for creating categorical representations of continuous time series data.
}
\details{
The function creates states by dividing the time series into quantile-based bins.
For example, with \code{num_quantiles = 4}, the function creates quartile-based
states: Q1 (0-25%), Q2 (25-50%), Q3 (50-75%), and Q4 (75-100%).
}
\examples{
\dontrun{
# Generate sample data
set.seed(123)
ts_data <- rnorm(100, mean = 50, sd = 15)

# Create quartile-based states
states_4 <- compute_quantile_states(ts_data, num_quantiles = 4)
table(states_4)

# Create tercile states with custom labels
states_3 <- compute_quantile_states(ts_data, num_quantiles = 3, 
                                   labels = c("Low", "Medium", "High"))
table(states_3)

# Use custom quantile probabilities
states_custom <- compute_quantile_states(ts_data, 
                                        quantile_probs = c(0, 0.1, 0.9, 1.0),
                                        labels = c("Very_Low", "Normal", "Very_High"))
table(states_custom)
}

}
