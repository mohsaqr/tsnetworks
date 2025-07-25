% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visibility_graph.R
\name{map_to_state}
\alias{map_to_state}
\title{Map time series values to discrete states}
\usage{
map_to_state(x, state_map)
}
\arguments{
\item{x}{Numeric vector of time series values}

\item{state_map}{State mapping specification. Can be:
\itemize{
  \item A named vector mapping ranges to states (e.g., c("0-10"="low", "11-50"="medium"))
  \item A function that takes a numeric value and returns a state
  \item A numeric vector specifying break points for equal-width bins
  \item A single number specifying the number of equal-width bins to create
}}
}
\value{
Character vector of state values corresponding to time series values
}
\description{
This function converts continuous time series values to discrete states
using various mapping strategies.
}
\details{
The function supports several mapping strategies:
\itemize{
  \item \strong{Named vector}: Use range strings like "0-10", "<5", ">100", "=50"
  \item \strong{Function}: Custom mapping function
  \item \strong{Numeric breaks}: Explicit break points for binning
  \item \strong{Number of bins}: Creates equal-width bins
}
}
\keyword{internal}
