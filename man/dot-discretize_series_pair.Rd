% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/discretization.R
\name{.discretize_series_pair}
\alias{.discretize_series_pair}
\title{Discretize a pair of time series}
\usage{
.discretize_series_pair(series1, series2, n_bins, method = "quantile")
}
\arguments{
\item{series1}{Numeric vector. The first time series.}

\item{series2}{Numeric vector. The second time series.}

\item{n_bins}{Integer. The number of bins to use for discretization.}

\item{method}{Character string. The discretization method to use (e.g., "quantile", "equal_width").}
}
\value{
A list containing two discretized numeric vectors, `d1` and `d2`.
}
\description{
This internal helper function discretizes two time series simultaneously
using a specified method. This is useful for information-theoretic measures
that require joint discretization.
}
\keyword{internal}
