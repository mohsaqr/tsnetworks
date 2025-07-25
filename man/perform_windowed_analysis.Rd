% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculator.R
\name{perform_windowed_analysis}
\alias{perform_windowed_analysis}
\title{Perform windowed analysis}
\usage{
perform_windowed_analysis(
  time_series,
  method,
  window_size,
  step_size,
  output_format,
  pairwise,
  ...
)
}
\arguments{
\item{time_series}{The time series data.}

\item{method}{The distance calculation method.}

\item{window_size}{The size of the sliding windows.}

\item{step_size}{The step size for the windows.}

\item{output_format}{The format of the output, either "matrix" or "edges".}

\item{pairwise}{If TRUE, all pairwise distances are calculated; if FALSE, only consecutive windows are compared.}

\item{...}{Additional arguments passed to the distance functions.}
}
\value{
A list containing the results of the windowed analysis.
}
\description{
Perform windowed analysis
}
