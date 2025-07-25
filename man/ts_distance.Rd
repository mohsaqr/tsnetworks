% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ts_distance.R
\name{ts_distance}
\alias{ts_distance}
\title{Time Series Distance Calculation}
\usage{
ts_distance(
  ts_data,
  method = "euclidean",
  window_size = NULL,
  step_size = NULL,
  pairwise = TRUE,
  symmetric = TRUE,
  ...
)
}
\arguments{
\item{ts_data}{Time series data (vector)}

\item{method}{Distance method ("euclidean", "manhattan", "cosine", "correlation", "dtw", etc.)}

\item{window_size}{Size of sliding windows (if NULL, uses full series)}

\item{step_size}{Step size for windows (defaults to window_size)}

\item{pairwise}{If TRUE, calculates all pairwise distances; if FALSE, consecutive windows only}

\item{symmetric}{If TRUE, makes matrix symmetric (both directions); if FALSE, keeps directional}

\item{...}{Additional arguments passed to distance calculation}
}
\value{
Distance matrix (always square matrix, regardless of pairwise setting)
}
\description{
Calculate distances between time series windows using various methods
}
