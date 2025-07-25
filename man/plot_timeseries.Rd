% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot_timeseries}
\alias{plot_timeseries}
\title{Plot Time Series with State Mapping}
\usage{
plot_timeseries(
  data,
  id_col = NULL,
  value_col = "value",
  time_col = NULL,
  state_col = "state",
  selected = NULL,
  overlay = "v",
  points = FALSE,
  ncol = NULL,
  max_series = 10,
  trend = FALSE,
  scales = "free"
)
}
\arguments{
\item{data}{A data.frame containing time series data}

\item{id_col}{Character name of the ID column (optional, creates single series if NULL)}

\item{value_col}{Character name of the value column}

\item{time_col}{Character name of the time column (optional, uses row order if NULL)}

\item{state_col}{Character name of the state column}

\item{selected}{Character vector of specific series IDs to plot}

\item{overlay}{Character string: "v" (vertical), "h" (horizontal), or NULL (no overlay)}

\item{points}{Logical, whether to show points colored by state}

\item{ncol}{Integer, number of columns for facets}

\item{max_series}{Integer, maximum number of series to plot}

\item{trend}{Logical, whether to add trend line}

\item{scales}{Character, facet scales: "fixed", "free", "free_x", "free_y"}
}
\value{
A ggplot object
}
\description{
A robust, self-contained function to visualize time series data with state overlays
}
