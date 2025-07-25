% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visibility_graph.R
\name{visibility_graph}
\alias{visibility_graph}
\title{Compute visibility graph from time series}
\usage{
visibility_graph(
  data,
  ts_value = NULL,
  state_value = NULL,
  method = "nvg",
  directed = FALSE,
  limit = NULL,
  penetrable = 0,
  decay_factor = 0,
  state_map = NULL,
  weight_agg = "sum",
  output_format = "matrix",
  default_view = NULL,
  num_cores = 1
)
}
\arguments{
\item{data}{A data frame or numeric vector containing the time series data.
If a data frame, the \code{ts_value} parameter must specify which column
contains the time series values.}

\item{ts_value}{Character string or NULL. If \code{data} is a data frame, this
specifies the column name containing time series values. If \code{data} is 
a numeric vector, this parameter is ignored.}

\item{state_value}{Character string, vector, or NULL. If \code{data} is a data frame,
this should be a column name containing state labels. If \code{data} is a vector,
this should be a vector of the same length containing state labels. If provided,
the output will include state-based graph representations.}

\item{method}{Character string specifying the visibility graph construction method:
\itemize{
  \item \code{"nvg"} (default): Natural Visibility Graph
  \item \code{"hvg"}: Horizontal Visibility Graph
}}

\item{directed}{Logical value indicating whether to create a directed graph.
If \code{TRUE}, edges have direction showing temporal order. Default is \code{FALSE}.}

\item{limit}{Positive integer or NULL. Maximum temporal distance (in time steps)
for visibility connections. If NULL (default), no limit is applied.}

\item{penetrable}{Non-negative integer. Number of points allowed to penetrate
the visibility line for Limited Penetrable Visibility Graphs (LPVG).
Default is 0 (standard visibility graph).}

\item{decay_factor}{Non-negative numeric value controlling temporal decay of
edge weights. Default is 0 (no decay). Higher values cause faster decay
with increasing temporal distance.}

\item{state_map}{State mapping specification used when \code{state_value} is NULL.
Can be:
\itemize{
  \item NULL (default): No state mapping
  \item Named vector: Maps value ranges to states (e.g., c("0-25"="low", ">50"="high"))
  \item Function: Takes numeric values and returns state labels
  \item Numeric vector: Break points for equal-width bins
  \item Single number: Number of equal-width bins to create
}}

\item{weight_agg}{Character string specifying how to aggregate edge weights
when creating state-based graphs:
\itemize{
  \item \code{"sum"} (default): Sum all weights between state pairs
  \item \code{"mean"}: Average of weights between state pairs
  \item \code{"max"}: Maximum weight between state pairs
  \item \code{"count"}: Count of connections between state pairs
}}

\item{output_format}{Character string specifying the output format:
\code{"matrix"} (default) or \code{"edgelist"}. This affects the primary
output slots but all representations are always included.}

\item{default_view}{Character string specifying which graph view to use
for the primary output slots:
\itemize{
  \item \code{"time"}: Time-point based graph (default when no states)
  \item \code{"state"}: State-based graph (default when states provided)
  \item \code{"value"}: Value-based graph (groups by unique values)
}}

\item{num_cores}{Positive integer specifying number of CPU cores to use
for parallel processing. Default is 1 (no parallelization). Parallelization
is only used for time series with more than 100 points.}
}
\value{
A list containing multiple graph representations:
  \describe{
    \item{\code{matrix}}{Primary adjacency matrix (format determined by \code{default_view})}
    \item{\code{edge_list}}{Primary edge list (format determined by \code{default_view})}
    \item{\code{time_matrix}}{Time-point based adjacency matrix}
    \item{\code{time_edge_list}}{Time-point based edge list}
    \item{\code{values_matrix}}{Value-based adjacency matrix (groups by unique values)}
    \item{\code{values_edge_list}}{Value-based edge list}
    \item{\code{state_matrix}}{State-based adjacency matrix (if states provided)}
    \item{\code{state_edge_list}}{State-based edge list (if states provided)}
    \item{\code{original_values}}{Original time series values}
    \item{\code{states}}{State labels (if provided)}
    \item{\code{method}}{Visibility method used}
    \item{\code{directed}}{Whether the graph is directed}
    \item{\code{parameters}}{List of all parameters used}
  }
}
\description{
Creates visibility graphs from time series data, where nodes represent time points
or states, and edges represent visibility between points according to the selected method.
The function supports both Natural Visibility Graphs (NVG) and Horizontal Visibility 
Graphs (HVG), with options for limited penetrable visibility, temporal decay, and 
state-based aggregation.
}
\details{
\subsection{Visibility Graph Methods}{
\strong{Natural Visibility Graph (NVG):} Two points can "see" each other if
a straight line connecting them doesn't intersect with any intermediate points.
This preserves more geometric information from the original time series.

\strong{Horizontal Visibility Graph (HVG):} Two points can "see" each other if
all intermediate points are strictly lower than both endpoints. This creates
sparser graphs that capture different temporal patterns.
}

\subsection{Advanced Features}{
\strong{Limited Penetrable Visibility (LPVG):} Allows a specified number of
points to penetrate the visibility line, creating more connected graphs.

\strong{Temporal Decay:} Edge weights decrease exponentially with temporal
distance, emphasizing local temporal relationships.

\strong{State-based Analysis:} Converts continuous time series to discrete
states for categorical analysis of temporal patterns.
}
}
\section{Performance Considerations}{

\itemize{
  \item Time complexity is O(n²) for the basic algorithm
  \item Memory usage scales with the number of edges created
  \item Parallel processing is automatically used for n > 100 when num_cores > 1
  \item State-based graphs can significantly reduce complexity for large time series
}
}

\examples{
\dontrun{
# Basic usage with a simple time series
set.seed(123)
ts_data <- rnorm(50, mean = 50, sd = 15)

# Create natural visibility graph
nvg <- visibility_graph(ts_data, method = "nvg")
print(dim(nvg$matrix))  # Time-based adjacency matrix

# Create horizontal visibility graph with temporal limit
hvg <- visibility_graph(ts_data, method = "hvg", limit = 10)

# Using data frame input
df <- data.frame(
  time = 1:20,
  value = c(5, 15, 25, 75, 45, 10, 30, 20, 60, 35,
            40, 80, 55, 25, 70, 45, 35, 65, 50, 30),
  category = rep(c("A", "B"), each = 10)
)

# Create state-based visibility graph
state_vg <- visibility_graph(df, ts_value = "value", state_value = "category")
print(state_vg$state_matrix)  # State-based adjacency matrix

# Using quantile-based states
quantile_states <- compute_quantile_states(ts_data, num_quantiles = 4)
qvg <- visibility_graph(ts_data, state_value = quantile_states)

# Advanced features
lpvg <- visibility_graph(ts_data, penetrable = 2, decay_factor = 0.1)

# Custom state mapping
custom_map <- c("0-25" = "low", "25-75" = "medium", "75-100" = "high")
custom_vg <- visibility_graph(ts_data, state_map = custom_map)

# Parallel processing for large time series
large_ts <- rnorm(1000)
large_vg <- visibility_graph(large_ts, num_cores = 4)
}

}
\references{
\itemize{
  \item Lacasa, L., Luque, B., Ballesteros, F., Luque, J., & Nuño, J. C. (2008).
    From time series to complex networks: The visibility graph. PNAS, 105(13), 4972-4975.
  \item Luque, B., Lacasa, L., Ballesteros, F., & Luque, J. (2009).
    Horizontal visibility graphs: Exact results for random time series. Physical Review E, 80(4), 046103.
}
}
\seealso{
\code{\link{compute_quantile_states}} for creating quantile-based states
}
