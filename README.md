# tsnetworks: An R Package for Time Series Network Analysis

`tsnetworks` is a comprehensive R package that provides a suite of tools for analyzing time series data through network-based approaches. It integrates three powerful methodologies:

1.  **Visibility Graphs (VG):** Construct networks based on the line-of-sight visibility between time series data points.
2.  **Dynamic Complexity and Regime Detection:** Use rolling window metrics to calculate the complexity of a time series and identify distinct regimes or state changes.
3.  **Time Series Proximity Networks:** Build networks based on the similarity between sliding windows of a time series, ideal for pattern discovery.

This package is designed to be a one-stop solution for researchers and analysts looking to explore the structural and dynamic properties of time series data.

## Installation

To install the package from the local source file, run the following command in your R console:

```r
# Make sure your working directory is set to where the .tar.gz file is located
install.packages("tsnetworks_1.0.0.tar.gz", repos = NULL, type = "source")
```

## Core Features

- **Multiple Network Models:** Go beyond simple correlations with Visibility Graphs, Proximity Networks, and Markov Transition Networks.
- **Rich Discretization Toolkit:** Convert continuous time series into discrete states using a wide variety of methods, including kmeans, quantiles, entropy, and Gaussian mixture models.
- **Dynamic Analysis:** Analyze how the complexity and structure of a time series evolve over time with rolling window-based measures.
- **Regime Detection:** Automatically identify and classify different states or regimes within your time series data.
- **Integrated Visualization:** A suite of plotting functions to visualize raw time series, regime changes, and network structures.

---

## Main Functions

The package is organized around its three core methodologies.

### 1. Discretization and State Analysis (`stna`)

This is the main entry point for discretizing a time series into states and analyzing the transitions between them.

- `stna(time_series, value_column, num_states, method, ...)`: The primary function. It takes a time series and a discretization method to produce state assignments, a Markov transition matrix, and detailed statistics.

#### Key Discretization Methods for `stna()`

- `kmeans`: Partitions data into *k* clusters.
- `quantile`: Bins data based on quantile ranges.
- `equal_width`: Bins data into intervals of equal width.
- `entropy`: Uses entropy maximization to find optimal bin widths.
- `mixture`: Uses Gaussian Mixture Models to identify underlying distributions as states.
- `hierarchical`: Uses hierarchical clustering to create states.
- `change_points`: Identifies states based on the largest shifts in the data.
- `proxy_windowed`: A powerful method that performs clustering on sliding windows of the time series using any distance metric from the `proxy` package.

### 2. Visibility Graphs

- `visibility_graph(data, ...)`: Constructs a Natural Visibility Graph (NVG) or Horizontal Visibility Graph (HVG) from a time series. Supports advanced features like limited penetrable visibility and temporal decay.
- `compute_quantile_states(x, ...)`: A utility to discretize a numeric vector into states based on quantiles, often used as a precursor to state-based graph analysis.

### 3. Dynamic Complexity & Regime Detection

- `rolling_measures(data, ...)`: Calculates various rolling window metrics, including `complexity`, `fluctuation`, and `distribution`, for one or more time series.
- `detect_regime(complexity_data, ...)`: Takes the output of `rolling_measures` (or any complexity score) and identifies distinct regimes using methods like `smart` (a combination of techniques), `changepoint`, `variance_shift`, and more.

### 4. Time Series Proximity Networks

- `ts_distance(ts_data, ...)`: The first step in this workflow. It creates a distance matrix by calculating the distance between sliding windows of a time series.
- `ts_network(distance_matrix, ...)`: Takes the distance matrix and builds a network. You can control the network density using methods like `knn` (k-nearest neighbors) or `percentile` thresholding.

### 5. Visualization

- `plot_timeseries(data, ...)`: A flexible `ggplot2`-based function to plot a time series with states shaded as vertical or horizontal bands.
- `create_combined_plots(data, ...)`: Creates a two-panel plot showing the original time series and its complexity scores, with regimes shaded for easy comparison.
- `plot_tsn(network_object, ...)`: A wrapper for `qgraph` to easily visualize the proximity networks created by `ts_network()`.

---

## Quick Start Example: A Full Workflow

This example shows a common workflow: calculating dynamic complexity, detecting regimes, and visualizing the results.

```r
# Load the package
library(tsnetworks)

# 1. Load or create your data
# We will use the built-in 'saqrsteps' dataset for this example
data(saqrsteps)

# 2. Calculate rolling complexity
# We will analyze the 'Steps' column with a window of 7 days
complexity_results <- rolling_measures(
  data = saqrsteps,
  ts_cols = "Steps",
  measures = c("complexity", "fluctuation"),
  window_width = 7
)

# 3. Detect regimes from the complexity scores
# The function automatically finds the complexity column and preserves all original data
regime_data <- detect_regime(
  complexity_results,
  method = "smart",
  sensitivity = "medium"
)

# 4. Visualize the results
# This plot shows the original time series and its complexity, with regimes shaded
create_combined_plots(
  data = regime_data,
  original_col = "Steps",
  complexity_col = "complexity_Steps", # The column created by rolling_measures
  state_col = "regime_stability"      # The column created by detect_regime
)

# 5. Explore the regime summary
# See how many data points fall into each stability category
print(table(regime_data$regime_stability))

```
