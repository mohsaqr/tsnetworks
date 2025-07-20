# AI Coder Plan: Add Resilience Metrics to TSnetwork Package

## Project Overview
Extend the existing `TSnetwork` package by adding comprehensive resilience metrics functionality. The new features will integrate seamlessly with TSnetwork's existing architecture, utilize its utilities, and follow its design patterns for time series analysis and state classification.

## Pre-Implementation Analysis

### Task 0.1: Analyze TSnetwork Package Structure
**Required Actions:**
1. Study TSnetwork's existing codebase to understand:
   - Current function naming conventions
   - Data structure patterns (how time series are handled)
   - State classification methodologies already implemented
   - Utility functions available for reuse
   - Documentation style and roxygen2 patterns
   - Testing framework structure

2. Identify existing functions to avoid duplication:
   - Time series preprocessing utilities
   - Rolling window calculations
   - State classification frameworks
   - Visualization functions
   - Data validation routines

3. Map TSnetwork's current capabilities to resilience metrics needs:
   - Which resilience metrics complement existing network analysis
   - How to integrate with existing state definitions
   - Where resilience analysis fits in the current workflow

## Implementation Plan

### Phase 1: Integration Planning and Setup

#### Task 1.1: Create Resilience Module Structure
**File Structure within TSnetwork:**
```
R/
├── resilience_metrics.R      # Core resilience metric calculations
├── resilience_capacity.R     # Absorptive, restorative, adaptive functions
├── resilience_classification.R # State classification for resilience
└── resilience_utils.R        # Utilities specific to resilience analysis
```

#### Task 1.2: Extend TSnetwork's Existing Data Structures
**Modify existing TSnetwork data objects to support resilience analysis:**

```r
# Extend existing TSnetwork functions to handle resilience metrics
# Follow TSnetwork's existing pattern for data handling

#' Add Resilience Analysis to TSnetwork Object
#' 
#' @param ts_object Existing TSnetwork time series object
#' @param capacity_types Vector of capacity types to analyze c("absorptive", "restorative", "adaptive")
#' @param window_size Rolling window size (inherit from TSnetwork defaults)
#' @return Enhanced TSnetwork object with resilience metrics
#' @export
add_resilience_analysis <- function(ts_object, 
                                   capacity_types = c("absorptive", "restorative", "adaptive"),
                                   window_size = NULL) {
  
  # Use TSnetwork's existing data validation
  # Leverage TSnetwork's window size defaults
  # Follow TSnetwork's object structure patterns
}
```

### Phase 2: Core Resilience Metrics Integration

#### Task 2.1: Absorptive Capacity Metrics (Leverage TSnetwork Utils)
**File: `R/resilience_capacity.R`**

```r
#' Calculate Variance Stability Index
#' 
#' Integrates with TSnetwork's existing rolling window functions
#' @param ts_data TSnetwork time series data object
#' @param baseline_method Method for baseline calculation (use TSnetwork's existing methods)
#' @return Variance stability index values
#' @export
#' @family resilience_metrics
calculate_variance_stability <- function(ts_data, baseline_method = "auto") {
  
  # Leverage TSnetwork's existing rolling window utilities
  # Use TSnetwork's baseline calculation methods if they exist
  # Follow TSnetwork's naming convention for metric functions
  
  # Get baseline variance using TSnetwork's baseline functions
  baseline_var <- get_baseline_variance(ts_data, method = baseline_method)
  
  # Use TSnetwork's rolling calculation utilities
  rolling_var <- apply_rolling_function(ts_data, 
                                       function(x) var(x, na.rm = TRUE),
                                       window_size = get_default_window(ts_data))
  
  # Calculate VSI following TSnetwork's metric calculation pattern
  vsi <- 1 - abs(rolling_var - baseline_var) / baseline_var
  vsi[vsi < 0] <- 0
  
  # Return in TSnetwork's standard metric format
  return(create_metric_object(vsi, metric_name = "variance_stability", 
                             ts_data = ts_data))
}

#' Calculate Resilience ARCH Effects
#' 
#' Extends TSnetwork's volatility analysis capabilities
#' @param ts_data TSnetwork time series data object
#' @return ARCH-based resilience metric
#' @export
#' @family resilience_metrics
calculate_resilience_arch <- function(ts_data) {
  
  # Check if TSnetwork already has ARCH analysis
  if (has_volatility_analysis(ts_data)) {
    # Extend existing volatility metrics for resilience
    volatility_metrics <- get_volatility_metrics(ts_data)
    # Adapt for resilience context
  } else {
    # Implement ARCH test using TSnetwork's regression utilities
  }
  
  # Use TSnetwork's existing AR model fitting if available
  ar_residuals <- get_ar_residuals(ts_data)
  
  # Calculate ARCH statistic using TSnetwork's regression functions
  arch_stat <- calculate_heteroscedasticity(ar_residuals)
  
  return(create_metric_object(arch_stat, metric_name = "resilience_arch", 
                             ts_data = ts_data))
}

#' Calculate Adaptive Coefficient of Variation
#' 
#' Extends TSnetwork's existing CV calculations for resilience context
#' @param ts_data TSnetwork time series data object
#' @return Adaptive coefficient of variation
#' @export
#' @family resilience_metrics
calculate_adaptive_cv <- function(ts_data) {
  
  # Check if TSnetwork already calculates CV
  if (has_coefficient_variation(ts_data)) {
    # Enhance existing CV for resilience analysis
    existing_cv <- get_coefficient_variation(ts_data)
    # Add resilience-specific adaptations
  } else {
    # Calculate using TSnetwork's rolling statistics utilities
    rolling_cv <- apply_rolling_function(ts_data,
                                        function(x) sd(x, na.rm = TRUE) / abs(mean(x, na.rm = TRUE)),
                                        window_size = get_adaptive_window(ts_data))
  }
  
  return(create_metric_object(rolling_cv, metric_name = "adaptive_cv", 
                             ts_data = ts_data))
}
```

#### Task 2.2: Restorative Capacity Metrics (Integrate with TSnetwork's Event Detection)
**File: `R/resilience_capacity.R` (continued)**

```r
#' Calculate Recovery Time Metrics
#' 
#' Leverages TSnetwork's existing change point detection and event analysis
#' @param ts_data TSnetwork time series data object
#' @param event_detection_method Use TSnetwork's existing event detection
#' @return Recovery time metrics
#' @export
#' @family resilience_metrics
calculate_recovery_metrics <- function(ts_data, event_detection_method = "auto") {
  
  # Use TSnetwork's existing event/disruption detection if available
  if (has_event_detection(ts_data)) {
    disruption_events <- get_detected_events(ts_data)
  } else {
    # Use TSnetwork's change point detection utilities
    disruption_events <- detect_change_points(ts_data, method = event_detection_method)
  }
  
  recovery_metrics <- list()
  
  # Calculate return time using TSnetwork's AR modeling utilities
  if (has_ar_modeling(ts_data)) {
    ar_params <- get_ar_parameters(ts_data)
    recovery_metrics$return_time <- calculate_return_time_from_ar(ar_params)
  }
  
  # Calculate recovery slope using TSnetwork's trend analysis
  if (has_trend_analysis(ts_data)) {
    recovery_metrics$recovery_slope <- calculate_recovery_slope_with_trends(ts_data, disruption_events)
  }
  
  # Use TSnetwork's exponential fitting utilities for half-life
  recovery_metrics$half_life <- calculate_half_life_recovery(ts_data, disruption_events)
  
  return(create_metric_object(recovery_metrics, metric_name = "recovery_metrics", 
                             ts_data = ts_data))
}

#' Calculate Engineering Resilience
#' 
#' Integrates with TSnetwork's stability analysis
#' @param recovery_time_metric Output from calculate_recovery_metrics
#' @return Engineering resilience values
#' @export
#' @family resilience_metrics
calculate_engineering_resilience <- function(recovery_time_metric) {
  
  # Use TSnetwork's metric transformation utilities
  engineering_resilience <- transform_metric(recovery_time_metric,
                                           transformation = function(x) 1/x,
                                           metric_name = "engineering_resilience")
  
  return(engineering_resilience)
}
```

#### Task 2.3: Adaptive Capacity Metrics (Extend TSnetwork's Complexity Analysis)
**File: `R/resilience_capacity.R` (continued)**

```r
#' Calculate Resilience Sample Entropy
#' 
#' Extends TSnetwork's existing entropy/complexity analysis
#' @param ts_data TSnetwork time series data object
#' @param pattern_length Pattern length for entropy calculation
#' @return Sample entropy for resilience analysis
#' @export
#' @family resilience_metrics
calculate_resilience_sample_entropy <- function(ts_data, pattern_length = 2) {
  
  # Check if TSnetwork has existing entropy calculations
  if (has_entropy_analysis(ts_data)) {
    # Extend existing entropy metrics for resilience
    base_entropy <- get_entropy_metrics(ts_data)
    # Adapt parameters for resilience analysis
  }
  
  # Use TSnetwork's pattern matching utilities if available
  tolerance <- get_adaptive_tolerance(ts_data)
  
  # Calculate sample entropy using TSnetwork's rolling analysis framework
  sampen_values <- apply_rolling_function(ts_data,
                                         function(x) calculate_sample_entropy_core(x, pattern_length, tolerance),
                                         window_size = get_entropy_window(ts_data))
  
  return(create_metric_object(sampen_values, metric_name = "resilience_sample_entropy", 
                             ts_data = ts_data))
}

#' Calculate DFA for Resilience Analysis
#' 
#' Leverages TSnetwork's existing fractal/scaling analysis if available
#' @param ts_data TSnetwork time series data object
#' @return DFA alpha values for adaptive capacity
#' @export
#' @family resilience_metrics
calculate_resilience_dfa <- function(ts_data) {
  
  # Check if TSnetwork has existing DFA or scaling analysis
  if (has_scaling_analysis(ts_data)) {
    # Use existing scaling analysis infrastructure
    scaling_metrics <- get_scaling_metrics(ts_data)
    # Adapt for resilience context
  }
  
  # Use TSnetwork's detrending utilities
  dfa_values <- apply_rolling_function(ts_data,
                                      function(x) calculate_dfa_alpha(x),
                                      window_size = get_dfa_window(ts_data))
  
  return(create_metric_object(dfa_values, metric_name = "resilience_dfa", 
                             ts_data = ts_data))
}

#' Calculate Adaptive Capacity Ratio
#' 
#' Integrates with TSnetwork's response analysis framework
#' @param ts_data TSnetwork time series data object
#' @param external_driver External driver data (optional)
#' @return Adaptive capacity ratio
#' @export
#' @family resilience_metrics
calculate_adaptive_capacity_ratio <- function(ts_data, external_driver = NULL) {
  
  # Use TSnetwork's existing driver-response analysis if available
  if (has_driver_analysis(ts_data) && !is.null(external_driver)) {
    # Leverage existing driver-response framework
    response_variance <- get_response_variance(ts_data)
    driver_variance <- get_driver_variance(external_driver)
  } else {
    # Use first differences as proxy for disturbance
    disturbance_proxy <- calculate_first_differences(ts_data)
    response_variance <- apply_rolling_function(ts_data, var)
    driver_variance <- apply_rolling_function(disturbance_proxy, var)
  }
  
  ac_ratio <- response_variance / (driver_variance + get_small_constant())
  
  return(create_metric_object(ac_ratio, metric_name = "adaptive_capacity_ratio", 
                             ts_data = ts_data))
}
```

### Phase 3: State Classification Integration

#### Task 3.1: Extend TSnetwork's State Classification System
**File: `R/resilience_classification.R`**

```r
#' Classify Resilience States in TSnetwork Framework
#' 
#' Extends TSnetwork's existing state classification with resilience-specific states
#' @param ts_data TSnetwork time series object with calculated resilience metrics
#' @param classification_method Method compatible with TSnetwork's classification system
#' @return TSnetwork object with added resilience state classifications
#' @export
classify_resilience_states_tsnetwork <- function(ts_data, classification_method = "threshold") {
  
  # Check existing TSnetwork state classifications
  existing_states <- get_existing_states(ts_data)
  
  # Get all resilience metrics calculated for this time series
  resilience_metrics <- get_resilience_metrics(ts_data)
  
  # Calculate composite capacity scores using TSnetwork's scoring framework
  capacity_scores <- calculate_capacity_scores_tsnetwork(resilience_metrics)
  
  # Use TSnetwork's classification utilities
  resilience_states <- apply_classification_rules(
    capacity_scores, 
    method = classification_method,
    state_definitions = get_resilience_state_definitions()
  )
  
  # Integrate with existing TSnetwork states
  integrated_states <- integrate_with_existing_states(existing_states, resilience_states)
  
  # Add to TSnetwork object using existing framework
  ts_data <- add_state_classification(ts_data, integrated_states, 
                                     classification_type = "resilience")
  
  return(ts_data)
}

#' Calculate Composite Capacity Scores for TSnetwork
#' 
#' Uses TSnetwork's existing metric aggregation framework
#' @param resilience_metrics List of calculated resilience metrics
#' @return Composite scores for each capacity type
calculate_capacity_scores_tsnetwork <- function(resilience_metrics) {
  
  # Use TSnetwork's metric normalization utilities
  normalized_metrics <- normalize_metrics_tsnetwork(resilience_metrics)
  
  # Absorptive capacity score
  absorptive_score <- aggregate_metrics(
    list(
      vsi = normalized_metrics$variance_stability,
      arch = normalized_metrics$resilience_arch,
      cv = normalized_metrics$adaptive_cv
    ),
    weights = c(0.4, 0.3, 0.3),
    method = get_aggregation_method("absorptive")
  )
  
  # Restorative capacity score
  restorative_score <- aggregate_metrics(
    list(
      return_time = normalized_metrics$return_time,
      recovery_slope = normalized_metrics$recovery_slope,
      engineering_resilience = normalized_metrics$engineering_resilience
    ),
    weights = c(0.4, 0.3, 0.3),
    method = get_aggregation_method("restorative")
  )
  
  # Adaptive capacity score
  adaptive_score <- aggregate_metrics(
    list(
      sample_entropy = normalized_metrics$resilience_sample_entropy,
      dfa = normalized_metrics$resilience_dfa,
      ac_ratio = normalized_metrics$adaptive_capacity_ratio
    ),
    weights = c(0.4, 0.3, 0.3),
    method = get_aggregation_method("adaptive")
  )
  
  return(list(
    absorptive = absorptive_score,
    restorative = restorative_score,
    adaptive = adaptive_score
  ))
}

#' Define Resilience States for TSnetwork Integration
#' 
#' Creates state definitions compatible with TSnetwork's state system
#' @return List of resilience state definitions
get_resilience_state_definitions <- function() {
  
  # Follow TSnetwork's state definition structure
  state_definitions <- list(
    "resilient_stable" = list(
      description = "High absorptive and restorative capacity",
      conditions = list(
        absorptive = "> 0.8",
        restorative = "> 0.8"
      ),
      color = "green",
      priority = 1
    ),
    
    "resilient_robust" = list(
      description = "Balanced capacity across all dimensions",
      conditions = list(
        absorptive = "> 0.6",
        restorative = "> 0.6",
        adaptive = "> 0.6"
      ),
      color = "blue",
      priority = 2
    ),
    
    "resilient_recovering" = list(
      description = "Positive recovery trends",
      conditions = list(
        restorative = "> 0.5",
        trend = "improving"
      ),
      color = "orange",
      priority = 3
    ),
    
    "resilient_vulnerable" = list(
      description = "Low absorptive capacity, at-risk state",
      conditions = list(
        absorptive = "< 0.3"
      ),
      color = "red",
      priority = 4
    ),
    
    "resilient_fluctuating" = list(
      description = "Variable performance across capacities",
      conditions = list(
        default = TRUE
      ),
      color = "purple",
      priority = 5
    )
  )
  
  return(state_definitions)
}
```

### Phase 4: Integration with TSnetwork Visualization and Reporting

#### Task 4.1: Extend TSnetwork's Visualization Functions
**File: `R/resilience_visualization.R`**

```r
#' Plot Resilience Analysis Using TSnetwork's Plotting Framework
#' 
#' Extends TSnetwork's existing plotting capabilities with resilience-specific visualizations
#' @param ts_data TSnetwork object with resilience analysis
#' @param plot_type Type of resilience plot to create
#' @return Plot object compatible with TSnetwork's visualization system
#' @export
plot_resilience_tsnetwork <- function(ts_data, plot_type = "capacity_timeline") {
  
  # Use TSnetwork's existing plotting utilities and themes
  base_plot_config <- get_tsnetwork_plot_config()
  
  switch(plot_type,
    "capacity_timeline" = {
      # Use TSnetwork's timeline plotting framework
      create_capacity_timeline_plot(ts_data, base_plot_config)
    },
    
    "resilience_states" = {
      # Extend TSnetwork's state visualization
      create_resilience_state_plot(ts_data, base_plot_config)
    },
    
    "capacity_radar" = {
      # Create radar chart using TSnetwork's multi-metric visualization
      create_capacity_radar_plot(ts_data, base_plot_config)
    },
    
    "metric_correlation" = {
      # Use TSnetwork's correlation visualization utilities
      create_resilience_correlation_plot(ts_data, base_plot_config)
    }
  )
}

#' Add Resilience Metrics to TSnetwork Dashboard
#' 
#' Integrates resilience analysis into TSnetwork's existing dashboard framework
#' @param dashboard_object Existing TSnetwork dashboard
#' @param resilience_data TSnetwork object with resilience analysis
#' @return Enhanced dashboard with resilience components
#' @export
add_resilience_to_dashboard <- function(dashboard_object, resilience_data) {
  
  # Use TSnetwork's dashboard extension framework
  resilience_panel <- create_resilience_dashboard_panel(resilience_data)
  
  # Add to existing dashboard using TSnetwork's panel system
  enhanced_dashboard <- add_dashboard_panel(dashboard_object, resilience_panel,
                                          position = "resilience_analysis")
  
  return(enhanced_dashboard)
}
```

#### Task 4.2: Extend TSnetwork's Reporting System
**File: `R/resilience_reporting.R`**

```r
#' Generate Resilience Report Using TSnetwork's Reporting Framework
#' 
#' Creates resilience analysis reports integrated with TSnetwork's existing reporting system
#' @param ts_data TSnetwork object with resilience analysis
#' @param report_template Template compatible with TSnetwork's reporting system
#' @return Report object following TSnetwork's format
#' @export
generate_resilience_report_tsnetwork <- function(ts_data, report_template = "standard") {
  
  # Use TSnetwork's existing report generation framework
  base_report <- initialize_tsnetwork_report(ts_data, template = report_template)
  
  # Add resilience-specific sections
  resilience_summary <- create_resilience_summary_section(ts_data)
  capacity_analysis <- create_capacity_analysis_section(ts_data)
  state_evolution <- create_state_evolution_section(ts_data)
  
  # Integrate with existing TSnetwork report structure
  enhanced_report <- add_report_sections(base_report, list(
    "resilience_summary" = resilience_summary,
    "capacity_analysis" = capacity_analysis,
    "state_evolution" = state_evolution
  ))
  
  return(enhanced_report)
}
```

### Phase 5: Testing and Documentation Integration

#### Task 5.1: Extend TSnetwork's Testing Framework
**File: `tests/testthat/test-resilience-metrics.R`**

```r
# Follow TSnetwork's existing testing patterns and utilities

test_that("Resilience metrics integrate properly with TSnetwork objects", {
  
  # Use TSnetwork's test data creation utilities
  test_ts_data <- create_test_tsnetwork_object()
  
  # Test resilience metric calculations
  resilience_ts_data <- add_resilience_analysis(test_ts_data)
  
  # Use TSnetwork's assertion utilities
  expect_tsnetwork_object_valid(resilience_ts_data)
  expect_metrics_calculated(resilience_ts_data, "resilience")
  
  # Test state classification integration
  classified_data <- classify_resilience_states_tsnetwork(resilience_ts_data)
  expect_states_valid(classified_data, "resilience")
})

test_that("Resilience metrics work with existing TSnetwork workflows", {
  
  # Test integration with existing TSnetwork analysis pipelines
  test_data <- create_test_pipeline_data()
  
  # Run standard TSnetwork analysis
  standard_analysis <- run_tsnetwork_analysis(test_data)
  
  # Add resilience analysis
  enhanced_analysis <- add_resilience_analysis(standard_analysis)
  
  # Verify compatibility
  expect_analysis_compatible(enhanced_analysis)
  expect_no_conflicts_with_existing_metrics(enhanced_analysis)
})
```

#### Task 5.2: Integrate Documentation with TSnetwork's Documentation System
**Documentation Standards:**

1. **Follow TSnetwork's roxygen2 documentation patterns**
2. **Use consistent parameter naming with existing TSnetwork functions**
3. **Include examples that work with TSnetwork's example datasets**
4. **Add resilience analysis vignettes to TSnetwork's vignette collection**
5. **Update TSnetwork's main documentation to reference resilience capabilities**

**Example Documentation Pattern:**
```r
#' Calculate Resilience Metrics for TSnetwork Objects
#' 
#' This function extends TSnetwork's time series analysis capabilities by adding
#' comprehensive resilience metrics across three capacity dimensions.
#' 
#' @param ts_data A TSnetwork time series object created with \code{\link{create_tsnetwork}}
#' @param capacity_types Character vector specifying which capacity types to analyze.
#'   Options: "absorptive", "restorative", "adaptive". Default is all three.
#' @param window_size Numeric. Rolling window size for metric calculations. 
#'   If NULL, uses TSnetwork's default window sizing algorithm.
#' @param baseline_method Character. Method for establishing baseline performance.
#'   Compatible with TSnetwork's baseline calculation methods.
#' 
#' @return A TSnetwork object with added resilience metrics. All existing
#'   TSnetwork functionality remains available.
#' 
#' @examples
#' # Create TSnetwork object
#' ts_data <- create_tsnetwork(example_timeseries)
#' 
#' # Add comprehensive resilience analysis
#' resilient_ts <- add_resilience_analysis(ts_data)
#' 
#' # Plot resilience capacity timeline
#' plot_resilience_tsnetwork(resilient_ts, "capacity_timeline")
#' 
#' # Classify resilience states
#' classified_ts <- classify_resilience_states_tsnetwork(resilient_ts)
#' 
#' # Generate integrated report
#' report <- generate_resilience_report_tsnetwork(classified_ts)
#' 
#' @seealso \code{\link{classify_resilience_states_tsnetwork}}, 
#'   \code{\link{plot_resilience_tsnetwork}}, \code{\link{create_tsnetwork}}
#' 
#' @family resilience_analysis
#' @family tsnetwork_extensions
#' @export
```

## Implementation Priorities and Dependencies

### Priority 1: Core Integration (Week 1-2)
- Analyze TSnetwork's existing codebase thoroughly
- Identify reusable utilities and patterns
- Create resilience metric functions that leverage existing TSnetwork infrastructure
- Ensure seamless integration with TSnetwork's data objects

### Priority 2: Metric Implementation (Week 3-4)
- Implement absorptive capacity metrics using TSnetwork's rolling analysis framework
- Implement restorative capacity metrics leveraging TSnetwork's event detection
- Implement adaptive capacity metrics extending TSnetwork's complexity analysis
- Test compatibility with existing TSnetwork workflows

### Priority 3: State Classification (Week 5)
- Extend TSnetwork's state classification system with resilience-specific states
- Ensure resilience states integrate well with existing TSnetwork state definitions
- Test state transition analysis and validation

### Priority 4: Visualization and Reporting (Week 6)
- Extend TSnetwork's plotting capabilities with resilience visualizations
- Integrate resilience analysis into TSnetwork's dashboard framework
- Enhance TSnetwork's reporting system with resilience sections

### Priority 5: Documentation and Testing (Week 7)
- Follow TSnetwork's documentation standards and patterns
- Extend TSnetwork's testing framework for resilience metrics
- Create vignettes showing resilience analysis integrated with TSnetwork workflows

## Success Criteria

1. **Seamless Integration**: Resilience metrics work transparently with all existing TSnetwork functionality
2. **No Duplication**: All resilience functions leverage TSnetwork's existing utilities where possible
3. **Consistent API**: New functions follow TSnetwork's naming conventions and parameter patterns
4. **Enhanced Capability**: TSnetwork users can perform comprehensive resilience analysis without learning new interfaces
5. **Backward Compatibility**: All existing TSnetwork functionality remains unchanged and available