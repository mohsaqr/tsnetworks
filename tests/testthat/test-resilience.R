# Test file for resilience analysis functions
# Tests all new resilience functionality for compatibility and correctness

test_that("calculate_resilience_metrics works with basic input", {
  # Create test data
  test_data <- data.frame(
    time = 1:50,
    ts1 = cumsum(rnorm(50)),
    ts2 = sin(seq(0, 2*pi, length.out = 50)) + rnorm(50, 0, 0.1)
  )
  
  # Test basic functionality
  result <- calculate_resilience_metrics(
    data = test_data,
    ts_cols = "ts1",
    window_width = 7
  )
  
  # Check that result is a data frame
  expect_s3_class(result, "data.frame")
  expect_s3_class(result, "resilience_data")
  
  # Check that original data is preserved
  expect_equal(nrow(result), nrow(test_data))
  expect_true("ts1" %in% names(result))
  
  # Check that resilience metrics are added
  resilience_cols <- grep("^(absorptive|restorative|adaptive)_", names(result), value = TRUE)
  expect_true(length(resilience_cols) > 0)
  
  # Check specific metric columns
  expect_true("absorptive_vsi_ts1" %in% names(result))
  expect_true("restorative_recovery_time_ts1" %in% names(result))
  expect_true("adaptive_sample_entropy_ts1" %in% names(result))
})

test_that("calculate_resilience_metrics works with multiple time series", {
  # Create test data with multiple time series
  test_data <- data.frame(
    time = 1:30,
    ts1 = cumsum(rnorm(30)),
    ts2 = cumsum(rnorm(30, sd = 1.5))
  )
  
  result <- calculate_resilience_metrics(
    data = test_data,
    ts_cols = c("ts1", "ts2"),
    window_width = 5
  )
  
  # Check that metrics exist for both time series
  expect_true("absorptive_vsi_ts1" %in% names(result))
  expect_true("absorptive_vsi_ts2" %in% names(result))
  expect_true("restorative_recovery_time_ts1" %in% names(result))
  expect_true("restorative_recovery_time_ts2" %in% names(result))
})

test_that("individual capacity functions work correctly", {
  test_data <- data.frame(
    time = 1:40,
    ts1 = cumsum(rnorm(40))
  )
  
  # Test absorptive capacity
  absorptive_result <- calculate_absorptive_capacity(
    data = test_data,
    ts_cols = "ts1",
    window_width = 7
  )
  expect_s3_class(absorptive_result, "data.frame")
  expect_true("absorptive_vsi_ts1" %in% names(absorptive_result))
  expect_true("absorptive_arch_ts1" %in% names(absorptive_result))
  expect_true("absorptive_cv_ts1" %in% names(absorptive_result))
  
  # Test restorative capacity
  restorative_result <- calculate_restorative_capacity(
    data = test_data,
    ts_cols = "ts1",
    window_width = 7
  )
  expect_s3_class(restorative_result, "data.frame")
  expect_true("restorative_recovery_time_ts1" %in% names(restorative_result))
  expect_true("restorative_eng_resilience_ts1" %in% names(restorative_result))
  
  # Test adaptive capacity
  adaptive_result <- calculate_adaptive_capacity(
    data = test_data,
    ts_cols = "ts1",
    window_width = 7
  )
  expect_s3_class(adaptive_result, "data.frame")
  expect_true("adaptive_sample_entropy_ts1" %in% names(adaptive_result))
  expect_true("adaptive_dfa_ts1" %in% names(adaptive_result))
})

test_that("classify_resilience_states works with different methods", {
  # Create test data with resilience metrics
  test_data <- data.frame(
    time = 1:30,
    ts1 = cumsum(rnorm(30))
  )
  
  resilience_data <- calculate_resilience_metrics(
    data = test_data,
    ts_cols = "ts1",
    window_width = 5
  )
  
  # Test threshold method
  threshold_result <- classify_resilience_states(
    data = resilience_data,
    method = "threshold"
  )
  expect_s3_class(threshold_result, "resilience_states")
  expect_true("resilience_state" %in% names(threshold_result))
  expect_true("resilience_state_absorptive_score" %in% names(threshold_result))
  
  # Test quantile method
  quantile_result <- classify_resilience_states(
    data = resilience_data,
    method = "quantile",
    n_states = 3
  )
  expect_s3_class(quantile_result, "resilience_states")
  expect_true("resilience_state" %in% names(quantile_result))
  
  # Check that states are valid
  states <- unique(threshold_result$resilience_state)
  states <- states[!is.na(states)]
  valid_states <- c("resilient_stable", "resilient_robust", "resilient_recovering", 
                   "resilient_vulnerable", "resilient_critical")
  expect_true(all(states %in% valid_states))
})

test_that("create_resilience_transition_matrix works correctly", {
  # Create test data with states
  test_data <- data.frame(
    time = 1:25,
    ts1 = cumsum(rnorm(25))
  )
  
  resilience_data <- calculate_resilience_metrics(test_data, ts_cols = "ts1", window_width = 5)
  state_data <- classify_resilience_states(resilience_data)
  
  # Create transition matrix
  transition_matrix <- create_resilience_transition_matrix(
    data = state_data,
    state_col = "resilience_state"
  )
  
  expect_true(is.matrix(transition_matrix))
  expect_true(all(rowSums(transition_matrix, na.rm = TRUE) <= 1.01))  # Allow for small rounding errors
  expect_equal(rownames(transition_matrix), colnames(transition_matrix))
})

test_that("get_resilience_state_summary provides comprehensive statistics", {
  # Create test data
  test_data <- data.frame(
    time = 1:30,
    ts1 = cumsum(rnorm(30))
  )
  
  resilience_data <- calculate_resilience_metrics(test_data, ts_cols = "ts1", window_width = 5)
  state_data <- classify_resilience_states(resilience_data)
  
  # Get summary
  summary_result <- get_resilience_state_summary(
    data = state_data,
    state_col = "resilience_state",
    ts_cols = "ts1"
  )
  
  expect_s3_class(summary_result, "resilience_state_summary")
  expect_true("n_states" %in% names(summary_result))
  expect_true("state_statistics" %in% names(summary_result))
  expect_true("transition_matrix" %in% names(summary_result))
  expect_true(is.numeric(summary_result$n_states))
  expect_true(summary_result$n_states > 0)
})

test_that("integration functions work with existing TSnetworks results", {
  skip_if_not_installed("tsnetworks")
  
  # Create test data
  test_data <- data.frame(
    time = 1:30,
    Steps = cumsum(rnorm(30, mean = 1000, sd = 100))
  )
  
  # Test with data frame input
  enhanced_result <- add_resilience_to_tsnetworks(
    tsnetworks_result = test_data,
    window_width = 5
  )
  
  expect_s3_class(enhanced_result, "enhanced_tsnetworks_data")
  expect_true("resilience_state" %in% names(enhanced_result))
  
  # Test extraction of resilience results
  resilience_metrics <- extract_resilience_results(enhanced_result, "metrics")
  expect_true(is.data.frame(resilience_metrics))
  expect_true(ncol(resilience_metrics) > 0)
  
  resilience_states <- extract_resilience_results(enhanced_result, "states")
  expect_true(is.data.frame(resilience_states))
  expect_true("resilience_state" %in% names(resilience_states))
})

test_that("resilience_analysis_pipeline works correctly", {
  # Create test data
  test_data <- data.frame(
    time = 1:25,
    ts1 = cumsum(rnorm(25))
  )
  
  # Run pipeline with resilience only
  pipeline_result <- resilience_analysis_pipeline(
    data = test_data,
    ts_cols = "ts1",
    analysis_types = "resilience",
    resilience_config = list(window_width = 5)
  )
  
  expect_s3_class(pipeline_result, "resilience_pipeline")
  expect_true("resilience_data" %in% names(pipeline_result))
  expect_true("resilience_summary" %in% names(pipeline_result))
  expect_equal(pipeline_result$ts_cols, "ts1")
  expect_equal(pipeline_result$analysis_types, "resilience")
})

test_that("format_for_tsnetworks maintains compatibility", {
  # Create test data with resilience analysis
  test_data <- data.frame(
    time = 1:20,
    ts1 = cumsum(rnorm(20))
  )
  
  resilience_data <- calculate_resilience_metrics(test_data, ts_cols = "ts1", window_width = 5)
  state_data <- classify_resilience_states(resilience_data)
  
  # Format for TSnetworks compatibility
  formatted_data <- format_for_tsnetworks(
    resilience_data = state_data,
    target_format = "generic"
  )
  
  expect_true(attr(formatted_data, "tsnetworks_compatible"))
  expect_equal(attr(formatted_data, "formatted_for"), "generic")
  expect_true("resilience_state" %in% names(formatted_data))
})

test_that("error handling works correctly", {
  # Test with invalid input
  expect_error(
    calculate_resilience_metrics(data = "not_a_dataframe"),
    "'data' must be a data frame"
  )
  
  # Test with empty data
  expect_error(
    calculate_resilience_metrics(data = data.frame()),
    "'data' cannot be empty"
  )
  
  # Test with invalid window width
  test_data <- data.frame(time = 1:10, ts1 = rnorm(10))
  expect_error(
    calculate_resilience_metrics(data = test_data, ts_cols = "ts1", window_width = 1),
    "'window_width' must be an integer >= 2"
  )
  
  # Test state classification without resilience metrics
  expect_error(
    classify_resilience_states(data = test_data),
    "No resilience metrics found in data"
  )
})

test_that("print methods work correctly", {
  # Test resilience_data print method
  test_data <- data.frame(time = 1:15, ts1 = cumsum(rnorm(15)))
  resilience_data <- calculate_resilience_metrics(test_data, ts_cols = "ts1", window_width = 5)
  
  expect_output(print(resilience_data), "Resilience Analysis Data")
  expect_output(print(resilience_data), "Capacity types analyzed")
  
  # Test resilience_states print method
  state_data <- classify_resilience_states(resilience_data)
  expect_output(print(state_data), "Resilience States Data")
  expect_output(print(state_data), "State distribution")
  
  # Test resilience_state_summary print method
  summary_data <- get_resilience_state_summary(state_data, ts_cols = "ts1")
  expect_output(print(summary_data), "Resilience State Summary")
  expect_output(print(summary_data), "Number of states")
  
  # Test resilience_pipeline print method
  pipeline_result <- resilience_analysis_pipeline(
    data = test_data,
    ts_cols = "ts1",
    analysis_types = "resilience"
  )
  expect_output(print(pipeline_result), "Resilience Analysis Pipeline Result")
  expect_output(print(pipeline_result), "Analysis types")
})

test_that("scaling integration works correctly", {
  # Test with different scaling methods
  test_data <- data.frame(
    time = 1:20,
    ts1 = cumsum(rnorm(20)) * 1000  # Large scale data
  )
  
  # Test auto scaling
  result_auto <- calculate_resilience_metrics(
    data = test_data,
    ts_cols = "ts1",
    scaling_method = "auto",
    window_width = 5
  )
  expect_s3_class(result_auto, "resilience_data")
  
  # Test standardize scaling
  result_std <- calculate_resilience_metrics(
    data = test_data,
    ts_cols = "ts1",
    scaling_method = "standardize",
    window_width = 5
  )
  expect_s3_class(result_std, "resilience_data")
  
  # Test no scaling
  result_none <- calculate_resilience_metrics(
    data = test_data,
    ts_cols = "ts1",
    scaling_method = "none",
    window_width = 5
  )
  expect_s3_class(result_none, "resilience_data")
})

test_that("grouped analysis works correctly", {
  # Test with grouped data
  test_data <- data.frame(
    time = rep(1:15, 2),
    group = rep(c("A", "B"), each = 15),
    ts1 = c(cumsum(rnorm(15)), cumsum(rnorm(15, sd = 1.5)))
  )
  
  # Test grouped resilience analysis
  result <- calculate_resilience_metrics(
    data = test_data,
    ts_cols = "ts1",
    id_col = "group",
    window_width = 5
  )
  
  expect_s3_class(result, "resilience_data")
  expect_equal(nrow(result), nrow(test_data))
  
  # Test grouped state classification
  state_result <- classify_resilience_states(result)
  expect_s3_class(state_result, "resilience_states")
  
  # Test grouped transition matrix
  transition_matrix <- create_resilience_transition_matrix(
    data = state_result,
    state_col = "resilience_state",
    id_col = "group"
  )
  expect_true(is.matrix(transition_matrix))
})
