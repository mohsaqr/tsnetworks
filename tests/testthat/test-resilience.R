# COMPREHENSIVE SCIENTIFIC VALIDATION TESTS FOR RESILIENCE FUNCTIONS
# These tests ensure mathematical accuracy and scientific rigor

library(testthat)

test_that("VSI calculation is mathematically correct", {
  # Test with known variance patterns
  set.seed(42)
  stable_series <- rep(5, 50)  # Constant series - should have high VSI
  volatile_series <- c(rep(1, 25), rep(10, 25))  # Change in variance - should have lower VSI
  
  data_stable <- data.frame(ts = stable_series)
  data_volatile <- data.frame(ts = volatile_series)
  
  result_stable <- calculate_resilience_metrics(data_stable, ts_cols = "ts", window_width = 20)
  result_volatile <- calculate_resilience_metrics(data_volatile, ts_cols = "ts", window_width = 20)
  
  # VSI should be higher for stable series (accessing as column in data frame)
  mean_vsi_stable <- mean(result_stable$absorptive_vsi_ts, na.rm = TRUE)
  mean_vsi_volatile <- mean(result_volatile$absorptive_vsi_ts, na.rm = TRUE)
  
  expect_true(mean_vsi_stable > mean_vsi_volatile, 
              info = "VSI should correctly detect variance stability")
  expect_true(all(result_stable$absorptive_vsi_ts >= 0 & 
                  result_stable$absorptive_vsi_ts <= 1, na.rm = TRUE),
              info = "VSI should be bounded between 0 and 1")
})

test_that("DFA alpha calculation follows scientific methodology", {
  # Test with known scaling behaviors
  set.seed(123)
  
  # White noise should have alpha ≈ 0.5
  white_noise <- rnorm(100)
  data_wn <- data.frame(ts = white_noise)
  
  # Brownian motion should have alpha ≈ 1.5  
  brownian <- cumsum(rnorm(100))
  data_bm <- data.frame(ts = brownian)
  
  result_wn <- calculate_resilience_metrics(data_wn, ts_cols = "ts", window_width = 40)
  result_bm <- calculate_resilience_metrics(data_bm, ts_cols = "ts", window_width = 40)
  
  alpha_wn <- mean(result_wn$adaptive_dfa_ts, na.rm = TRUE)
  alpha_bm <- mean(result_bm$adaptive_dfa_ts, na.rm = TRUE)
  
  # Relaxed bounds for small sample sizes
  expect_true(alpha_wn > 0.2 && alpha_wn < 1.0, 
              info = "White noise DFA alpha should be in reasonable range")
  expect_true(alpha_bm > 0.8 && alpha_bm < 2.0, 
              info = "Brownian motion DFA alpha should be in reasonable range")
  expect_true(alpha_bm > alpha_wn, 
              info = "Brownian motion should have higher alpha than white noise")
})

test_that("ARCH effects detection is statistically valid", {
  # Test with homoscedastic vs heteroscedastic series
  set.seed(456)
  
  # Homoscedastic: constant variance
  homo_series <- rnorm(100, 0, 1)
  
  # Heteroscedastic: ARCH(1) process
  hetero_series <- numeric(100)
  hetero_series[1] <- rnorm(1)
  for (i in 2:100) {
    sigma_t <- sqrt(0.1 + 0.8 * hetero_series[i-1]^2)
    hetero_series[i] <- rnorm(1, 0, sigma_t)
  }
  
  data_homo <- data.frame(ts = homo_series)
  data_hetero <- data.frame(ts = hetero_series)
  
  result_homo <- calculate_resilience_metrics(data_homo, ts_cols = "ts", window_width = 30)
  result_hetero <- calculate_resilience_metrics(data_hetero, ts_cols = "ts", window_width = 30)
  
  arch_homo <- mean(result_homo$absorptive_arch_ts, na.rm = TRUE)
  arch_hetero <- mean(result_hetero$absorptive_arch_ts, na.rm = TRUE)
  
  expect_true(arch_hetero > arch_homo, 
              info = "ARCH effects should be higher for heteroscedastic series")
  expect_true(all(result_homo$absorptive_arch_ts >= 0 & 
                  result_homo$absorptive_arch_ts <= 1, na.rm = TRUE),
              info = "ARCH effects should be bounded between 0 and 1")
})

test_that("Sample entropy calculation is algorithmically correct", {
  # Test with series of known complexity
  set.seed(789)
  
  # Regular pattern - low entropy
  regular_series <- rep(c(1, 2, 3), length.out = 60)
  
  # Random series - higher entropy
  random_series <- rnorm(60)
  
  data_regular <- data.frame(ts = regular_series)
  data_random <- data.frame(ts = random_series)
  
  result_regular <- calculate_resilience_metrics(data_regular, ts_cols = "ts", window_width = 20)
  result_random <- calculate_resilience_metrics(data_random, ts_cols = "ts", window_width = 20)
  
  entropy_regular <- mean(result_regular$adaptive_sample_entropy_ts, na.rm = TRUE)
  entropy_random <- mean(result_random$adaptive_sample_entropy_ts, na.rm = TRUE)
  
  expect_true(entropy_random > entropy_regular, 
              info = "Random series should have higher sample entropy than regular pattern")
  expect_true(all(result_random$adaptive_sample_entropy_ts >= 0, na.rm = TRUE),
              info = "Sample entropy should be non-negative")
})

test_that("Recovery time calculation is scientifically valid", {
  # Test with AR(1) processes of known half-life
  set.seed(101)
  
  # Fast recovery: φ = 0.3 → half-life ≈ 1.4
  fast_recovery <- numeric(100)
  fast_recovery[1] <- 1
  for (i in 2:100) {
    fast_recovery[i] <- 0.3 * fast_recovery[i-1] + rnorm(1, 0, 0.1)
  }
  
  # Slow recovery: φ = 0.8 → half-life ≈ 3.1
  slow_recovery <- numeric(100)
  slow_recovery[1] <- 1
  for (i in 2:100) {
    slow_recovery[i] <- 0.8 * slow_recovery[i-1] + rnorm(1, 0, 0.1)
  }
  
  data_fast <- data.frame(ts = fast_recovery)
  data_slow <- data.frame(ts = slow_recovery)
  
  result_fast <- calculate_resilience_metrics(data_fast, ts_cols = "ts", window_width = 30)
  result_slow <- calculate_resilience_metrics(data_slow, ts_cols = "ts", window_width = 30)
  
  recovery_fast <- mean(result_fast$restorative_recovery_time_ts, na.rm = TRUE)
  recovery_slow <- mean(result_slow$restorative_recovery_time_ts, na.rm = TRUE)
  
  expect_true(recovery_slow > recovery_fast, 
              info = "Higher AR coefficient should yield longer recovery time")
  expect_true(all(result_fast$restorative_recovery_time_ts > 0, na.rm = TRUE),
              info = "Recovery time should be positive")
})

test_that("Edge cases are handled properly", {
  # Test with problematic inputs
  
  # Constant series - should now work without warnings
  constant_data <- data.frame(ts = rep(5, 50))
  result_constant <- calculate_resilience_metrics(constant_data, ts_cols = "ts", window_width = 20)
  expect_true(is.data.frame(result_constant))
  expect_true("absorptive_vsi_ts" %in% names(result_constant))
  
  # Series with NAs
  na_series <- c(rnorm(25), rep(NA, 10), rnorm(15))
  na_data <- data.frame(ts = na_series)
  result_na <- calculate_resilience_metrics(na_data, ts_cols = "ts", window_width = 20)
  expect_true(is.data.frame(result_na))
  
  # Very short series
  short_data <- data.frame(ts = rnorm(15))
  result_short <- calculate_resilience_metrics(short_data, ts_cols = "ts", window_width = 10)
  expect_true(is.data.frame(result_short))
})

test_that("Statistical properties are preserved", {
  set.seed(12345)
  test_data <- data.frame(ts = rnorm(100, 10, 2))
  result <- calculate_resilience_metrics(test_data, ts_cols = "ts", window_width = 25)
  
  # Check that all metrics are finite where not NA
  resilience_cols <- grep("^(absorptive|restorative|adaptive)_", names(result), value = TRUE)
  
  for (col in resilience_cols) {
    finite_values <- result[[col]][!is.na(result[[col]])]
    expect_true(all(is.finite(finite_values)), 
                info = paste("All", col, "values should be finite"))
  }
  
  # Check dimensionality consistency
  expect_equal(nrow(result), 100,
               info = "Output should preserve input length")
})

# Performance test to ensure efficiency
test_that("Functions are computationally efficient", {
  skip_on_cran()  # Skip on CRAN due to time constraints
  
  set.seed(999)
  large_data <- data.frame(ts = rnorm(1000))
  
  # Should complete within reasonable time
  start_time <- Sys.time()
  result <- calculate_resilience_metrics(large_data, ts_cols = "ts", window_width = 50)
  end_time <- Sys.time()
  
  computation_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  expect_true(computation_time < 30, 
              info = "Computation should complete within 30 seconds for 1000 points")
  
  expect_equal(nrow(result), 1000,
               info = "Large dataset should be processed correctly")
})
