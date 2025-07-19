#' Enhanced Regime Detection with Multiple Methods (v2 Enhanced)
#'
#' Detects regime changes in complexity data using multiple sophisticated methods
#' including cumulative peaks, changepoint detection, variance shifts, threshold analysis, gradient changes, and entropy analysis.
#' Preserves all original input columns and provides comprehensive analysis with stability classification and pattern recognition.
#'
#' @param complexity_data Data frame with complexity measurements or numeric vector.
#'   If a data frame, it should ideally have a 'time_index' column and all original columns will be preserved in the output.
#' @param method Character string specifying detection method. Options:
#'   \itemize{
#'     \item \code{"cumulative_peaks"}: Detects cumulative complexity peaks using Z-tests.
#'     \item \code{"changepoint"}: Statistical change point detection (multi-window mean-shift test).
#'     \item \code{"threshold"}: Adaptive quartile-based regime classification.
#'     \item \code{"variance_shift"}: Detects changes in variance patterns.
#'     \item \code{"gradient"}: Detects changes in gradient patterns (slope of rolling linear models).
#'     \item \code{"entropy"}: Detects changes in the Shannon entropy of the complexity series, calculated in rolling windows.
#'     \item \code{"smart"} (default): Combines gradient, peaks, and changepoint methods.
#'     \item \code{"all"}: Applies all individual methods (cumulative_peaks, changepoint, threshold, variance_shift, gradient, entropy) and uses ensemble voting.
#'   }
#' @param sensitivity Character string for detection sensitivity. Options: "low", "medium", "high".
#'   Default is "medium". This affects thresholds and window sizes within detection methods.
#' @param min_change_size Minimum number of observations between changes. If NULL, auto-detected
#'   (typically 10% of data length, min 10). Default is NULL.
#' @param complexity_col Column name for complexity values if \code{complexity_data} is a data frame.
#'   This column is assumed to contain the pre-computed complexity scores. Auto-detected if NULL. Default is NULL.
#' @param consolidate_similar Logical, whether to merge regimes that are statistically similar after initial detection.
#'   Default is FALSE.
#' @param similarity_threshold Numeric (0-1), threshold for considering regimes similar enough to merge if \code{consolidate_similar} is TRUE.
#'   Based on normalized mean difference. Default is 0.6.
#' @param window_size Integer, base window size for rolling calculations. This is further adjusted by \code{sensitivity}.
#'   Default is 10.
#' @param peak_threshold Numeric, base Z-score threshold for individual peak detection in \code{cumulative_peaks}.
#'   Adjusted by \code{sensitivity}. Default is 1.96.
#' @param cumulative_threshold Numeric (0-1), base proportion threshold for identifying cumulative peak regions.
#'   Adjusted by \code{sensitivity}. Default is 0.6.
#'
#' @return Data frame containing all original columns from \code{complexity_data} (if it was a data frame)
#'   plus new columns with regime information:
#'   \describe{
#'     \item{time_index}{Time index for each observation (ensured or created)}
#'     \item{complexity}{Original complexity values that were analyzed (if input was vector, this is it; if data.frame, this is the content of \code{complexity_col})}
#'     \item{regime_change}{Logical indicating regime change points}
#'     \item{regime_id}{Numeric regime identifier (re-calculated after any consolidation)}
#'     \item{regime_name}{Descriptive regime names (Pattern A, B, C, etc.)}
#'     \item{regime_stability}{Categorical stability: "Stable", "Transitional", "Unstable"}
#'     \item{complexity_pattern}{Pattern description: "Increasing", "Decreasing", "Stable", "Fluctuating"}
#'     \item{stability_score}{Numeric stability score (0-1)}
#'     \item{change_type}{Type of change detected by the method}
#'     \item{change_magnitude}{Magnitude of the change (method-specific interpretation)}
#'     \item{confidence}{Confidence in the detection (method-specific interpretation, typically 0-1)}
#'   }
#'
#' @details The function first processes the input and determines parameters based on the chosen \code{sensitivity}
#'   and other specific threshold arguments. The specified \code{complexity_col} is treated as the direct complexity score.
#'   It then applies the selected \code{method}.
#'   If \code{consolidate_similar} is TRUE, an attempt is made to merge regimes that are statistically
#'   similar based on their mean complexity values.
#'   Finally, stability markings and complexity patterns are calculated for each point.
#'   All original columns from the input data frame are preserved in the output.
#'
#'   Required packages: \code{stats} (for lm, median, quantile, sd, var, coef, hist). \code{zoo} is recommended for
#'   more robust rolling calculations if available (though basic fallbacks are used).
#'
#' @examples
#' # Create sample complexity data
#' set.seed(123)
#' n_obs <- 200
#' complexity_example_df <- data.frame(
#'   original_id = paste0("ID", 1:n_obs),
#'   time_index = 1:n_obs,
#'   category = sample(letters[1:3], n_obs, replace = TRUE),
#'   complexity_score_value = c(rnorm(50, 0.3, 0.05), rnorm(50, 0.8, 0.08),
#'                              rnorm(50, 0.5, 0.05), rnorm(50, 0.9, 0.1))
#' )
#'
#' # Detect regimes using the 'smart' method, preserving original columns
#' regimes_smart <- detect_regime(
#'   complexity_example_df,
#'   complexity_col = "complexity_score_value",
#'   method = "smart",
#'   sensitivity = "medium"
#' )
#' print(head(regimes_smart))
#' print(colnames(regimes_smart)) # Should include original_id, category
#'
#' # Detect regimes using the new 'entropy' method
#' regimes_entropy <- detect_regime(
#'   complexity_example_df,
#'   complexity_col = "complexity_score_value",
#'   method = "entropy",
#'   sensitivity = "medium",
#'   window_size = 15 # Entropy might benefit from a slightly larger window
#' )
#' print(head(regimes_entropy[, c("time_index", "complexity_score_value", "regime_id", "change_type")]))
#'
#' @importFrom stats median quantile sd coef lm var aggregate predict rnorm setNames na.pass
#' @importFrom zoo na.approx na.locf
#' @importFrom utils globalVariables
#' @export
detect_regime <- function(complexity_data,
                         method = "smart",
                         sensitivity = "medium",
                         min_change_size = NULL,
                         complexity_col = NULL,
                         consolidate_similar = FALSE,
                         similarity_threshold = 0.6,
                         window_size = 10, 
                         peak_threshold = 1.96, 
                         cumulative_threshold = 0.6) {
  
  # --- Store original data if it's a data frame ---
  original_df <- NULL
  if (is.data.frame(complexity_data)) {
    original_df <- complexity_data
  }
  
  # --- Input validation and processing ---
  # This extracts the complexity vector and time_index, potentially creating time_index
  processed_input <- .process_regime_input_v2(complexity_data, complexity_col)
  complexity_vector <- processed_input$complexity # This is the vector to be analyzed
  time_idx_vector <- processed_input$time_index # This is the time index vector for analysis
  
  # Determine the actual name of the complexity column used for analysis (for later reference if needed)
  # .process_regime_input_v2 returns the vector, not the col name it used if auto-detected.
  # We need to re-evaluate if complexity_col was NULL.
  analyzed_complexity_col_name <- complexity_col
  if (is.data.frame(complexity_data) && is.null(complexity_col)) {
    analyzed_complexity_col_name <- .detect_complexity_column_v2(complexity_data)
  } else if (is.numeric(complexity_data) && is.vector(complexity_data)) {
    analyzed_complexity_col_name <- "complexity_input_vector" # Placeholder name
  }
  # If complexity_col was provided, analyzed_complexity_col_name is already set.
  
  
  n <- length(complexity_vector)
  
  if (n == 0) {
    stop("Input complexity data is empty after processing.", call. = FALSE)
  }
  
  if(any(is.na(complexity_vector))) {
    warning("Complexity data contains NA values. Attempting imputation.", call. = FALSE)
    if(all(is.na(complexity_vector))) stop("All complexity values are NA. Cannot proceed.", call. = FALSE)
    
    if(requireNamespace("zoo", quietly = TRUE)){
      complexity_vector_imputed <- zoo::na.approx(complexity_vector, na.rm = FALSE)
      if(any(is.na(complexity_vector_imputed))) complexity_vector_imputed <- zoo::na.locf(complexity_vector_imputed, na.rm = FALSE)
      if(any(is.na(complexity_vector_imputed))) complexity_vector_imputed <- zoo::na.locf(complexity_vector_imputed, na.rm = FALSE, fromLast = TRUE)
      complexity_vector <- complexity_vector_imputed # Use imputed vector for analysis
    }
    if(any(is.na(complexity_vector))) { 
      mean_val <- mean(complexity_vector, na.rm = TRUE)
      if(is.na(mean_val)) mean_val <- median(complexity_vector, na.rm = TRUE) 
      if(is.na(mean_val)) mean_val <- 0 
      complexity_vector[is.na(complexity_vector)] <- mean_val # Modify the analysis vector
      warning("NA values imputed with mean/median or 0. Results might be affected.", call. = FALSE)
    }
  }
  
  if (is.null(min_change_size)) {
    min_change_size <- max(10, floor(n * 0.10)) 
    message("Auto-detected min_change_size: ", min_change_size)
  }
  min_change_size <- max(1, floor(min_change_size)) 
  
  params <- .get_sensitivity_params_v2(sensitivity, window_size, peak_threshold, cumulative_threshold)
  
  valid_methods_single <- c("cumulative_peaks", "changepoint", "threshold", "variance_shift", "gradient", "entropy")
  
  if (method == "smart") {
    detection_result <- .detect_smart_combination_v2(complexity_vector, time_idx_vector, params, min_change_size)
  } else if (method == "all") {
    detection_result <- .apply_all_methods_v2(complexity_vector, time_idx_vector, params, min_change_size, valid_methods_single)
  } else if (method %in% valid_methods_single) {
    detection_result <- .apply_single_method_v2(complexity_vector, time_idx_vector, method, params, min_change_size)
  } else {
    stop("Invalid method specified. Choose from: ",
         paste(c(valid_methods_single, "smart", "all"), collapse = ", "), call. = FALSE)
  }
  
  # --- Construct the results data frame ---
  # This df contains the core results based on the (potentially imputed) complexity_vector
  result_regime_info <- data.frame(
    time_index = time_idx_vector, # Use the time_index aligned with complexity_vector
    complexity = processed_input$complexity, # Report the original (pre-imputation) complexity values from input
    regime_change = detection_result$regime_change,
    regime_id = detection_result$regime_id,
    change_type = detection_result$change_type,
    change_magnitude = detection_result$change_magnitude,
    confidence = detection_result$confidence,
    stringsAsFactors = FALSE
  )
  
  # --- Consolidate similar regimes if requested ---
  if (consolidate_similar) {
    message("Attempting to consolidate similar regimes...")
    # Pass result_regime_info which has the current regime_id and original complexity
    consolidation_output <- .consolidate_similar_regimes_v2(result_regime_info, similarity_threshold, min_change_size)
    result_regime_info <- consolidation_output$result 
    if (nrow(result_regime_info) > 1) {
      result_regime_info$regime_change <- c(FALSE, diff(result_regime_info$regime_id) != 0)
    } else if (nrow(result_regime_info) == 1) { # Handle single row case
      result_regime_info$regime_change <- FALSE
    }
    result_regime_info$regime_change <- .apply_min_change_constraint(result_regime_info$regime_change, min_change_size)
    result_regime_info$regime_id <- .generate_regime_ids(result_regime_info$regime_change) 
    
    num_consolidated <- sum(consolidation_output$is_consolidated_event) # This was a named vector
    message(sum(num_consolidated > 0), " original regimes were involved in consolidation events.")
  }
  
  # --- Calculate stability and patterns ---
  # These will be added to result_regime_info
  stability_info <- .calculate_stability_markings_v2(result_regime_info, min_change_size)
  result_regime_info$stability_score <- stability_info$stability_score
  result_regime_info$regime_stability <- stability_info$regime_stability
  result_regime_info$complexity_pattern <- .describe_complexity_patterns_v2(result_regime_info) 
  
  unique_ids_for_naming <- sort(unique(result_regime_info$regime_id))
  if (length(unique_ids_for_naming) > 0) { 
    if (length(unique_ids_for_naming) > 26*26) { # Pattern Z9, then Pattern AA1
      name_map <- stats::setNames(paste0("Pattern", 
                                         LETTERS[((unique_ids_for_naming-1) %/% 260) %% 26 + 1], # For > 260 regimes, not perfect but extends
                                         LETTERS[((unique_ids_for_naming-1) %/% 10) %% 26 + 1], 
                                         (unique_ids_for_naming-1) %% 10 +1
      ), unique_ids_for_naming)
    } else if (length(unique_ids_for_naming) > 26) { 
      name_map <- stats::setNames(paste0("Pattern", LETTERS[((unique_ids_for_naming-1) %/% 10) +1], (unique_ids_for_naming-1) %% 10 + 1), unique_ids_for_naming)
      
    } else {
      name_map <- stats::setNames(paste("Pattern", LETTERS[unique_ids_for_naming]), unique_ids_for_naming)
    }
    result_regime_info$regime_name <- name_map[as.character(result_regime_info$regime_id)]
  } else { 
    result_regime_info$regime_name <- character(nrow(result_regime_info))
  }
  
  # --- Combine with original data frame if applicable ---
  final_output_df <- NULL
  if (!is.null(original_df)) {
    # Merge result_regime_info with original_df
    # Ensure time_index column in original_df has the same name as in result_regime_info for merging
    # .process_regime_input_v2 ensures result_regime_info$time_index is consistent.
    # If original_df doesn't have a 'time_index' column, but .process_regime_input_v2 created one,
    # we need to align them. .process_regime_input_v2 uses 1:n if 'time_index' is absent.
    
    # Check if original_df has a 'time_index' column. If not, add one based on row order.
    if (!"time_index" %in% names(original_df)) {
      if (nrow(original_df) == nrow(result_regime_info)) {
        original_df$time_index <- result_regime_info$time_index # Assume direct alignment
      } else {
        warning("Original data frame does not have 'time_index' and row count mismatch. Cannot reliably merge all original columns.", call. = FALSE)
        final_output_df <- result_regime_info # Fallback to just regime info
      }
    }
    
    if(is.null(final_output_df)){ # If merge is still possible
      # Remove columns from result_regime_info that are already in original_df to avoid _x, _y suffixes,
      # except for 'time_index' (key) and 'complexity' (which we want from result_regime_info as it's what was analyzed)
      cols_to_potentially_remove_from_results <- names(result_regime_info)[names(result_regime_info) %in% names(original_df)]
      cols_to_keep_in_results_despite_name_clash <- c("time_index", "complexity") 
      cols_to_remove <- setdiff(cols_to_potentially_remove_from_results, cols_to_keep_in_results_despite_name_clash)
      
      result_for_merge <- result_regime_info
      if(length(cols_to_remove) > 0) {
        result_for_merge <- result_regime_info[, !names(result_regime_info) %in% cols_to_remove, drop = FALSE]
      }
      
      # If 'complexity' col in result_for_merge has the same name as the original analyzed_complexity_col_name,
      # and that col is also in original_df, original_df's version will be kept by merge unless names are changed.
      # Ensure the 'complexity' column in the output is indeed processed_input$complexity.
      # The current `result_regime_info` already has `complexity = processed_input$complexity`.
      # If `original_df` has a column named `analyzed_complexity_col_name`, R's merge will create .x, .y.
      # To avoid this, if `analyzed_complexity_col_name` is different from "complexity" (the name in result_regime_info),
      # we can remove `analyzed_complexity_col_name` from `original_df` before merge.
      
      temp_original_df <- original_df
      if (!is.null(analyzed_complexity_col_name) && 
          analyzed_complexity_col_name %in% names(temp_original_df) && 
          analyzed_complexity_col_name != "complexity") { # if the original col name is not "complexity"
        temp_original_df[[analyzed_complexity_col_name]] <- NULL # Remove it to avoid clash with result_for_merge$complexity
      }
      
      final_output_df <- merge(temp_original_df, result_for_merge, by = "time_index", all.x = TRUE, sort = FALSE)
      # Ensure original order if `time_index` was just 1:n
      if (identical(time_idx_vector, 1:n)) {
        final_output_df <- final_output_df[order(final_output_df$time_index),]
      }
    }
  } else {
    # Input was a vector, so result_regime_info is the final output
    final_output_df <- result_regime_info
  }
  
  message("Detected ", length(unique(final_output_df$regime_id)), " final regimes using method: ", method)
  
  return(final_output_df)
}
