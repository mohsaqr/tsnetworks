# Input Validation Functions
# Validation utilities for time series distance calculations

#' Validate time series input
#' @param series Time series data
#' @param series_name Name for error messages
#' @param allow_empty Whether to allow empty series
#' @return Validated series
validate_time_series <- function(series, series_name = "series", allow_empty = FALSE) {
  if (is.null(series)) {
    stop(paste(series_name, "cannot be NULL"))
  }
  
  if (!is.numeric(series) && !is.integer(series)) {
    stop(paste(series_name, "must be numeric"))
  }
  
  if (!allow_empty && length(series) == 0) {
    stop(paste(series_name, "cannot be empty"))
  }
  
  if (any(is.infinite(series))) {
    warning(paste(series_name, "contains infinite values"))
  }
  
  return(as.numeric(series))
}

#' Validate window parameters
#' @param window_size Window size
#' @param step_size Step size  
#' @param series_length Length of the series
validate_window_params <- function(window_size, step_size = NULL, series_length = NULL) {
  if (!is.null(window_size)) {
    if (!is.numeric(window_size) || length(window_size) != 1 || window_size <= 0) {
      stop("window_size must be a positive integer")
    }
    
    if (!is.null(series_length) && window_size > series_length) {
      stop("window_size cannot be larger than series length")
    }
  }
  
  if (!is.null(step_size)) {
    if (!is.numeric(step_size) || length(step_size) != 1 || step_size <= 0) {
      stop("step_size must be a positive integer")
    }
  }
  
  return(list(window_size = as.integer(window_size), 
              step_size = if(!is.null(step_size)) as.integer(step_size) else NULL))
}

#' Validate method parameters
#' @param method_name Distance method name
#' @param available_methods Vector of available methods
validate_method <- function(method_name, available_methods = NULL) {
  if (is.null(method_name)) {
    stop("method_name must be specified")
  }
  
  if (!is.character(method_name) || length(method_name) != 1) {
    stop("method_name must be a single character string")
  }
  
  if (!is.null(available_methods) && !method_name %in% available_methods) {
    stop(paste("method_name must be one of:", paste(available_methods, collapse = ", ")))
  }
  
  return(method_name)
}

#' Validate data frame input for time series
#' @param data Data frame
#' @param series_col_name Column name for time series
#' @param id_col_name Column name for identifiers
validate_dataframe_input <- function(data, series_col_name = NULL, id_col_name = NULL) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }
  
  if (nrow(data) == 0) {
    stop("data frame cannot be empty")
  }
  
  if (!is.null(series_col_name)) {
    if (!series_col_name %in% names(data)) {
      stop(paste("Column", series_col_name, "not found in data frame"))
    }
    
    if (!is.numeric(data[[series_col_name]])) {
      stop(paste("Column", series_col_name, "must be numeric"))
    }
  }
  
  if (!is.null(id_col_name)) {
    if (!id_col_name %in% names(data)) {
      stop(paste("Column", id_col_name, "not found in data frame"))
    }
  }
  
  return(TRUE)
}

#' Get Available Distance Methods
#'
#' Returns a character vector of distance methods that are available for use
#' within the package's validation routines.
#'
#' @return A character vector of available distance method names.
#' @keywords internal
get_available_methods <- function() {
  c("dist_cor", "dist_ccf", "dist_dtw", "dist_nmi", "dist_voi", "dist_mic", "dist_es", "dist_vr")
}

#' Check Required Packages
#'
#' Checks if specified packages are installed. If not, it stops execution
#' with an informative error message.
#'
#' @param pkgs Character vector of package names to check.
#' @param fn_name Character string. The name of the function calling this check,
#'   used for more informative error messages.
#' @keywords internal
check_required_packages <- function(pkgs, fn_name) {
  missing_pkgs <- character(0)
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_pkgs <- c(missing_pkgs, pkg)
    }
  }
  
  if (length(missing_pkgs) > 0) {
    stop(paste0("Function '", fn_name, "' requires the following package(s) to be installed: ",
                paste(missing_pkgs, collapse = ", "),
                ". Please install them using install.packages()."), call. = FALSE)
  }
}

#' Ensure two time series have the same length
#'
#' This internal helper function checks if two time series vectors have the same length.
#' If their lengths differ, it stops execution with an error message.
#'
#' @param series1 Numeric vector. The first time series.
#' @param series2 Numeric vector. The second time series.
#' @param fn_name Character string. The name of the calling function, used for error messages.
#' @keywords internal
.ensure_same_length <- function(series1, series2, fn_name) {
  if (length(series1) != length(series2)) {
    stop(paste0(fn_name, "Time series must have the same length."), call. = FALSE)
  }
}

#' Validate input for STNA function
#'
#' This internal helper function validates the input data for the `stna` function,
#' ensuring that the necessary columns are present and have the correct types.
#'
#' @param time_series A data frame containing the time series data.
#' @param value_column Character string. The name of the column containing the numeric time series values.
#' @param user_column Character string or NULL. The name of the column containing user or group identifiers.
#' @return A list containing the validated time series data, the extracted series vector, and user IDs (if applicable).
#' @keywords internal
.validate_input <- function(time_series, value_column, user_column) {
  if (!is.data.frame(time_series)) {
    stop("'time_series' must be a data frame.", call. = FALSE)
  }
  if (nrow(time_series) == 0) {
    stop("'time_series' data frame cannot be empty.", call. = FALSE)
  }
  
  if (!is.character(value_column) || length(value_column) != 1) {
    stop("'value_column' must be a single character string.", call. = FALSE)
  }
  if (!value_column %in% names(time_series)) {
    stop(paste0("Value column '", value_column, "' not found in time_series data frame."), call. = FALSE)
  }
  if (!is.numeric(time_series[[value_column]])) {
    stop(paste0("Value column '", value_column, "' must be numeric."), call. = FALSE)
  }
  
  series_vector <- time_series[[value_column]]
  
  user_ids <- NULL
  if (!is.null(user_column)) {
    if (!is.character(user_column) || length(user_column) != 1) {
      stop("'user_column' must be a single character string or NULL.", call. = FALSE)
    }
    if (!user_column %in% names(time_series)) {
      stop(paste0("User column '", user_column, "' not found in time_series data frame."), call. = FALSE)
    }
    user_ids <- time_series[[user_column]]
  }
  
  return(list(data = time_series, series = series_vector, users = user_ids))
}
