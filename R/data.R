#' Daily Step Count Data
#'
#' A dataset containing daily step count measurements over approximately 2.4 years.
#' This dataset is used for demonstrating time series network analysis, regime detection,
#' and dynamic complexity calculations in the tsnetworks package.
#'
#' @format A data frame with 861 rows and 4 variables:
#' \describe{
#'   \item{date}{Character. Date in YYYY-MM-DD format}
#'   \item{Steps}{Numeric. Daily step count}
#'   \item{Seq}{Integer. Sequential day number (1 to 861)}
#'   \item{ID}{Character. Identifier for the individual ("Saqr")}
#' }
#' @source Personal fitness tracking data collected from November 2020 to March 2023
#' @examples
#' data(saqrsteps)
#' head(saqrsteps)
#' plot(saqrsteps$Seq, saqrsteps$Steps, type = "l", 
#'      xlab = "Day", ylab = "Steps", main = "Daily Step Count")
"saqrsteps"

#' VAR(9) Time Series Dataset
#'
#' A multivariate time series dataset generated from a Vector Autoregressive model
#' of order 9 (VAR(9)). This dataset is used for demonstrating multivariate time
#' series network analysis and testing various distance measures in the tsnetworks package.
#'
#' @format A data frame with time series observations and multiple variables
#' @source Simulated data from a VAR(9) model for testing purposes
#' @examples
#' data(var9)
#' head(var9)
#' # Basic summary of the dataset
#' summary(var9)
"var9"
