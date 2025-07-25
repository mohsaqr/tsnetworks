% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scaling_utils.R
\name{auto_detect_scaling}
\alias{auto_detect_scaling}
\title{Auto-detect Best Scaling Method}
\usage{
auto_detect_scaling(x, prefer_robust = TRUE)
}
\arguments{
\item{x}{Numeric vector to analyze}

\item{prefer_robust}{Whether to prefer robust methods for outlier-heavy data}
}
\value{
Recommended scaling method
}
\description{
Automatically selects the most appropriate scaling method based on
data characteristics such as distribution, outliers, and range.
}
\examples{
# Normal-like data
normal_data <- rnorm(100)
auto_detect_scaling(normal_data)

# Data with outliers
outlier_data <- c(rnorm(95), 100, -100)
auto_detect_scaling(outlier_data)

}
