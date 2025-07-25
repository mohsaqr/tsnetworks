% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scaling_utils.R
\name{reverse_scaling}
\alias{reverse_scaling}
\title{Reverse Scaling Transformation}
\usage{
reverse_scaling(scaled_data, scaling_params)
}
\arguments{
\item{scaled_data}{Scaled data (numeric vector or scaled_data object)}

\item{scaling_params}{Scaling parameters (from apply_scaling result)}
}
\value{
Numeric vector with original scale restored
}
\description{
Reverses a scaling transformation to return data to its original scale.
This is useful for interpreting results or converting predictions back
to the original scale.
}
\examples{
# Scale and then reverse
data <- c(1, 5, 10, 15, 20, 100)
scaled_result <- apply_scaling(data, method = "standardize")
reversed_data <- reverse_scaling(scaled_result$data, scaled_result$scaling_params)
all.equal(data, reversed_data)  # Should be TRUE

}
