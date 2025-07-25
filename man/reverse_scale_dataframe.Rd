% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scaling_utils.R
\name{reverse_scale_dataframe}
\alias{reverse_scale_dataframe}
\title{Reverse Scaling for Data Frame}
\usage{
reverse_scale_dataframe(scaled_df, cols = NULL)
}
\arguments{
\item{scaled_df}{Scaled data frame (from scale_dataframe)}

\item{cols}{Columns to reverse (default: all scaled columns)}
}
\value{
Data frame with original scales restored
}
\description{
Reverses scaling transformations applied to a data frame.
}
\examples{
# Scale and reverse a data frame
df <- data.frame(x = 1:10, y = (1:10)^2)
scaled_df <- scale_dataframe(df, cols = c("x", "y"), method = "minmax")
reversed_df <- reverse_scale_dataframe(scaled_df)
all.equal(df, reversed_df[, c("x", "y")])  # Should be TRUE

}
