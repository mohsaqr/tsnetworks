% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/similarity_utils.R
\name{distance_to_similarity}
\alias{distance_to_similarity}
\title{Distance to Similarity Conversion and Normalization Utilities}
\usage{
distance_to_similarity(
  distance_matrix,
  method = "normalized_inverse",
  sigma = NULL,
  remove_self_sim = TRUE
)
}
\arguments{
\item{distance_matrix}{Square distance matrix}

\item{method}{Conversion method:
- "normalized_inverse": similarity = 1 - (normalized_distance) [DEFAULT]
- "inverse": similarity = 1 / (1 + distance) 
- "negative_exp": similarity = exp(-distance/sigma)
- "max_minus": similarity = (max_distance - distance) / max_distance
- "gaussian": similarity = exp(-distance^2 / (2*sigma^2))
- "reciprocal": similarity = 1 / distance (for distance > 0)}

\item{sigma}{Scale parameter for exponential/gaussian methods (auto-detected if NULL)}

\item{remove_self_sim}{Remove self-similarities (set diagonal to 0)}
}
\value{
Similarity matrix where higher values = more similar
}
\description{
Functions for converting distances to similarities and normalizing for network visualization
}
\details{
Convert Distance Matrix to Similarity Matrix

Transform distance values to similarity values using various methods.
Higher similarity values indicate more similar (closer) time series.
}
\examples{
# Create example distance matrix
distances <- matrix(c(0, 10, 20, 10, 0, 15, 20, 15, 0), nrow = 3)

# Convert to similarities
sim_norm_inv <- distance_to_similarity(distances, "normalized_inverse")
sim_inverse <- distance_to_similarity(distances, "inverse")
sim_exp <- distance_to_similarity(distances, "negative_exp")
sim_gaussian <- distance_to_similarity(distances, "gaussian")
}
