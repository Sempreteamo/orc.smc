#' Function to simulate observations for binomial distributions
#'
#' @param x_particles State which to evaluate observations at
#' @param M Number of experiments for the binomial distribution
#'
#' @return Observations generated at the specific state
#' @export
#'
simulate_observation_bin <- function(x_particles, M) {
  
  
  if (!is.matrix(x_particles)) {
    stop("x_particles must be a matrix: rows = particles, cols = dimensions")
  }
  
  # Dimensions
  n_particles <- nrow(x_particles)
  dim_x <- ncol(x_particles)
  
  # Compute probabilities p = 1 / (1 + exp(-x))
  p <- 1 / (1 + exp(-x_particles))
  
  # Clamp probabilities to avoid extremes
  p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
  
  # Expand M if scalar
  if (length(M) == 1) {
    M <- matrix(M, n_particles, dim_x)
  } else if (length(M) == dim_x) {
    M <- matrix(rep(M, each = n_particles), nrow = n_particles)
  } else if (!all(dim(M) == dim(p))) {
    stop("M must be scalar, vector of length = dimensions, or matrix matching x_particles.")
  }
  
  # Sample y ~ Bin(M, p) element-wise
  y <- matrix(stats::rbinom(length(p), size = as.vector(M), prob = as.vector(p)), 
              nrow = n_particles, ncol = dim_x)
  
  return(y)
}
#' @import stats
