#' Function to simulate observations for binomial distributions
#'
#' @param x_particles State which to evaluate observations at
#' @param M Number of experiments for the binomial distribution
#'
#' @return Observations generated at the specific state
#' @export
#'
simulate_observation_bin <- function(x, params) {
  M <- params[[1]]

  # Calculate success probability
  p <- 1 / (1 + exp(-x))
  p <- pmin(pmax(p, 1e-10), 1 - 1e-10)

  if (is.matrix(x)) {
    # Return a matrix when input is a particle matrix (during filtering)
    y_vals <- stats::rbinom(length(p), size = M, prob = as.vector(p))
    return(matrix(y_vals, nrow = nrow(x), ncol = ncol(x)))
  } else {
    # Return a vector when input is a state vector (during sample_obs)
    return(stats::rbinom(length(p), size = M, prob = p))
  }
}
#' @import stats
