#' Compute log-likelihood for Binomial observations with logistic-transformed latent states
#'
#'
#' @param x Hidden states
#' @param datum Observations
#' @param M Number of experiments for the binomial distribution
#' @param model the model that used
#'
#' @return Log-likelihood
#' @export
evaluate_likelihood_bin <- function(x, datum, obs_params) {
  # M is the number of neurons
  M <- obs_params[[1]]

  # Logistic link function mapping state to probability
  p <- 1 / (1 + exp(-x))
  p <- pmin(pmax(p, 1e-10), 1 - 1e-10) # Prevent log(0)

  # Return the sum of log-densities
  loglik <- stats::dbinom(datum, size = M, prob = p, log = TRUE)
  return(sum(loglik))
}
#' @import stats
