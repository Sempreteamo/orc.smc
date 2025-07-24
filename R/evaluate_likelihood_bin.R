#' Compute log-likelihood for Binomial observations with logistic-transformed latent states
#'
#'
#' @param x Hidden states
#' @param y Observations
#' @param M Number of experiments for the binomial distribution
#' @param model the model that used
#'
#' @return Log-likelihood
#' @export
evaluate_likelihood_bin <- function(x, y, M, model) {
  M <- model$obs_p

  if (length(y) != length(x)) {
    stop("y and x must be of the same length")
  }

  # Compute p_t = 1 / (1 + exp(-x))
  p <- 1 / (1 + exp(-x))

  # Avoid numerical issues at p = 0 or 1
  p <- pmin(pmax(p, 1e-10), 1 - 1e-10)

  # Compute log-likelihood
  loglik <- stats::dbinom(y, size = M, prob = p, log = TRUE)

  return(sum(loglik)) # Return total log-likelihood
}
#' @import stats
