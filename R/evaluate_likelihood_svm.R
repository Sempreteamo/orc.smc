#' Function to evaluate the log-likelihood of the observations for svm models
#'
#'This function evaluates the log-likelihood of the observations for stochastic volatility obs_paramss.
#'It is expected to be updated to evaluate log-density of general potential functions.
#'
#' @param x State at which to evaluate
#' @param datum Data point at which to evaluate
#' @param obs_params List containing parameters information of observation density
#'
#' @return Log-density of the observation density
#' @export
#'
evaluate_likelihood_svm <- function(x, datum, obs_params) {



  obs_cov <- obs_params[[2]]
  d <- length(x)

  likelihood <- stats::dnorm(datum, mean = 0, sd = sqrt(obs_cov*exp(x)), log = TRUE)

  return(likelihood)
}
#' @import stats

