#' Function to simulate observations for svm model
#'
#' @param state State which to evaluate observations at
#' @param params Parameters of the observation distribution
#'
#' @return Observations generated at the specific state
#' @export
#'
simulate_observation_svm <- function(state, params){

  den_cov <- params[[2]]
  obs <- stats::rnorm(1, 0, sqrt(den_cov*exp(state)))
  return(obs)
}
#' @import stats
