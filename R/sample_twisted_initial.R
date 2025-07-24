#' Function to provide samples from the twisted initial distribution
#'
#'This function samples from a twisted initial distribution specified as a Gaussian mixture.
#'
#' @param params List of Gaussian mixture initial parameters
#' @param psi Parameters of the twisting function
#' @param N Number of samples to provide
#'
#' @return Samples from the twisted initial distribution
#' @export
#'
sample_twisted_initial <- function(params, psi, N){

  output <- compute_twisted_params(params, psi)

  mu <- output[[1]]

  cov <- output[[2]]

  samples <- mvnfast::rmvn(N, mu, cov)

  return(samples)
}

#' @import mvnfast
