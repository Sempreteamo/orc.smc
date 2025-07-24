#' Function to compute the parameters of a twisting Gaussian mixture distribution
#'
#'This function takes the parameters of a Gaussian mixture distribution and
#'a Gaussian twisting function and returns the parameters of the corresponding twisted distribution
#'obtained by multiplying the two together and renormalizing.

#' @param params List of Gaussian mixture parameters
#' @param psi Parameters of the twisting function
#'
#' @return List of twisted Gaussian mixture parameters
#' @export
#'
compute_twisted_params <- function(params, psi){
  d <- length(psi)/2
  dist_mu <- params$mean
  dist_cov <- params$cov
  psi_mu <- psi[1:d]
  psi_cov <- diag(psi[(d+1):(d+d)], d, d)
  
  dist_prec <- solve(dist_cov)
  psi_prec <- solve(psi_cov)
  
  combined_prec <- dist_prec + psi_prec
  combined_cov <- solve(combined_prec)
  
  # Combined mean
  combined_mean <- combined_cov %*% (dist_prec %*% dist_mu + psi_prec %*% psi_mu)
  
  return(params = list(mu = as.vector(combined_mean), cov = combined_cov))
}

