#' Function to return the Gaussian parameters of psi into regression coefficients
#' 
#' @param psi Vector. Length of 2d with first d's are mean and following d's are diagnasis of covariance matrix 
#' @return A list of coefficient vectors of length 2d
#' 
#' @export
#'
psi_to_coeff <- function(psi) {
  # NA means initial time t0
  # In the coefficient space, log(psi) = 0 means at initial time t0, psi = 1
  if (any(is.na(psi))) {
    return(rep(0, length(psi)))
  }
  
  d <- length(psi) / 2
  mu <- psi[1:d]
  sigma2 <- psi[(d + 1):(2 * d)]
  
  # -(x-mu)^2 / (2*sigma2) = -(1/2sigma2)x^2 + (mu/sigma2)x + const
  coeff_linear <- mu / sigma2             # coefficient of x 
  coeff_quadratic <- -1 / (2 * sigma2)   # coefficient of x^2 
  
  return(c(coeff_linear, coeff_quadratic))
}
