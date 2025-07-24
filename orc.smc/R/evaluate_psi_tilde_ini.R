#' Function to evaluate a twisting function psi_tilde initial
#'
#'This function evaluates the one-step-ahead expectation of a twisting function with
#'the specified parameters at one or more points.
#'
#' @param psi_pa Parameters of twisting function
#' @param model List containing model parameters
#'
#' @return Twisting function psi-tilde evaluated at specified points
#' @export
#'
evaluate_psi_tilde_ini <- function(psi_pa, model){



  ini <- model$ini_mu
  ini_cov <- model$ini_cov
  d <- length(ini)

  if(all(is.na(psi_pa))){

    psi_tilde <- 0

  }else{
    dif <- as.vector(ini - psi_pa[1:d])

    full_covariance <- diag(psi_pa[(d + 1):(d + d)], nrow=d, ncol=d) + ini_cov

    psi_tilde <- (-d / 2) * log(2 * pi) - (1 / 2) * log(det(full_covariance)) -
      (1 / 2) * t(dif) %*% solve(full_covariance)%*% dif

  }

  return(psi_tilde)
}
