#' Function to learn twisting log_psi function parameters
#'
#'This function takes a collection of particle locations and calculates twisting log_psi function arguments.
#'
#' @param psi_t previous psi functions at each time t
#' @param H_prev a list of particles containing necessary generated from the simultaion process
#' @param model List containing model parameters
#'
#' @return Twisting log_psi function parameters
#' @export
#'
#'
#'
learn_psi <- function(psi_t, H_prev, model){
  N <- nrow(H_prev$X)
  d <- ncol(H_prev$X)
  X_prev <- H_prev$X
  log_likelihoods <- H_prev$log_li

  log_psi <- numeric(N)
  psi_pa <- numeric(2*d)

  if(all(is.na(psi_t))){

    log_psi <- log_likelihoods

  }else{
    for(i in 1:N){
      log_psi[i] <- log_likelihoods[i] + evaluate_psi_tilde(X_prev[i,], psi_t, model)

    }
  }

  psi_pa <- optimize_psi(X_prev, log_psi)

  return(psi_pa = psi_pa)

}


