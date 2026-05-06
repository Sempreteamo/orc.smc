#' Function to resample
#'
#'Carry out resampling by drawing the indices according to a specified vector of log weights
#'and the specified resampling scheme.
#'
#' @param logW Log weights
#' @param target_N Check whether the number of particles change during the iteration
#' @param mode Which resampling scheme to use. Now we support multivariate resampling and residual resampling.
#'
#' @return List of particle indices
#' @export
#'
resample <- function(logW, target_N = NULL, mode = 'res'){
  
  N_in <- length(logW)
  
  if (is.null(target_N)) target_N <- N_in 
  
  w_ <- normalise_weights_in_log_space(logW)[[1]]
  
  if(mode == 'multi'){
    # change the number of particles as target_N ---
    indices <- sample(1:N_in, size = target_N, replace = TRUE, prob = w_)
    return(indices)
    
  } else {
    
    Ntm <- as.integer(target_N * w_)
    
    indices <- unlist(lapply(1:N_in, function(i) {
      if(Ntm[i] > 0) rep(i, Ntm[i]) else NULL
    }))
    
    #
    mr <- target_N - length(indices)
    
    if (mr > 0) {
      
      w_hat <- (w_ * target_N - Ntm) / mr
      
      
      w_hat[w_hat < 0] <- 0 
      
      indices <- c(sample(1:N_in, size = mr, replace = TRUE, prob = w_hat), indices)
    }
    
    
    return(indices)
  }
}
