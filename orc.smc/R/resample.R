#' Function to resample
#'
#'Carry out resampling by drawing the indices according to a specified vector of log weights
#'and the specified resampling scheme.
#'
#' @param logW Log weights
#' @param mode Which resampling scheme to use. Now we support multivariate resampling and residual resampling.
#'
#' @return List of particle indices
#' @export
#'
resample <- function(logW, mode = 'res'){

  N <- length(logW)
  w_ <- normalise_weights_in_log_space(logW)[[1]]
  #w_ <- logW

  if(mode == 'multi'){

    indices <- sample(1:N, N, replace = TRUE, prob = w_)
    return(indices)

  }else{
    Ntm <- as.integer(N*w_)

    indices <- unlist(lapply(1:N, function(i) {rep(i, Ntm[i])}))
    mr <- N - sum(Ntm)

    w_hat <- w_ - Ntm/N
    w_hat <- w_hat*N/mr

    indices <- c(sample(1:N, mr, replace = TRUE, prob = w_hat), indices)
    return(indices)
  }
}
