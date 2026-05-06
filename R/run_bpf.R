run_bpf <- function(data, model, N) {
  Time <- nrow(data$obs)
  d <- length(model$ini_mu)
  obs <- data$obs
  
  A <- model$tran_mu
  B <- model$tran_cov
  
  logZ <- numeric(Time)
 
  # --- t = 1 ---
  
  X <- mvnfast::rmvn(N, model$ini_mu, model$ini_cov) 
  logW <- rep(NA, N)
  for(i in 1:N){
    logW[i] <- model$eval_likelihood(X[i,], obs[1,, drop = FALSE], model$obs_params)
  }
  
  res_norm <- normalise_weights_in_log_space(logW)
  logZ[1] <- res_norm[[2]]

  indices <- resample(logW, target_N = N, mode = 'res')
  X <- X[indices, ]
  
  # --- t = 2 to Time ---
  for (t in 2:Time) {
    
    for (i in 1:N) {
      
      X[i, ] <- mvnfast::rmvn(1, A %*% X[i, ], B)
      logW[i] <- model$eval_likelihood(X[i,], obs[t,, drop = FALSE], model$obs_params)
      
    }
    
    res_norm <- normalise_weights_in_log_space(logW)
    logZ[t] <- logZ[t-1] + res_norm[[2]]
    
    
    indices <- resample(logW, target_N = N, mode = 'res')
    X <- X[indices, ]
  }
  
  return(list(logZ = logZ[Time]))
}
