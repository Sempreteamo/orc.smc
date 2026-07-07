#' Controlled Sequential Monte Carlo (CSMC) - Full implementation
#'
#' @param data 
#' @param Napf 
#' @param K
#' @param model 
#' 
#' @return 
#' @export
run_CSMC <- function(data, Napf, K, model) {
  Time <- nrow(data$obs)
  d    <- length(model$ini_mu)
  
  
  logZ_history <- numeric(K)
  X0 <- matrix(0, nrow = Napf, ncol = d)
  w0 <- matrix(log(1/Napf), 1, Napf)
  logZ_vec <- numeric(Time)
  
  
  psi_pa <- matrix(NA, nrow = Time + 1, ncol = 2 * d)
  
  
  H <- vector('list', Time + 1)
  
  H[[1]] <- list(X = X0, logW = w0, logZ = 0)
  
  # 2. Main Loop: for k = 1, ..., K
  for (k in 1:K) {
    message("Iteration: ", k)
    
    H_prev <- NULL 
    for (t in 1:Time) {
      
      if(k == 1){
        res <- run_psi_APF_rolling(data = data, 
                                   t = t, 
                                   psi_t = psi_pa[t,, drop = FALSE ], 
                                   H[[t]], 
                                   model = model, 
                                   init = TRUE, 
                                   target_N = Napf)
        H[[t+1]] <- res$H
       
      }else{
        res <- run_psi_APF_rolling(data = data, 
                                   t = t, 
                                   psi_t = psi_pa[t, ], 
                                   H[[t]], 
                                   model = model, 
                                   init = FALSE, 
                                   target_N = Napf)
        H[[t+1]] <- res$H
       
      }
      
    }
    
    
    logZ_history[k] <- H[[Time+1]]$logZ
    
    
    for (t in Time:1) {
      
      psi_pa[t, ] <- learn_psi(psi_t = psi_pa[t+1, , drop = FALSE], 
                               H_prev = H[[t+1]], 
                               model = model)
    }

  }
  
  
  return(list(
    logZ_final = logZ_history[K],
    logZ_history = logZ_history,
    psi_final = psi_pa,
    H_final = H[[Time]]
  ))
}
