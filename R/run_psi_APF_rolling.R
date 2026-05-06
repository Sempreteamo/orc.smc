#' Function of psi_APF
#'
#' This function performs sampling and resampling procedures of
#' psi-Auxiliary Particle Filter with $/kappa$-adaptive resampling
#'
#' @param t time
#' @param data List containing required data information on which to run the filter
#' @param H_prev previous information
#' @param psi_t Parameters of the twisting function
#' @param init Index to decide whether the algorithm is at its initialized block
#' @param model A list containing information of model
#' @param target_N Check whether the number of particles change during the iteration
#'
#' @return A list containing:
#' 
#'
#' @export
#'
run_psi_APF_rolling <- function(data, t, psi_t, H_prev, model, init, target_N) { 
  N_prev <- nrow(H_prev$X) #achieve the particles input
  d <- ncol(H_prev$X)
  X_prev <- H_prev$X
  logW_prev <- H_prev$logW
  logZ_t <- H_prev$logZ
  
  ini_mu <- model$ini_mu
  ini_cov <- model$ini_cov
  A <- model$tran_mu
  B <- model$tran_cov
  obs <- data$obs
  kappa <- model$parameters$kappa
  obs_params <- model$obs_params
  
  # Step 1: Compute v^n in log domain 
  if(t == 1){
    log_psi_tilde <- rep(evaluate_psi_tilde_ini(psi_t, model), N_prev)
  }else{
    
    log_psi_tilde <- sapply(1:N_prev, function(i) evaluate_psi_tilde(X_prev[i,], psi_t, model))
  }
  
  log_v <- logW_prev + log_psi_tilde  
  logZ_v <- log_sum_exp(log_v)
  logV_normalized <- log_v - logZ_v 
  
  
  ESS <- exp(-log_sum_exp(2 * logV_normalized))
  
  
  # if N changes, then resample target_N new index
  if (ESS < kappa * N_prev || N_prev != target_N) { 
    
    ancestors <- resample(log_v, target_N) 
    logV_for_step5 <- rep(-log(target_N), target_N) 
  } else {
    ancestors <- 1:N_prev 
    logV_for_step5 <- logV_normalized 
  }
  
  # Step 4: Sample new states 
  X_new <- matrix(NA, target_N, d) 
  log_likelihoods <- rep(NA, target_N)
  log_psi_t <- rep(NA, target_N)
  
  for(i in 1:target_N){ 
    if (init == TRUE && t == 1) {
      X_new[i, ] <- mvnfast::rmvn(1, ini_mu, ini_cov)
    } else if (init == TRUE && t != 1) {
      
      X_new[i, ] <- mvnfast::rmvn(1, ini_mu + A %*% (X_prev[ancestors[i], ] - ini_mu), B)
    } else if (init == FALSE && t == 1) {
      X_new[i, ] <- sample_twisted_initial(list(mean = ini_mu, cov = as.matrix(ini_cov)), psi_t, 1)
    } else {
      
      X_new[i, ] <- sample_twisted_transition(X_prev[ancestors[i], ], model, psi_t, 1)
    }
    
    log_likelihoods[i] <- model$eval_likelihood(X_new[i,], obs[t,, drop = FALSE], obs_params)
    log_psi_t[i] <- evaluate_psi(X_new[i,], psi_t)
  }
  
  # Step 5: Compute log w_t 
  log_w <- logV_for_step5 + log_likelihoods - log_psi_t
  
  # Step 6
  logZ_w <- log_sum_exp(log_w) 
  logW <- log_w - logZ_w 
  
  logZ_new <- logZ_t + logZ_w + logZ_v  
  
  return(list(H = list(X = X_new, logW = logW, logZ = logZ_new, 
                       log_li = log_likelihoods, anc = ancestors, 
                       N_adaptive = target_N), 
              current_ess = ESS))
}

#' @import mvnfast
