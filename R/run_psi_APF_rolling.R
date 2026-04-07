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
#'
#' @return A list containing:
#' 
#'
#' @export
#'
run_psi_APF_rolling <- function(data, t, psi_t, H_prev, model, init) {
  N <- nrow(H_prev$X)
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
  
  #model$ini_cov <- B
  
  X_new <- matrix(NA, N, d)
  
  # Step 1: Compute v^n in log domain
  if(t == 1){
    log_psi_tilde <- rep(evaluate_psi_tilde_ini(psi_t, model), N)
  }else{
    log_psi_tilde <- sapply(1:N, function(i) evaluate_psi_tilde(X_prev[i,], psi_t, model))
  }
  
 
  log_v <- logW_prev + log_psi_tilde  # log(W_t-1) + log(f_t)
  
  
  # Step 2: Compute logZ_t and normalized logV
  logZ_v <- log_sum_exp(log_v)
    
  logV <- log_v - logZ_v 
  #log(normalise_weights_in_log_space(log_v)[[1]])
    #log_v - logZ_v#log(normalise_weights_in_log_space(log_v)[[1]]) ?
  
  # Step 3: Compute ESS and resample if necessary
  ESS <- exp(-log_sum_exp(2 * logV))
  
  if (ESS < kappa * N) {
    ancestors <- resample(log_v)
    cat('re at ', t)
    logV <- rep(-log(N), N)  # Reset logV after resampling
    
  } else {
    ancestors <- 1:N  # No resampling needed
  
  }
  
  # Step 4: Sample new states using f_t^ψ
  log_likelihoods <- rep(NA, N)
  
  log_psi_t <- rep(NA, N)
  
  
  for(i in 1:N){
    if (init == TRUE && t == 1) {
     
      X_new[i, ] <- mvnfast::rmvn(1, ini_mu, ini_cov)
     
      
    } else if (init == TRUE && t != 1) {
     
      X_new[i, ] <- mvnfast::rmvn(1, ini_mu + A %*% (X_prev[ancestors[i], ] - ini_mu), B)
      
    } else if (init == FALSE && t == 1) {
    
      X_new[i, ] <- sample_twisted_initial(list(mean = ini_mu, cov = as.matrix(ini_cov)), psi_t, 1)
      
      
     
    } else {
      
      X_new[i, ] <- sample_twisted_transition(X_prev[ancestors[i], ], model, psi_t, 1)
    }
    
   
    log_likelihoods[i] <-  model$eval_likelihood(X_new[i,], obs[t,, drop = FALSE], obs_params)
    
    log_psi_t[i] <- evaluate_psi(X_new[i,], psi_t)
    
  }
  
  
  # Step 5: Compute log w_t
  
    log_w <- logV + log_likelihoods - log_psi_t
  
  
  
  # Step 6: Compute log Z_t and normalize weights
  logZ_w <- log_sum_exp(log_w) #normalise_weights_in_log_space(log_w)[[2]]
    #log_sum_exp(log_w)
  
  logW <- log_w - logZ_w #log(normalise_weights_in_log_space(log_w)[[1]]) #log_w - logZ_w
    
    # # Self-normalized log weights
  
  logZ_new <- logZ_t + logZ_w + logZ_v 
  
  # Return updated particle system
  return(list(H = list(X = X_new, logW = logW, logZ = logZ_new, log_li = log_likelihoods, anc = ancestors), current_ess = ESS))
}


#' @import mvnfast
