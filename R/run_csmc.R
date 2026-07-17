#' Controlled Sequential Monte Carlo (CSMC)
#'
#' @param lag Rolling window length (passed per your specification, though full T pass is used here)
#' @param data List containing the observations (e.g., data$obs)
#' @param model List containing model parameters, likelihoods, and dimensions
#' @param Napf Number of particles (N)
#'
#' @return A list containing the approximated log marginal likelihood, smoothed trajectories, and final twisting functions.
#' 
#' @export
run_CSMC <- function(data, model, Napf) {
  
  # Extract basic configurations
  obs <- data$obs
  Time <- nrow(obs)
  d <- ncol(obs)
  K <- model$parameters$k
  
  # --- 1: Initialization ---
  # Set psi_1 = ... = psi_{T+1} = 1
  # In log-space, this means evaluating to 0. Passing `NA` perfectly handles this 
  # based on your `evaluate_psi` and `evaluate_psi_tilde` helper logic.
  psi <- vector("list", Time + 1)
  for (t in 1:(Time + 1)) {
    psi[[t]] <- rep(NA, 2 * d)
  }
  
  # --- 2: CSMC Iterations (k = 1, ..., K) ---
  for (k in 1:K) {
    
    # Input: Initial particle system H0
    H_prev <- list(
      X = matrix(0, nrow = Napf, ncol = d),
      logW = rep(-log(Napf), Napf),
      logZ = 0,
      log_li = rep(0, Napf),
      anc = 1:Napf
    )
    
    H_history <- vector("list", Time)
    
    # The 'init' flag handles the difference between the unguided first iteration 
    # and the subsequently guided iterations.
    init_flag <- (k == 1)
    
    # --- 3: Forward Pass (t = 1, ..., T) ---
    for (t in 1:Time) {
      
      # 4: sample Ht <- psi-apf(t, psi_t, H_{t-1})
      apf_out <- run_psi_APF_rolling(data = data, 
                                     t = t, 
                                     psi_t = psi[[t]], 
                                     H_prev = H_prev, 
                                     model = model, 
                                     init = init_flag,
                                     target_N = Napf)
      
      H_prev <- apf_out$H
      H_history[[t]] <- H_prev
    }
    
    # --- 5: Backward Pass (t = T, ..., 1) ---
    for (t in Time:1) {
      
      # 6: set psi_t <- learn-psi(t, psi_{t+1}, H_t)
      # Calculates the twisting log_psi function arguments using the particles and log-likelihoods
      psi[[t]] <- learn_psi(psi_t = psi[[t + 1]], 
                            H_prev = H_history[[t]], 
                            model = model)
    }
  }
  
  # --- Output: Approximations ---
  # 1. p(y_1:T) approximate Z_T
  final_logZ <- H_history[[Time]]$logZ
  final_logW <- H_history[[Time]]$logW
  
  # 2. Extract particle lineages recursively: X^{(n)}_{1:T}
  X_trajectories <- array(0, dim = c(Napf, d, Time))
  ancestry <- 1:Napf
  
  for (t in Time:1) {
    X_trajectories[, , t] <- H_history[[t]]$X[ancestry, ]
    # Update ancestry for the previous time step
    if (t > 1) {
      ancestry <- H_history[[t]]$anc[ancestry]
    }
  }
  
  return(list(
    log_marginal_likelihood = final_logZ,
    log_weights = final_logW,
    trajectories = X_trajectories,
    psi = psi,
    H_history = H_history
  ))
}
