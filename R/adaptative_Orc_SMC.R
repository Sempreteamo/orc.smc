#' Adaptive Online Rolling Controlled Sequential Monte Carlo (Algorithm 5)
#'
#' @param L_max Maximum lag length (sliding window capacity)
#' @param K_max Maximum number of iterations for policy refinement per time step
#' @param R Patience parameter (required consecutive stable iterations for early stopping)
#' @param eps_K Convergence threshold for twisting function parameters (xi in Algo 5)
#' @param eps_L Threshold for lag adaptation (detecting significant information gain)
#' @param data List containing observation matrix 'obs' (Time x d)
#' @param model List containing model parameters and density functions
#' @param N Number of particles
#'
#' @return A list containing logZ estimates, filtering means, ESS history, and lag evolution.
#' @export
adaptative_Orc_SMC <- function(L_max, K_max, R, eps_K, eps_L, data, model, N) {
  obs <- data$obs
  Time <- nrow(obs)
  d <- ncol(obs)
  
  # --- Initialization of storage ---
  logZ_vec <- numeric(Time)
  ess_history <- numeric(Time)
  t0_history <- numeric(Time) 
  filtering_estimates <- matrix(0, Time, d)
  
  # Storage for twisting function parameters (psi)
  # Assuming 2*d parameters (e.g., mean and variance for Gaussian twists)
  psi_pa <- matrix(NA, Time + 1, 2*d) 
  
  # H stores the final confirmed particle systems
  # H_tilde stores intermediate particle systems during refinement iterations
  H <- vector('list', Time + 1)
  H_tilde <- vector('list', Time + 1)
  
  # Initial state at t=0
  X0 <- matrix(rnorm(N * d), nrow = N, ncol = d) # Or sampled from model$init
  w0 <- rep(log(1/N), N)
  H[[1]] <- list(X = X0, logW = w0, logZ = 0)
  H_tilde[[1]] <- list(X = X0, logW = w0, logZ = 0)
  
  # --- Distance Metric (xi) for Algorithm 5 ---
  # Computes the difference between twisting function parameters
  calc_xi <- function(psi_old, psi_new) {
    if(any(is.na(psi_old)) || any(is.na(psi_new))) return(Inf)
    return(max(abs(psi_old - psi_new)))
  }
  
  # Initial window start
  t0 <- 1
  
  # --- Main Sequential Loop (t = 1 to T) ---
  for (t in 1:Time) {
    cat(sprintf("\rTime Step: %d/%d | Window Start t0: %d", t, Time, t0))
    
    # [CRITICAL] Snapshot psi before starting refinement for this time step
    # This corresponds to psi^(t,0) in Algorithm 5 Line 15
    psi_pa_at_start_of_step <- psi_pa 
    
    # --- 1. Initial Pass (Algorithm 5 Lines 6-7) ---
    # Run a forward pass for the new observation using existing psi
    output_init <- run_psi_APF_rolling(data, t, psi_pa[t, ], H[[t]], model, init = TRUE)
    H_tilde[[t+1]] <- output_init$H
    
    # --- 2. Policy Refinement (Algorithm 5 Line 8) ---
    Kt <- 0
    consecutive_stable <- 0
    
    while(Kt < K_max && consecutive_stable < R) {
      Kt <- Kt + 1
      psi_pa_before_iter <- psi_pa # Used for convergence check (eps_K)
      
      # (a) Backward pass: Learn twisting functions (Algorithm 5 Lines 9-10)
      # Sweep from current time t back to window start t0
      for (s in t:t0) {
        psi_pa[s, ] <- learn_psi(psi_pa[s+1, , drop = FALSE], H_tilde[[s+1]], model)
      }
      
      # (b) Forward pass: Re-simulate particles using updated psi (Algorithm 5 Lines 11-12)
      for (s in t0:t) {
        # If s is the start of window, use fixed history H; otherwise use updated H_tilde
        prev_H <- if(s == t0) H[s] else H_tilde[s]
        output_iter <- run_psi_APF_rolling(data, s, psi_pa[s, ], prev_H, model, init = FALSE)
        H_tilde[[s+1]] <- output_iter$H
      }
      
      # (c) Convergence check (Algorithm 5 Line 13)
      # Check the maximum change in psi within the block [t0, t]
      max_dist_K <- 0
      for(s in t0:t) {
        dist <- calc_xi(psi_pa_before_iter[s, ], psi_pa[s, ])
        if(dist > max_dist_K) max_dist_K <- dist
      }
      
      if(max_dist_K <= eps_K) {
        consecutive_stable <- consecutive_stable + 1
      } else {
        consecutive_stable <- 0
      }
    }
    
    # --- 3. Finalization and Results Storage ---
    # After Kt iterations, H_tilde[[t+1]] contains the optimized distribution
    H[[t+1]] <- H_tilde[[t+1]]
    
    # Calculate filtering estimates (Weighted Mean)
    logW <- H[[t+1]]$logW
    W <- exp(logW - max(logW))
    W <- W / sum(W)
    filtering_estimates[t, ] <- colSums(H[[t+1]]$X * as.vector(W))
    
    # Record Evidence (Log-Likelihood) and ESS
    logZ_vec[t] <- H[[t+1]]$logZ
    ess_history[t] <- 1 / sum(W^2)
    
    # --- 4. Lag Adaptation (Algorithm 5 Lines 15-16) ---
    # Find the earliest time s where psi changed significantly during the Kt iterations
    t_prime_0 <- t
    for(s in t0:t) {
      dist_L <- calc_xi(psi_pa[s, ], psi_pa_at_start_of_step[s, ])
      if(dist_L > eps_L) {
        t_prime_0 <- s
        break # Smallest s found; stop searching
      }
    }
    
    # Update t0 for the next time step, constrained by L_max
    t0 <- max(t_prime_0, t - L_max + 1)
    t0_history[t] <- t0
  }
  
  cat("\nComputation Finished.\n")
  return(list(
    logZ = logZ_vec, 
    f_means = filtering_estimates, 
    ess_history = ess_history,
    t0_history = t0_history,
    final_H = H
  ))
}
