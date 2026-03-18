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
#' @return A list containing logZ estimates, filtering means, forward history, and ESS history.
#' @export
adaptative_Orc_SMC <- function(L_max, K_max, R, eps_K, eps_L, data, model, N) {
  obs <- data$obs
  Time <- nrow(obs)
  d <- ncol(obs)

  # --- Initialization ---
  X0 <- matrix(0, nrow = N, ncol = d)
  w0 <- matrix(log(1/N), 1, N)

  logZ_vec <- numeric(Time)
  ess_history <- numeric(Time)

  # Storage for filtering results and twisting function parameters
  filtering_estimates <- matrix(NA, Time, d)
  psi_pa <- matrix(NA, Time + 1, 2*d)

  # H stores final estimates; H_tilde stores intermediate results during refinement
  H <- vector('list', Time + 1)
  H_tilde <- vector('list', Time + 1)

  H[[1]] <- list(X = X0, logW = w0, logZ = 0)
  H_tilde[[1]] <- list(X = X0, logW = w0, logZ = 0)

  # --- Distance Metric (xi) for Algorithm 5 ---
  # Calculates the difference between twisting function parameters across iterations
  calc_xi <- function(psi_old, psi_new) {
    if(any(is.na(psi_old)) || any(is.na(psi_new))) return(Inf)
    max(abs(psi_old - psi_new), na.rm = TRUE)
  }

  # Initialize window start (t0) at the first observation
  t0 <- 1

  # --- Sequential Loop ---
  for (t in 1:Time) {
    # Log progress and current adaptive window status
    print(paste("Time step:", t, " | Current lag start (t0):", t0))

    # Initial pass at current time t with identity twisting function (psi = 1)
    output <- run_psi_APF_rolling(data, t, psi_pa[t, ], H_tilde[[t]], model, init = TRUE)
    H_tilde[[t+1]] <- output$H

    # Capture initial state of psi for lag adaptation (psi^(t,0) in Algo 5)
    psi_pa_initial <- psi_pa

    # --- Policy Refinement (Adaptive Iterations) ---
    Kt <- 0
    consecutive_stable <- 0

    while(Kt < K_max && consecutive_stable < R) {
      Kt <- Kt + 1

      # Snapshot of psi before this iteration to check for convergence
      psi_pa_old <- psi_pa

      # Backward pass: Learn twisting functions from t down to t0 (Algo 5: Lines 9-10)
      for (s in t:t0) {
        psi_pa[s,] <- learn_psi(psi_pa[s+1,, drop = FALSE], H_tilde[[s+1]], model)
      }

      # Forward pass: Re-simulate particles using updated psi (Algo 5: Lines 11-12)
      for (s in t0:t) {
        output <- run_psi_APF_rolling(data, s, psi_pa[s,, drop = FALSE], H_tilde[[s]], model, init = FALSE)
        H_tilde[[s+1]] <- output$H
      }

      # Convergence check: Monitor maximum change in psi within the current block [t0, t]
      max_dist <- 0
      for(s in t0:t) {
        dist <- calc_xi(psi_pa_old[s,], psi_pa[s,])
        if(dist > max_dist) max_dist <- dist
      }

      # If change is below eps_K, increment stability counter; otherwise, reset
      if(max_dist <= eps_K) {
        consecutive_stable <- consecutive_stable + 1
      } else {
        consecutive_stable <- 0
      }
    }

    # --- Final Pass ---
    # Generate final filtering distributions and evidence (logZ) using optimized psi
    for (s in t0:t) {
      output <- run_psi_APF_rolling(data, s, psi_pa[s,, drop = FALSE], H[[s]], model, init = FALSE)
      H[[s+1]] <- output$H

      if(s == t){
        ess_history[t] <- output$current_ess
      }
    }

    logZ_vec[t] <- H[[t + 1]]$logZ

    # --- Lag Adaptation (Algorithm 5: Lines 15-16) ---
    # Determine the earliest time step s (t_prime_0) where psi changed significantly
    t_prime_0 <- t
    for(s in t0:t) {
      # Compare latest psi parameters against the snapshot before this time step's iterations
      dist_L <- calc_xi(psi_pa[s,], psi_pa_initial[s,])
      if(dist_L > eps_L) {
        t_prime_0 <- s
        break # Smallest s found; exit search
      }
    }

    # Update t0 for the next time step, constrained by L_max
    t0 <- max(t_prime_0, t - L_max + 1)
  }

  return(list(logZ = logZ_vec, f_means = filtering_estimates, H_forward = H, ess_history = ess_history))
}
