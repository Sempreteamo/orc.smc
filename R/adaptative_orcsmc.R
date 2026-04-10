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
#' @importFrom stats rnorm
#' @export
#' Adaptive Online Rolling Controlled Sequential Monte Carlo (Algorithm 5)

      # Backward pass: t down to t0
      for (s in t:t0) {
        psi_pa[s,] <- learn_psi(psi_pa[s+1,, drop = FALSE], H_tilde[[s+1]], model)
      }

      # Forward pass: t0 to t
      for (s in t0:t) {
        # [CORRECTION] When s=t0, we MUST start from the fixed history H
        prev_state <- if(s == t0) H[[s]] else H_tilde[[s]]
        output <- run_psi_APF_rolling(data, s, psi_pa[s,, drop = FALSE], prev_state, model, init = FALSE)
        H_tilde[[s+1]] <- output$H
      }

      # Convergence check
      max_dist <- 0
      for(s in t0:t) {
        dist <- calc_xi(psi_pa_old_iter[s,], psi_pa[s,])
        if(dist > max_dist) max_dist <- dist
      }
      if(max_dist <= eps_K) { consecutive_stable <- consecutive_stable + 1 } else { consecutive_stable <- 0 }
    }

    # --- Final Pass (Finalizing this window) ---
    for (s in t0:t) {
      output <- run_psi_APF_rolling(data, s, psi_pa[s,, drop = FALSE], H[[s]], model, init = FALSE)
      H[[s+1]] <- output$H
    }

    # [CORRECTION] Fill filtering estimates using weighted mean of particles
    # This is essential for your MSE calculations
    current_W <- exp(H[[t+1]]$logW - max(H[[t+1]]$logW))
    current_W <- current_W / sum(current_W)
    filtering_estimates[t, ] <- colSums(H[[t+1]]$X * as.vector(current_W))

    logZ_vec[t] <- H[[t+1]]$logZ
    ess_history[t] <- 1 / sum(current_W^2)

    # --- Lag Adaptation ---
    t_prime_0 <- t
    for(s in t0:t) {
      # Compare against the snapshot taken at the very start of step t
      dist_L <- calc_xi(psi_pa[s,], psi_pa_initial_step[s,])
      if(dist_L > eps_L) {
        t_prime_0 <- s
        break
      }
    }
    t0 <- max(t_prime_0, t - L_max + 1)
    t0_history[t] <- t0
  }

  return(list(
    logZ = logZ_vec,
    f_means = filtering_estimates,
    H_forward = H,
    ess_history = ess_history,
    t0_history = t0_history
  ))
}
#' @import stats
#' @import rnorm
#' @import stats
#' @import rnorm
