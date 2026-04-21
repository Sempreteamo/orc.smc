#' @export
adaptative_Orc_SMC <- function(B_max, K_min, K_max, K1, gamma, data, model, N1) {
  obs <- data$obs
  T_total <- nrow(obs)
  d <- ncol(obs)
  

  kappa <- rep(0, K_max)
  beta  <- rep(0, B_max)
  psi_pa <- matrix(1, nrow = T_total + 1, ncol = 2 * d)
  
  logZ_vec <- numeric(T_total)
  filtering_estimates <- matrix(NA, T_total, d)
  ess_history <- numeric(T_total)
  

  K_hist <- numeric(T_total); B_hist <- numeric(T_total); N_hist <- numeric(T_total)
  
  H <- vector('list', T_total + 1)
  H_tilde <- vector('list', T_total + 1)
  
  Bt <- 1
  Kt <- K1
  Nt <- N1
  
  # 初始化 H
  H[[1]] <- H_tilde[[1]] <- list(
    X = matrix(0, Nt, d), 
    logW = rep(-log(Nt), Nt), 
    logZ = 0
  )
  
  for (t in 1:T_total) {
    print(t)
    t0 <- max(1, t - Bt + 1)
    
    
    H_prev_trial <- H_tilde[[t]]
    
    
    out_init <- run_psi_APF_rolling(data, t, psi_pa[t, ], H_prev_trial, model, 
                                    init = TRUE, target_N = Nt) 
    H_tilde[[t+1]] <- out_init$H
    
    psi_prev_t <- psi_pa
    xi_k_array <- numeric(max(1, Kt))
    
    # Line 6: Policy Refinement
    if (Kt > 0) {
      for (k in 1:Kt) {
        psi_before_k <- psi_pa
        
        # Line 7-8: Backward Pass
        for (s in t:t0) {
          psi_pa[s, ] <- learn_psi(psi_pa[s+1,, drop = FALSE], H_tilde[[s+1]], model)
        }
        
        # Line 9-10: Forward Pass
        
        for (s in t0:t) {
          prev_s_tilde <- H_tilde[[s]] 
          
          out_tilde <- run_psi_APF_rolling(data, s, psi_pa[s, ], prev_s_tilde, model, 
                                           init = FALSE, target_N = Nt)
          H_tilde[[s+1]] <- out_tilde$H
        }
        
        xi_k_array[k] <- mean(abs(psi_pa[t0:t, , drop=FALSE] - psi_before_k[t0:t, , drop=FALSE]) / 
                                (abs(psi_pa[t0:t, , drop=FALSE]) + abs(psi_before_k[t0:t, , drop=FALSE]) + 1e-8))
      }
    }
    
    # Line 11-12: Final Forward Pass
    for (s in t0:t) {
      prev_s_H <- H[[s]]
      
      out_H <- run_psi_APF_rolling(data, s, psi_pa[s, ], prev_s_H, model, 
                                   init = FALSE, target_N = Nt)
      H[[s+1]] <- out_H$H
    }
    
    
    if (Kt > 0) {
      for (k in 1:Kt) kappa[k] <- 0.9 * kappa[k] + 0.1 * xi_k_array[k]
    }
    CK <- which(diff(kappa[1:K_max]) > 0) + 1
    Kt_next <- if (length(CK) > 0) min(CK[CK >= K_min & CK <= K_max]) else K_max
    
    for (b in 1:Bt) {
      idx_s <- t - b + 1
      if (idx_s < 1) next
      xi_b <- mean(abs(psi_pa[idx_s,] - psi_prev_t[idx_s,]) / 
                     (abs(psi_pa[idx_s,]) + abs(psi_prev_t[idx_s,]) + 1e-8))
      beta[b] <- 0.9 * beta[b] + 0.1 * xi_b
    }
    CB <- which(diff(beta[1:B_max]) > 0) + 1
    Bt_next <- if (length(CB) > 0) min(CB[CB <= B_max]) else B_max
    Bt_next <- min(Bt_next, t + 1)
    
    
    K_hist[t] <- Kt; B_hist[t] <- Bt; N_hist[t] <- Nt
    
    Nt_next <- floor(gamma / (Bt_next * (Kt_next + 1)))
    Nt_next <- max(Nt_next, 20)
    
    
    res_final <- H[[t+1]]
    w_final <- exp(res_final$logW - max(res_final$logW))
    w_final <- w_final / sum(w_final)
    filtering_estimates[t, ] <- colSums(as.matrix(res_final$X) * as.vector(w_final))
    logZ_vec[t] <- res_final$logZ
    ess_history[t] <- 1 / sum(w_final^2)
    
    Bt <- Bt_next; Kt <- Kt_next; Nt <- Nt_next
  }
  
  return(list(
    logZ = logZ_vec, 
    f_means = filtering_estimates, 
    H_forward = H, 
    ess_history = ess_history,
    K_history = K_hist, 
    B_history = B_hist, 
    N_history = N_hist
  ))
}
