#' Function to run the Adaptive Online Rolling Controlled SMC (ORCSMC) algorithm
#'
#' It automatically tunes the number of iterations (K), block size/lag length (B), 
#' and the number of particles (N) under a fixed computational complexity budget (gamma).
#' 
#' @param B_max Integer. Maximum permissible block size (Bt).
#' @param K_min Integer. Minimum number of iterations (Kt).
#' @param K_max Integer. Maximum number of iterations (Kt).
#' @param K1 Integer. Initial number of iterations at time t=1.
#' @param gamma Numeric. Approximate computational complexity budget (gamma).
#' @param data List. Contains the observation data (e.g., data$obs).
#' @param model List. Contains model parameters, state transition, and likelihood functions.
#' @param N1 Integer. Initial number of particles at time t=1.
#'
#' @return A list containing the log-normalising constant estimates, filtering means, 
#'         forward particle systems, ESS history, and the adaptation histories for K, B, and N.
#'
#' @export
adaptative_Orc_SMC <- function(B_max, K_min, K_max, K1, gamma, data, model, N1) {
  obs <- data$obs
  T_total <- nrow(obs)
  d <- ncol(obs)
  
  kappa <- seq(1, 0, length.out = K_max) 
  beta  <- seq(1, 0, length.out = B_max) 
  
  #initialize the parameters of the psi when psi === 1
  psi_pa <- matrix(0, nrow = T_total + 1, ncol = 2 * d) 
  
  logZ_vec <- numeric(T_total)
  filtering_estimates <- matrix(NA, T_total, d)
  ess_history <- numeric(T_total)
  
  K_hist <- numeric(T_total); B_hist <- numeric(T_total); N_hist <- numeric(T_total)
  
  H <- vector('list', T_total + 1)
  H_tilde <- vector('list', T_total + 1)
  
  Bt <- 1
  Kt <- K1
  Nt <- N1
  
  H[[1]] <- H_tilde[[1]] <- list(
    X = matrix(0, Nt, d), 
    logW = rep(-log(Nt), Nt), 
    logZ = 0
  )
  
  for (t in 1:T_total) {
    if (t %% 10 == 0) cat("Time step:", t, "| B:", Bt, "| K:", Kt, "| N:", Nt, "\n") # 方便监控
    
    t0 <- max(1, t - Bt + 1)
    
    H_prev_trial <- H_tilde[[t]]
    
    #sample \tilde{H}_t <- psi-apf(t, psi_t^{(t, 0)}, \tilde{H}_{t-1})
    out_init <- run_psi_APF_rolling(data, t, psi_pa[t, ], H_prev_trial, model, 
                                    init = TRUE, target_N = Nt) 
    H_tilde[[t+1]] <- out_init$H
    
    #\psi_{t0:(t-1)}^{(t,0)} <- \psi_{t0:(t-1)}^{(t-1, K_{t-1})}
    psi_prev_t <- psi_pa
    xi_k_array <- numeric(max(1, Kt))
    
    if (Kt > 0) {
      for (k in 1:Kt) {
        psi_before_k <- psi_pa
        
        for (s in t:t0) {
          #set \psi_s^{(t,k)} <- learn-\psi(s, \psi_{s+1}^{(t,k)}, \tilde{H}_s)
          psi_pa[s, ] <- learn_psi(psi_pa[s+1,, drop = FALSE], H_tilde[[s+1]], model)
        }
        
        for (s in t0:t) {
          prev_s_tilde <- H_tilde[[s]] 
          
          #sample \tilde{H}_s <- \psi-apf(s, \psi_s^{(t,k)}, \tilde{H}_{s-1})
          #time 1 H[[s]] stores H_0, and H_tilde[[s]] actually stores information from previous time
          out_tilde <- run_psi_APF_rolling(data, s, psi_pa[s, ], prev_s_tilde, model, 
                                           init = FALSE, target_N = Nt)
          H_tilde[[s+1]] <- out_tilde$H
        }
        
        # Calculate metric \xi for updating \kappa
        xi_k_array[k] <- mean(abs(psi_pa[t0:t, , drop=FALSE] - psi_before_k[t0:t, , drop=FALSE]) / 
                                (abs(psi_pa[t0:t, , drop=FALSE]) + abs(psi_before_k[t0:t, , drop=FALSE]) + 1e-8))
      }
    }
    
   
    for (s in t0:t) {
      prev_s_H <- H[[s]]
      
      #sample H_s <- \psi-apf(s, \psi_s^{(t,Kt)}, H_{s-1})
      out_H <- run_psi_APF_rolling(data, s, psi_pa[s, ], prev_s_H, model, 
                                   init = FALSE, target_N = Nt)
      H[[s+1]] <- out_H$H
    }
    
    #Set \kappa_k <- 0.9\kappa_k + 0.1(\xi(...) - \kappa_k)
    if (Kt > 0) {
      for (k in 1:Kt) kappa[k] <- 0.8 * kappa[k] + 0.1 * xi_k_array[k]
    }
    CK <- which(diff(kappa[1:K_max]) > 0) + 1
    
    #set Kt+1 <- min C_K else set Kt+1 <- Kmax
    Kt_next <- if (length(CK) > 0) min(CK[CK >= K_min & CK <= K_max]) else K_max
    
    
    #Set \beta_b <- 0.9\beta_b + 0.1(\xi(...) - \beta_b)
    for (b in 1:Bt) {
      idx_s <- t - b + 1
      if (idx_s < 1) next
      xi_b <- mean(abs(psi_pa[idx_s,] - psi_prev_t[idx_s,]) / 
                     (abs(psi_pa[idx_s,]) + abs(psi_prev_t[idx_s,]) + 1e-8))
      beta[b] <- 0.8 * beta[b] + 0.1 * xi_b
    }
    CB <- which(diff(beta[1:B_max]) > 0) + 1
    
    #set Bt+1 <- min C_B else Bt+1 <- Bmax
    Bt_next <- if (length(CB) > 0) min(CB[CB <= B_max]) else B_max
    
    
    Bt_next <- min(Bt_next, t + 1)
    
    
    K_hist[t] <- Kt; B_hist[t] <- Bt; N_hist[t] <- Nt
    
    #N_{t+1} := floor(gamma / (B_{t+1} * (K_{t+1} + 1)))
    Nt_next <- floor(gamma / (Bt_next * (Kt_next + 1)))
    #Nt_next <- max(Nt_next, 150) 
   
    #approximations
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
