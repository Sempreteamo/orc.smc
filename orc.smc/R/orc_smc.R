#' Rolling window function
#'
#'
#'
#' @param lag lag that user specified
#' @param data observation data
#' @param model the type of model used
#' @param N number of particles
#'
#' @return A matrix of observations with dimensions Time x d.
#' @export
Orc_SMC <- function(lag, data, model, N) {
  obs <- data$obs
  Time <- nrow(obs)
  d <- ncol(obs)

  #A <- model$tran_mu
  #B <- model$tran_cov
  #obs_params <- model$obs_params
  K <- model$parameters$k

  X0 <- matrix(0, nrow = N, ncol = d)
  w0 <- matrix(log(1/N), 1, N)

  logZ_vec <- numeric(Time)

  ess_history <- numeric(Time)

  X <- X_apf <- array(NA, c(Time, N, d))
  log_W_apf <- log_W <- matrix(NA, Time, N)
  log_likelihoods <- log_likelihoods_apf <- matrix(NA, Time, N)
  filtering_estimates <- matrix(NA, Time, d)

  psi_pa <- matrix(NA, Time + 1, 2*d)

  H <- vector('list', Time + 1)
  H_tilde <- vector('list', Time + 1)

  H[[1]] <- list(X = X0, logW = w0, logZ = 0)
  H_tilde[[1]] <- list(X = X0, logW = w0, logZ = 0)

  for (t in 1:Time) {
    print(t)
    t0 <- max(t - lag + 1, 1)

    # Step 3: Init pass with psi ≡ 1
    #psi_pa[t, ] <- rep(NA, 2*d)
    #psi_pa[t+1, ] <- rep(NA, 2*d)

    output <- run_psi_APF_rolling(data, t,  psi_pa[t, ] , H_tilde[[t]] , model, init = TRUE)
    H_tilde[[t+1]] <- output$H

    # Step 4: Policy Refinement


    for (k in 1:K) {

      for (s in t:t0) {

       psi_pa[s,] <- learn_psi(psi_pa[s+1,, drop = FALSE ], H_tilde[[s+1]], model)


      }

      for (s in t0:t) {
        output <- run_psi_APF_rolling(data, s,
                          psi_pa[s,, drop = FALSE], H_tilde[[s]], model, init = FALSE)

        H_tilde[[s+1]] <- output$H

      }
    }

    # Step 5: Final pass to get filtering distributions + logZ

    for (s in t0:t) {

      output <- run_psi_APF_rolling(data, s,
                            psi_pa[s,, drop = FALSE], H[[s]], model, init = FALSE)

      H[[s+1]] <- output$H
      if(s == t){
        ess_history[t] <- output$current_ess
      }

    }


    logZ_vec[t] <- H[[t + 1]]$logZ

  }

  return(list(logZ = logZ_vec, f_means = filtering_estimates, H_forward  = H,  ess_history = ess_history))
}
