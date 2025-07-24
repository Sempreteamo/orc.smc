#' Function to sample observations
#'
#'This function samples synthetic observations for each time step given the state space model.
#'
#' @param model List containing model parameters
#' @param Time Total time length
#' @param d The dimension of observations
#'
#' @return A list containing the observation sequence
#' @export
sample_obs <- function(model, Time, d){
    ini <- model$ini_mu
    ini_c <- model$ini_cov
    tran_mu <- model$tran_mu
    #rep(as.matrix(model$tran_mu)[1,1], d)
    tran_cov <- model$tran_c

    X <- matrix(0, nrow = Time, ncol = d)
    data <- matrix(0, nrow = Time, ncol = d)

    X[1,] <- mvnfast::rmvn(1, ini, ini_c)

    if(Time >= 2){
      for(t in 2:Time){
        X[t,] <- mvnfast::rmvn(1, tran_mu%*%X[t-1,], tran_cov)
        #stats::rnorm(d, ini + tran_mu*(X[t-1,] - ini), sqrt(as.matrix(tran_cov)[1,1]))
      }
    }


    for(t in 1:Time){
      data[t,] <- model$simu_observation(X[t,], model$obs_params)
    

    }

  return(data)
}

#' @import mvnfast
#' @import stats
