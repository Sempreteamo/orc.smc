#' Function to evaluate a twisting function psi
#'
#'This function evaluates the Gaussian density of
#'a twisting function with the specified parameters at one or more points.

#' @param x Points to evalute the function at
#' @param psi_pa Parameters of twisting function
#'
#' @return Psi function evaluated at specified points
#' @export
#'
#'
evaluate_psi <- function(x, psi_pa){
  if(all(is.na(psi_pa))){
    psi_x <- 0
  }else{
    d <- length(x)
    dif <- x - psi_pa[1:d]
    psi_x <- -(d/2)*log(2*pi) - (1/2)*log(det(diag(psi_pa[(d+1):(d+d)], nrow=d,ncol=d))) +
      (-1/2)*t(dif)%*%diag((psi_pa[(d+1):(d+d)])^(-1), nrow=d,ncol=d)%*%dif
    
  }
  
 
  return(psi_x)
}

