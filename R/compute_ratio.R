#' Function to compute ratio between the normalising constant and the estimates
#'
#' @param Z Estimate results given by the algorithm
#' @param fkf.obj Normalising constant of the Gaussian distributions computed by FKF
#'
#' @return Log ratio between real values and estimates
#' @export
#'
compute_ratio <- function(Z, fkf.obj){
  fkf.obj_Z <- fkf.obj$logLik
  cat('NC = ', fkf.obj_Z, 'est = ', Z, 'ratio = ', exp(Z-fkf.obj_Z))
  return(exp(Z-fkf.obj_Z))
}
