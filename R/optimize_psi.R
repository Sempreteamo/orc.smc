#' Function to optimize the Gaussian twisting function parameters
#'
#'This function takes a collection of particle locations and the associated value of the log target
#'function evaluated at those points and fits a Gaussian twisting function via the controlled SMC approach.
#'
#' @param x Locations
#' @param lfn Log function value
#'
#' @return Twisting function parameters
#' @export
#'
optimize_psi <- function(x, lfn){
  x <- as.matrix(x)
  d <- ncol(x)
  
  valid_idx <- is.finite(lfn) & apply(is.finite(x), 1, all)
  
  if (sum(valid_idx) == 0) {
    return(c(rep(0, d), rep(1e10, d))) 
  }
  
  x_valid <- x[valid_idx, , drop = FALSE]
  lfn_valid <- lfn[valid_idx]
  
  fit_df <- data.frame(x2 = x_valid^2, x = x_valid)
  
  fit <- try(stats::lm(lfn_valid ~ ., data = fit_df), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(c(rep(0, d), rep(1e10, d))) 
  }
  
  coef <- -fit$coefficients
  
  a <- coef[2:(1+d)]
  b <- coef[(2+d):length(coef)]
  
  if (any(a <= 0, na.rm = TRUE)) {
    warning("optimize_psi: a <= 0 detected, clamped to 1e-10")
    a[a <= 0] <- 1e-10
  }
  
  params <- c(b/(-2*a), 1/(2*a))
  
  return(params)
}
#' @import stats
