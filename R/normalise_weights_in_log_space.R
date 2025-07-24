#' Function to calculate the normalised weights and the normalising constant
#'
#'This function takes a vector of unnormalised weights and normalises them in the usual numerically
#'stable way. It returns the normalised weights as well as the normalising constant.
#'The inputs and outputs are all in log-space.
#'
#' @param W_un Vector of unnormalised weights
#'
#' @return A list with two elements:
#' logW is log of the vector of self-normalised weights
#' logZ is log of the normalising constant, i.e. of the sum of the unnormalised weights
#'
#' @export
#'
normalise_weights_in_log_space <- function(W_un){
  mx <- max(W_un)
  w_centered <- exp(W_un - mx)
  sum_wc <- sum(w_centered)
  logW <- w_centered/sum_wc
  logZ <- log(sum_wc) - log(length(W_un)) + mx
  return(list(logW, logZ))
}
