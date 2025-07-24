
#'
#'
#' Log-Sum-Exp function to prevent underflow
#'
#' @param log_weights the log of weights of particles
#'
#' @return A list containing:
#'
#'
#' @export
#'
log_sum_exp <- function(log_weights) {
  max_logW <- max(log_weights)
  w_centered <- exp(log_weights - max_logW)
  sum_wc <- sum(w_centered)
  log_w <- max_logW + log(sum_wc)
  return(log_w)
}

# Main function implementing ψ-APF
