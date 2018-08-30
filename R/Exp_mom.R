#' First stage estimate for lambda in the Exponential distribution
#'
#' Computes the the method of moments estimate for the rate parameter in the Exponential distribution.
#'
#' @param x A numeric vector of i.i.d. observations presumed to be draws from an Exponential distribution.
#' @export


exp_rate <- function(x) {
  if (mean(x) >= 0) {
    return(mean(x)^{-1})
  } else {
    return(1)
  }
}
