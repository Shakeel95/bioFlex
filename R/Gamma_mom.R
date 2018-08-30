#' First stage estimate for alpha in the Gamma distribution
#'
#' Computes the method of moments estimate for the shape parameter in Gamma distribution.
#' @param x A numeric vector of i.i.d. observations presumed to be draws from a Gamma distribution.
#' @export

gamma_shape <- function(x) {
  return(mean(x)^2/var(x))
}


#' First stage estimate for beta in the Gamma distribution
#'
#' Computes the method of moments estimate for the rate parameter in Gamma distribution.
#' @param x A numeric vector of i.i.d. observations presumed to be draws from a Gamma distribution.
#' @export
gamma_rate <- function(x) {
  return(mean(x)/var(x))
}
