#'First stage estimate for mu in the Log-Normal distribution
#'
#'Computes the method of moments estimate for the effective location parameter in the Log-Normal distribution.
#'@param x A numeric vector of i.i.d. observations presumed to be draws from a Log-Normal distribution.
#'@export
#'@references Faith Ginos, Brenda. "Parameter Estimation For The Lognormal Distribution." Brigham Young University Scholars Archive (2018): 1-111. Web. 10 Aug. 2018.


lnorm_mean <- function(x) {
  x <- c(as.numeric(x))
  return(-0.5*log(sum(x)) + 2*log(sum(x)) - (3/2)*log(length(x)))
}

#'First stage estimate for sigma squared in the Log-Normal distribution
#'
#'Computes the method of moments estimate for the effective dispersion parameter in the Log-Normal distru.
#'@param x A numeric vector of i.i.d. observations presumed to be draws from a Log-Normal distribution.
#'@export
#'@references Faith Ginos, Brenda. "Parameter Estimation For The Lognormal Distribution." Brigham Young University Scholars Archive (2018): 1-111. Web. 10 Aug. 2018.


lnorm_var <- function(x) {
  x <- c(as.numeric(x))
  return(log(sum(x^2))-2*log(sum(x))+log(length(x)))
}

