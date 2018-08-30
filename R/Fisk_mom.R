#'First stage estimate for alpha in the Fisk distribution
#'
#'Computes the method of moments estimate for the scale parameter in the Fisk dustribution.
#'@param x A numeric vector of i.i.d. observations presumed to be draws from a Fisk distribution.
#'@export
#'@references "Solving Log-Logistic Distribution Parameters From Moments." Cross Validate. N.p., 2017. Web. 11 Aug. 2018.

fisk_scale <- function(x) {
  return(median(x))
}

#'First stage estimate for beta in the Fisk distribution
#'
#'Computes the method of moments estimate for the shape parameter in the Fisk dustribution.
#'@param x A numeric vector of i.i.d. observations presumed to be draws from a Fisk distribution.
#'@param lower_bound Lower bound in solution pace of the trigonometric equation used to find the value of beta. Ideally zero since as the solution will always be between zero and pi/2, however since the function must be defined at both endpoints zero cannot be used. By default the parameter is set to 0.001.
#'@param upper_bound Upper bound in solution pace of the trigonometric equation used to find the value of beta.  By default the parameter is set to pi/2.
#'@export
#'@references "Solving Log-Logistic Distribution Parameters From Moments." Cross Validate. N.p., 2017. Web. 11 Aug. 2018.

fisk_shape <- function(x, lower_bound = 0.001, upper_bound = pi/2) {

  # obtain E(x)^2/E(x^2), which is a function of the shape parameter only
  m = mean(x)^{2}/mean(x^2)

  # define a function whose root is the shape parameter
  fun <- function(theta) {
    return(theta/tan(theta) - m)
  }

  # find the root
  roots <- uniroot(fun, lower = lower_bound, upper = upper_bound)

  # solve for shape parameter and return
  beta = pi/roots$root
  return(beta)
}
