#'First stage estimate for alpha in the Generalized Gamma distribution
#'
#'Computes the method of moments estimate for the first shape parameter in the Generalized Gamma dustribution. The parameterization proposed by Stacey (1962) is assumed.
#'@param x A numeric vector of i.i.d. observations presumed to be draws from a Generalized Gamma distribution.
#'@param lower Lower bound of the solution space in which the estimate may lie. By default the parameter is set to zero, since alpha is strictly positive.
#'@param upper Upper bound of the solution space in which the estimate may lie. By default the parameter is set to 100, although this choice is somewhat arbitrary as the estimate is generally small and the function defining the parameter space may only have one root.
#'@export
#'@references Alberto Achcar, Jorge, Pedro Luiz Ramos, and Edson Zangiacomi Martinez. "Some Computational Aspects To Find Accurate Estimates For The Parameters Of The Generalized Gamma Distribution." Pesquisa Operacional Vol.37(2) (2017): n. pag. Print.

ggamma_alpha <- function(x, lower = 0, upper = 100){
  # define k function
  k <- function(a) {
    (length(x) - (3/2) -a/(1+a))*(1/a) - log(a)*a/(1+a)^{2}
  }
  # define phi function
  phi <- function(a) {
    k(a)*sum(x^{a})/(length(x)*sum(x^{a}*log(x^{a})) - sum(x^{a})*sum(log(x^{a})))
  }
  # define h function
  h <- function(a) {
    length(x)*digamma(length(x)*phi(a)) + a*sum(log(x)) - (1/phi(a)) - digamma(phi(a)) - log(sum(x^{a}))
  }
  # Find minimum of function
  optim <- optimize(h, interval = c(lower , upper), maximum = FALSE)
  return(optim$minimum)
}

#'First stage estimate for phi in the Generalized Gamma distribution
#'
#'Computes the method of moments estimate for the second shape parameter in the Generalized Gamma dustribution. The parameterization proposed by Stacey (1962) is assumed.
#'@param x A numeric vector of i.i.d. observations presumed to be draws from a Generalized Gamma distribution.
#'@export
#'@references Alberto Achcar, Jorge, Pedro Luiz Ramos, and Edson Zangiacomi Martinez. "Some Computational Aspects To Find Accurate Estimates For The Parameters Of The Generalized Gamma Distribution." Pesquisa Operacional Vol.37(2) (2017): n. pag. Print.


ggamma_phi <- function(x) {
  # return estimated alpha
  a <- ggamma_alpha(x)
  # define k function
  k <- function(x) {
    (length(x) - (3/2) - a/(1+a))*(1/a) - log(a)*a/(1+a)^{2}
  }
  # solve for phi
  phi = k(x)*sum(x^{a})/(length(x)*sum(x^{a}*log(x^{a})) - sum(x^{a})*sum(log(x^{a})))
  return(phi)
}

#'First stage estimate for mu in the Generalized Gamma distribution
#'
#'Computes the method of moments estimate for the scale parameter in the Generalized Gamma dustribution. The parameterization proposed by Stacey (1962) is assumed; in many popular parameterizations the value returned is the inverse scale parameter.
#'@param x A numeric vector of i.i.d. observations presumed to be draws from a Generalized Gamma distribution.
#'@export
#'@references Alberto Achcar, Jorge, Pedro Luiz Ramos, and Edson Zangiacomi Martinez. "Some Computational Aspects To Find Accurate Estimates For The Parameters Of The Generalized Gamma Distribution." Pesquisa Operacional Vol.37(2) (2017): n. pag. Print.


ggamma_mu <- function(x) {
  # define alpha and phi
  a = ggamma_alpha(x)
  p = ggamma_phi(x)
  # greater than one
  if (a*p*length(x) >= 1) {
    u = ((length(x)*a*p - 1)/(a*sum(x^{a})))^{(1/a)}
  } else if ( a*p*length(x) < 1) {
    u = ((length(x)*p)/(sum(x^{a})))^{(1/a)}
  }
  return(u)
}
