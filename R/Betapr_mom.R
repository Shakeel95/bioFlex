#'First stage estimate for beta in the Beta Prime distribution
#'
#'Computes the method of moments estimate for the second shape parameter in the Beta Prime dustribution.
#'@param x A numeric vector of i.i.d. observations presumed to be draws from a Beta Prime distribution.
#'@export


betapr_shape2 <- function(x) {
  # Bits of quadratic
  a = var(x)/mean(x)
  b = -(3*a + mean(x) +1)
  c = 2*a + mean(x) + 1
  # Solve quadratic
  if (b^2 >= 4*a*c) {
    return((-b+sqrt(b^2 - 4*a*c))/(2*a))
  } else {
    return(1)
  }
}

#'First stage estimate for alpha in the Beta Prime distribution
#'
#'Computes the method of moments estimate for the first shape parameter in the Beta Prime dustribution.
#'@param x A numeric vector of i.i.d. observations presumed to be draws from a Beta Prime distribution.
#'@export

betapr_shape1 <- function(x) {
  # Define shape 2
  shape2 = betapr_shape2(x)
  # Solve for shape 2
  if (shape2 > 1){
    return(mean(x)*(shape2-1))
  } else if (shape2 < 1) {
    return(abs(mean(x)*(shape2-1)))
  } else {
    return(median(x)*(shape2 +1) + 1)
  }
}
