#'First stage estimate for alpha in the Lomax distribution
#'
#'Computes the method of moments estimate for the shape parameter in the Lomax distribution.
#'@param x A numeric vector of i.i.d. observations presumed to be draws from a Lomax distribution.
#'@export


lomax_shape <- function(x) {

  # Usual estimation procedure
  if (var(x) >= mean(x)^{2}) {
    # Usual estimation procedure
    m = var(x)/mean(x)^{2}
    shape = 2*m/(m-1)
    return(shape)

    # Use cscale estimate when usual estimation procedure fails
  } else if (var(x) < mean(x)^{2}) {
    scale = lomax_scale(x)
    return(scale/mean(x) + 1)
  }
}

#'First stage estimate for lambda in the Lomax distribution
#'
#'Computes the method of moments estimate for the scale parameter in the Lomax distribution.
#'@param x A numeric vector of i.i.d. observations presumed to be draws from a Lomax distribution.
#'@export

lomax_scale <- function(x) {

  # Usual estimation procedure
  if (var(x) >= mean(x)^{2}) {
    shape = lomax_shape(x)
    scale = mean(x)*(shape-1)
    return(abs(scale))

    # Assume a = 2 when usual procedure fails
  } else if (var(x) < mean(x)^{2}) {
    return(median(x)/(sqrt(2)-1))
  }
}
