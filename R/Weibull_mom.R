#'First stage estimate for lambda in the Weibull distribution
#'
#'Computes the method of moments estimate for the scale parameter in the Weibull dustribution.
#'@param x A numeric vector of i.i.d. observations presumed to be draws from a Weibull distribution.
#'@export
#'@references NIOLA, VINCENZO, ROSARIO OLIVIERO, and GIUSEPPE QUAREMBA. "A Method Of Moments For The Estimation Of Weibull Pdf Parameters." Proceedings of the 8th WSEAS Int. Conference on Automatic Control, Modeling and Simulation 8 (2014)


weibull_scale <- function(x) {

  # Tidy up data
  x <- c(as.numeric(x))

  # Generate empirical CDF
  CDF <- ecdf(x)

  # First step
  R = log(1-CDF(x))/x

  # Check for silly values
  for (i in 1:length(R)) {
    if (R[i] == Inf | R[i] == -Inf) {
      R[i] <- 0
    } else {
      R[i] <- R[i]
    }
  }

  # Second step
  A = -(length(x))^{-1}*sum(R)

  # Output scale parameter
  return(A^{-1})
}

#'First stage estimate for kappa in the Weibull distribution
#'
#'Computes the method of moments estimate for the shape parameter in the Weibull dustribution.
#'@param x A numeric vector of i.i.d. observations presumed to be draws from a Weibull distribution.
#'@export
#'@references NIOLA, VINCENZO, ROSARIO OLIVIERO, and GIUSEPPE QUAREMBA. "A Method Of Moments For The Estimation Of Weibull Pdf Parameters." Proceedings of the 8th WSEAS Int. Conference on Automatic Control, Modeling and Simulation 8 (2014)


weibull_shape <- function(x) {

  # Tidy up data
  x <- c(sort(as.numeric(x)))
  vec1 <- x[c(2:length(x))]
  vec2 <- x[c(1:length(x)-1)]

  # Define empirical CDF
  CDF <- ecdf(x)

  # LHS log operation, and chech for Inf
  LHS <- log(log((1-CDF(vec1))^{-1}))
  for (i in 1:length(LHS)) {
    if (LHS[i] == Inf | LHS[i] == -Inf) {
      LHS[i] <- 0
    } else {
      LHS[i] <- LHS[i]
    }
  }

  # RHS log operation, and chech for Inf
  RHS <- log(log((1-CDF(vec2))^{-1}))
  for (i in 1:length(RHS)) {
    if (RHS[i] == Inf | RHS[i] == -Inf) {
      RHS[i] <- 0
    } else {
      RHS[i] <- RHS[i]
    }
  }

  # Sum log opperations
  R = LHS - RHS
  R = sum(R)

  # Output answer
  ans = (R/(max(x)-min(x)))

  if (R >= 0) {
    return(R)
  } else {
    print("Warning: moment estimate of scale parameter is not positive, Weibull distribution may not be appropriate!")
    return(abs(R))
  }
}
