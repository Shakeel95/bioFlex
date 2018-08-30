# Generalised Gamma denisty from Stacy (1962)

dStacy_gamma <- function(x, mu, alpha, phi) {
  alpha <- abs(alpha)
  alpha/gamma(phi)*mu^{alpha*phi}*x^{alpha*phi - 1}*exp(-(mu*x)^{alpha})
}

# gb2 density

dgb2 <- function(x,shape1,scale,shape2,shape3){
  y <- (x/scale)^shape1
  dy_dx <- (shape1/scale)*(x/scale)^(shape1-1)
  z <- y/(1+y)
  z[z==1] <- 1-.Machine$double.eps
  dz_dy <- (1+y)^(-2)
  dens <- dbeta(z, shape2, shape3) * dz_dy * dy_dx
  v <- (x==Inf)
  dens[v] <- 0
  return(dens)
}

pgb2 <- function(x,shape1,scale,shape2,shape3){
  y <- (x/scale)^shape1
  z <- y/(1+y)
  prob <- pbeta(z, shape2, shape3)
  v <- (x==Inf)
  prob[v] <- 1
  return(prob)
}

# Lomax density

dLomax <- function(x, scale, shape) {
  return((shape/scale)*(1+(x/scale))^{-(shape+1)})
}
