auto_cholesky <- function(model, draws = 5) {
  
  C <- t(chol.default(vcov(model), pivot = TRUE))
  mod <- as.data.frame(tidy(model))
  
  if (draws == 0) {
    params <- as.vector(mod[,2])
    names(params) <- mod[,1]
    return(params)
  } else {
    params <- rep(0, nrow(mod))
    for (i in 1:draws) {
      params <- params + (mod[,2] + C %*% rnorm(nrow(mod)))
    }
    params <- as.vector(params)
    names(params) <- mod[,1]
    return(params/(draws))
  }
}
