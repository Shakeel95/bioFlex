weibull_flexprob <- function(K, model, features, visualise = TRUE, xlim = c(0,15), draws = 5) {
  
  mod <- as.data.frame(tidy(model))
  
  #================================================================#
  # linear predictor has intercept and is a function of covariates # 
  #================================================================#
  
  if (isTRUE("Intercept" %in% mod[,1])) {
    
    kappa <- -1
    while(isTRUE(kappa<0)) {
      params <- auto_cholesky(model = model, draws = draws)
      kappa <- params[1]
      Intercept <- params[2]
      betas <- params[3:length(params)]
    }
    lambda <- exp(Intercept + sum(betas*features))
    
    if(isTRUE(visualise)) {
      preview <<- function(x) {
        dweibull(x, shape = kappa, scale = lambda)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(pweibull(K, shape = kappa, scale = lambda))
    } else {
      return(pweibull(K, shape = kappa, scale = lambda))
    }
    
    #=====================================#
    # mu is not a function of covariates  #
    #=====================================# 
    
  } else if (!isTRUE("beta1" %in% mod[,1])) {
    
    kappa <- -1
    lambda <- -1
    while(isTRUE(kappa<0) | isTRUE(lambda<0)) {
      params <- auto_cholesky(model = model, draws = draws)
      kappa <- params[2]
      lambda <- params[1]
    }
    
    if(isTRUE(visualise)) {
      preview <<- function(x) {
        dweibull(x, shape = kappa, scale = lambda)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(pweibull(K, shape = kappa, scale = lambda))
    } else {
      return(pweibull(K, shape = kappa, scale = lambda))
    }
    
    #====================================================================#
    # linear predictor has no intercept but is a function of covariates  # 
    #====================================================================#
    
  } else {
    
    kappa <- -1 
    while(isTRUE(kappa<0)) {
      params <- auto_cholesky(model = model, draws = draws)
      kappa <- params[1]
      betas <- params[2:length(params)]
    }
    lambda <- exp(sum(betas*features))
    
    if(isTRUE(visualise)) {
      preview <<- function(x) {
        dweibull(x, shape = kappa, scale = lambda)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(pweibull(K, shape = kappa, scale = lambda))
    } else {
      return(pweibull(K, shape = kappa, scale = lambda))
    }
  }
}
  