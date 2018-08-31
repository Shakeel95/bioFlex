exp_flexprob <- function(K, model, features, visualise = TRUE, xlim = c(0,15), draws = 5) {
  
  mod <- as.data.frame(tidy(model))
  
  #================================================================#
  # linear predictor has intercept and is a function of covariates # 
  #================================================================#
  
  if (isTRUE("Intercept" %in% mod[,1])) {
    
    params <- auto_cholesky(model = model, draws = draws)
    Intercept <- params[1]
    betas <- params[2:length(params)]
    lambda <- exp(Intercept + sum(features*betas))
    
    if(isTRUE(visualise)) {
      preview <<- function(x) {
        dexp(x, rate = lambda)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(pexp(K, rate = lambda))
    } else {
      return(pexp(K, rate = lambda))
    }
    
    #=====================================#
    # mu is not a function of covariates  #
    #=====================================#
    
  } else if (!isTRUE("beta1" %in% mod[,1])) {
    
    lambda <- -1 
    while(isTRUE(lambda<0)) {
      lambda <- auto_cholesky(model = model, draws = draws)
    }
    
    if(isTRUE(visualise)) {
      preview <<- function(x) {
        dexp(x, rate = lambda)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(pexp(K, rate = lambda))
    } else {
      return(pexp(K, rate = lambda))
    }
    
    #====================================================================#
    # linear predictor has no intercept but is a function of covariates  # 
    #====================================================================#
    
  } else {
    
    betas <- auto_cholesky(model = model, draws = draws)
    lambda <- exp(sum(betas*features))
    
    if(isTRUE(visualise)) {
      preview <<- function(x) {
        dexp(x, rate = lambda)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(pexp(K, rate = lambda))
    } else {
      return(pexp(K, rate = lambda))
    }
  }
}
  