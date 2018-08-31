lomax_flexprob <- function(K, model, features, visualise = TRUE, xlim = c(0,15), draws = 5) {
  
  mod <- as.data.frame(tidy(model))
  
  #================================================================#
  # linear predictor has intercept and is a function of covariates # 
  #================================================================#
  
  if (isTRUE("Intercept" %in% mod[,1])) {
    
    alpha <- -1 
    while(isTRUE(alpha<0)) {
      params <- auto_cholesky(model = model, draws = draws)
      alpha <- params[1]
      Intercept <- params[2]
      betas <- params[3:length(params)] 
    }
    lambda <- exp(Intercept + sum(betas*features))
    
    if(isTRUE(visualise)) {
      preview <<- function(x) {
        dLomax(x, scale = lambda, shape = alpha)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(plomax(K, scale = lambda, shape = alpha))
    } else {
      return(plomax(K, scale = lambda, shape = alpha))
    }
    
    #=====================================#
    # mu is not a function of covariates  #
    #=====================================#
    
  } else if (!isTRUE("beta1" %in% mod[,1])) {
    
    alpha <- -1
    lambda <- -1
    while(isTRUE(alpha<0) | isTRUE(lambda<0)) {
      params <- auto_cholesky(model = model, draws = draws)
      alpha <- params[1]
      lambda <- params[2]
    }
    
    if(isTRUE(visualise)) {
      preview <<- function(x) {
        dLomax(x, scale = lambda, shape = alpha)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(plomax(K, scale = lambda, shape = alpha))
    } else {
      return(plomax(K, scale = lambda, shape = alpha))
    }
    
    #====================================================================#
    # linear predictor has no intercept but is a function of covariates  # 
    #====================================================================#
    
  } else {
    
    alpha <- -1 
    while(isTRUE(alpha<0)) {
      params <- auto_cholesky(model = model, draws = draws)
      alpha <- params[1]
      betas <- params[2:length(params)]
    }
    lambda <- exp(sum(betas*features))
    
    if(isTRUE(visualise)) {
      preview <<- function(x) {
        dLomax(x, scale = lambda, shape = alpha)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(plomax(K, scale = lambda, shape = alpha))
    } else {
      return(plomax(K, scale = lambda, shape = alpha))
    }
  }
}