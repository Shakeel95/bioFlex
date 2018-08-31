fisk_flexprob <- function(K, model, features, visualise = TRUE, xlim = c(0,15), draws = 5) {
  
  mod <- as.data.frame(tidy(model))
  
  #================================================================#
  # linear predictor has intercept and is a function of covariates # 
  #================================================================#
  
  if (isTRUE("Intercept" %in% mod[,1])) {
    
    beta <- -1 
    while(isTRUE(beta<0)){
      params <- auto_cholesky(model = model, draws = draws)
      beta <- params[1]
      Intercept <- params[2]
      betas <- params[3:length(params)]
    }
    alpha <- exp(Intercept + sum(betas*features))
    
    if(isTRUE(visualise)) {
      preview <<- function(x) {
        dllogis(x, scale = alpha, shape = beta)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(pllogis(K, scale = alpha, shape = beta))
    } else {
      return(pllogis(K, scale = alpha, shape = beta))
    }
    
    #=====================================#
    # mu is not a function of covariates  #
    #=====================================# 
    
  } else if (!isTRUE("beta1" %in% mod[,1])) {
    
    beta <- -1
    alpha <- -1
    while(isTRUE(beta<0) | isTRUE(alpha<0)) {
      params <- auto_cholesky(model = model, draws = draws)
      alpha <- params[1]
      beta <- params[2]
    }
    
    if(isTRUE(visualise)) {
      preview <<- function(x) {
        dllogis(x, scale = alpha, shape = beta)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(pllogis(K, scale = alpha, shape = beta))
    } else {
      return(pllogis(K, scale = alpha, shape = beta))
    }
    
    #====================================================================#
    # linear predictor has no intercept but is a function of covariates  # 
    #====================================================================#
    
  } else {
    
    beta <- -1 
    while(isTRUE(beta<0)) {
      params <- auto_cholesky(model = model, draws = draws)
      beta <- params[1]
      betas <- params[2:length(params)]
    }
    alpha <- exp(sum(betas*features))
    
    if(isTRUE(visualise)) {
      preview <<- function(x) {
        dllogis(x, scale = alpha, shape = beta)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(pllogis(K, scale = alpha, shape = beta))
    } else {
      return(pllogis(K, scale = alpha, shape = beta))
    }
  }
}
  