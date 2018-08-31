lnorm_flexprob <- function(K, model, features, visualise = TRUE, xlim = c(0,15), draws = 5) {
  
  mod <- as.data.frame(tidy(model))
  
  #================================================================#
  # linear predictor has intercept and is a function of covariates # 
  #================================================================#
  
  if (isTRUE("Intercept" %in% mod[,1])) {
    
    sigma <- -1 
    while(isTRUE(sigma <0)) {
      params <- auto_cholesky(model = model, draws = draws)
      sigma <- params[1]
      Intercept <- params[2]
      betas <- params[3:length(params)]
    }
    mu <- Intercept + sum(betas*features)
    
    if (isTRUE(visualise)) {
      preview <<- function(x) {
        dlnorm(x, meanlog = mu, sdlog = sigma)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(plnorm(K, meanlog = mu, sdlog = sigma))
    } else {
      return(plnorm(K, meanlog = mu, sdlog = sigma))
    }
    
    #=====================================#
    # mu is not a function of covariates  #
    #=====================================# 
    
  } else if (!isTRUE("beta1" %in% mod[,1])) {
    
    sigma <- -1
    while(isTRUE(sigma < 0)) {
      params <- auto_cholesky(model = model, draws = draws)
      sigma <- params[1]
      mu <- params[2]
    }
    
    if (isTRUE(visualise)) {
      preview <<- function(x) {
        dlnorm(x, meanlog = mu, sdlog = sigma)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(plnorm(K, meanlog = mu, sdlog = sigma))
    } else {
      return(plnorm(K, meanlog = mu, sdlog = sigma))
    }
    
    #====================================================================#
    # linear predictor has no intercept but is a function of covariates  # 
    #====================================================================#
    
  } else {
    
    sigma <- -1 
    while(isTRUE(sigma<0)) {
      params <- auto_cholesky(model = model, draws = draws)
      sigma <- params[1]
      betas <- params[2:length(params)]
    }
    mu <- sum(betas*features)
    
    if (isTRUE(visualise)) {
      preview <<- function(x) {
        dlnorm(x, meanlog = mu, sdlog = sigma)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(plnorm(K, meanlog = mu, sdlog = sigma))
    } else {
      return(plnorm(K, meanlog = mu, sdlog = sigma))
    }
    
  }
}