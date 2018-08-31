ggamma_flexprob <- function(K, model, features, visualise = TRUE, xlim = c(0,15), draws = 5) {
  
  mod <- as.data.frame(tidy(model))
  
  #================================================================#
  # linear predictor has intercept and is a function of covariates # 
  #================================================================#
  
  if (isTRUE("Intercept" %in% mod[,1])) {
    alpha <- -1 
    phi <- -1 
    while(isTRUE((alpha<0) | (phi<0))) {
      params <- auto_cholesky(model = model, draws = draws)
      alpha <- params[1]
      phi <- params[2]
      Intercept <- params[3]
      betas <- params[4:length(params)]
    }
    mu <- exp(Intercept + sum(betas*features))^{-1}
    
    if (isTRUE(visualise)) {
      preview <<- function(x) {
        dStacy_gamma(x, mu = mu, alpha = alpha, phi = phi)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      ans <- integrate(preview, lower = 0, upper = K)
      return(ans)
    } else {
      preview <<- function(x) {
        dStacy_gamma(x, mu = mu, alpha = alpha, phi = phi)
      }
      ans <- integrate(preview, lower = 0, upper = K)
      return(ans)
    }
    
    #=====================================#
    # mu is not a function of covariates  #
    #=====================================# 
    
  } else if (!isTRUE("beta1" %in% mod[,1])) {
    
    alpha <- -1 
    phi <- -1 
    mu <- -1 
    while(isTRUE((alpha<0) | (phi<0) | (mu<0))) {
      params <- auto_cholesky(model = model, draws = draws)
      alpha <- params[1]
      mu <- params[2]
      phi <- params[3]
    }
    
    if (isTRUE(visualise)) {
      preview <<- function(x) {
        dStacy_gamma(x, mu = mu, alpha = alpha, phi = phi)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      ans <- integrate(preview, lower = 0, upper = K)
      return(ans)
    } else {
      preview <<- function(x) {
        dStacy_gamma(x, mu = mu, alpha = alpha, phi = phi)
      }
      ans <- integrate(preview, lower = 0, upper = K)
      return(ans)
    }
    
    #====================================================================#
    # linear predictor has no intercept but is a function of covariates  # 
    #====================================================================#
    
  } else {
    
    alpha <- -1 
    phi <- -1 
    while(isTRUE((alpha<0) | (phi<0))) {
      params <- auto_cholesky(model = model, draws = draws)
      alpha <- params[1]
      phi <- params[2]
      betas <- params[3:length(params)]
    }
    mu <- exp(sum(features*betas))^{-1}
    
    if (isTRUE(visualise)) {
      preview <<- function(x) {
        dStacy_gamma(x, mu = mu, alpha = alpha, phi = phi)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      ans <- integrate(preview, lower = 0, upper = K)
      return(ans)
    } else {
      preview <<- function(x) {
        dStacy_gamma(x, mu = mu, alpha = alpha, phi = phi)
      }
      ans <- integrate(preview, lower = 0, upper = K)
      return(ans)
    }
  }
}
  