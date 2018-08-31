dagum_flexprob <- function(K, model, features, visualise = TRUE, xlim = c(0,15), draws = 5) {
  
  mod <- as.data.frame(tidy(model))
  
  #================================================================#
  # linear predictor has intercept and is a function of covariates # 
  #================================================================#
  
  if (isTRUE("Intercept" %in% mod[,1])) {
    
    a <- -1 
    p <- -1 
    while(isTRUE((a<0) | (p<0))) {
      params <- auto_cholesky(model = model, draws = draws)
      a <- params[1]
      p <- params[2]
      Intercept <- params[3]
      betas <- params[4:length(params)]
    }
    b <- exp(Intercept + sum(betas*features))
    
    if(isTRUE(visualise)){
      preview <<- function(x) {
        ddagum(x, scale = b, shape1.a = a, shape2.p = p)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(pdagum(K, scale = b, shape1.a = a, shape2.p = p))
    } else {
      return(pdagum(K, scale = b, shape1.a = a, shape2.p = p))
    }
    
    #=====================================#
    # mu is not a function of covariates  #
    #=====================================#
    
  } else if (!isTRUE("beta1" %in% mod[,1])) {
    
    a <- -1 
    b <- -1 
    p <- -1 
    while(isTRUE((a<0) | (b<0) | (p<0))) {
      params <- auto_cholesky(model = model, draws = draws)
      a <- params[1]
      b <- params[2]
      p <- params[3]
    }
    
    if(isTRUE(visualise)){
      preview <<- function(x) {
        ddagum(x, scale = b, shape1.a = a, shape2.p = p)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(pdagum(K, scale = b, shape1.a = a, shape2.p = p))
    } else {
      return(pdagum(K, scale = b, shape1.a = a, shape2.p = p))
    }
    
    #====================================================================#
    # linear predictor has no intercept but is a function of covariates  # 
    #====================================================================#
    
  } else {
    
    a <- -1 
    p <- -1 
    while(isTRUE((a<0) | (p<0))) {
      params <- auto_cholesky(model = model, draws = draws)
      a <- params[1]
      p <- params[2]
      betas <- params[3:length(params)]
    }
    b <- exp(sum(features*betas))
    
    if(isTRUE(visualise)){
      preview <<- function(x) {
        ddagum(x, scale = b, shape1.a = a, shape2.p = p)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(pdagum(K, scale = b, shape1.a = a, shape2.p = p))
    } else {
      return(pdagum(K, scale = b, shape1.a = a, shape2.p = p))
    }
  }
}