gb2_flexprob <- function(K, model, features, visualise = TRUE, xlim = c(0,15), draws = 5) {
  
  mod <- as.data.frame(tidy(model))
  
  #================================================================#
  # linear predictor has intercept and is a function of covariates # 
  #================================================================#
  
  if (isTRUE("Intercept" %in% mod[,1])) {
    
    a <- -1 
    p <- -1 
    q <- -1 
    while(isTRUE((a<0) | (p<0) | (q<0))) {
      params <- auto_cholesky(model = model, draws = draws)
      a <- params[1]
      p <- params[2]
      q <- params[3]
      Intercept <- params[4]
      betas <- params[5:length(params)]
    }
    b <- exp(Intercept + sum(betas*features))
    
    if(isTRUE(visualise)) {
      preview <<- function(x) {
        dgb2(x, scale = b, shape1 = a, shape2 = p, shape3 = q)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(pgb2(K, scale = b, shape1 = a, shape2 = p, shape3 = q))
    } else {
      return(pgb2(K, scale = b, shape1 = a, shape2 = p, shape3 = q))
    }
    
    #=====================================#
    # mu is not a function of covariates  #
    #=====================================#
    
  } else if (!isTRUE("beta1" %in% mod[,1])) {
    
    a <- -1 
    b <- -1 
    p <- -1 
    q <- -1 
    
    while(isTRUE((a<0) | (b<0) | (p<0) | (q<0))) {
      params <- auto_cholesky(model = model, draws = draws)
      a <- params[2]
      b <- params[1]
      p <- params[3]
      q <- params[4]
    }
    
    if(isTRUE(visualise)) {
      preview <<- function(x) {
        dgb2(x, scale = b, shape1 = a, shape2 = p, shape3 = q)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(pgb2(K, scale = b, shape1 = a, shape2 = p, shape3 = q))
    } else {
      return(pgb2(K, scale = b, shape1 = a, shape2 = p, shape3 = q))
    }
    
    #====================================================================#
    # linear predictor has no intercept but is a function of covariates  # 
    #====================================================================#
    
  } else {
    
    a <- -1 
    p <- -1 
    q <- -1 
    while(isTRUE((a<0) | (p<0) | (q<0))) {
      params <- auto_cholesky(model = model, draws = draws)
      a <- params[1]
      p <- params[2]
      q <- params[3]
      betas <- params[4:length(params)]
    }
    b <- exp(sum(betas*features))
    
    if(isTRUE(visualise)) {
      preview <<- function(x) {
        dgb2(x, scale = b, shape1 = a, shape2 = p, shape3 = q)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(pgb2(K, scale = b, shape1 = a, shape2 = p, shape3 = q))
    } else {
      return(pgb2(K, scale = b, shape1 = a, shape2 = p, shape3 = q))
    }
  }
}