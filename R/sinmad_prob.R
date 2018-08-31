sinmad_flexprob <- function(K, model, features, visualise = TRUE, xlim = c(0,15), draws = 5) {
  
  mod <- as.data.frame(tidy(model))
  
  #================================================================#
  # linear predictor has intercept and is a function of covariates # 
  #================================================================#
  
  if (isTRUE("Intercept" %in% mod[,1])) {
    
    a <- -1 
    q <- -1 
    while(isTRUE((a<0) | (q<0))) {
      params <- auto_cholesky(model = model, draws = draws)
      a <- params[1]
      q <- params[2]
      Intercept <- params[3]
      betas <- params[4:length(params)]
    }
    b <- exp(Intercept + sum(features*betas))
    
    if(isTRUE(visualise)) {
      preview <<- function(x){
        dsinmad(x, scale = b, shape1.a = a, shape3.q = q)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(psinmad(K, scale = b, shape1.a = a, shape3.q = q))
    } else {
      return(psinmad(K, scale = b, shape1.a = a, shape3.q = q))
    }
    
    #=====================================#
    # mu is not a function of covariates  #
    #=====================================#
    
  } else if (!isTRUE("beta1" %in% mod[,1])) {
    
    a <- -1 
    b <- -1 
    q <- -1 
    while(isTRUE((a<0) | (b<0) | (q<0))) {
      params <- auto_cholesky(model = model, draws = draws)
      a <- params[1]
      b <- params[2]
      q <- params[3]
    }
    
    if(isTRUE(visualise)) {
      preview <<- function(x){
        dsinmad(x, scale = b, shape1.a = a, shape3.q = q)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(psinmad(K, scale = b, shape1.a = a, shape3.q = q))
    } else {
      return(psinmad(K, scale = b, shape1.a = a, shape3.q = q))
    }
    
    #====================================================================#
    # linear predictor has no intercept but is a function of covariates  # 
    #====================================================================#
    
  } else {
    
    a <- -1 
    q <- -1 
    while(isTRUE((a<0) | (q<0))) {
      params <- auto_cholesky(model = model, draws = draws)
      a <- params[1]
      q <- params[2]
      betas <- params[3:length(params)]
    }
    b <- exp(sum(betas*features))
    
    if(isTRUE(visualise)) {
      preview <<- function(x){
        dsinmad(x, scale = b, shape1.a = a, shape3.q = q)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(psinmad(K, scale = b, shape1.a = a, shape3.q = q))
    } else {
      return(psinmad(K, scale = b, shape1.a = a, shape3.q = q))
    }
  }
}