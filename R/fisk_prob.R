#' Conditional probability for Fisk distribution
#' 
#' fisk_flexprob returns the conditional probability P(Y<k | X) of a model fitted via the function fisk_flexfit; where α has been specified to be a function of covariates the required value should be specified using the ‘features’ parameter. The function includes a procedure for visualizing the conditional probability. fisk_flexprob also allows for the correlation of estimated parameters via the Cholesky decomposition of the variance-covariance matrix.
#'@param K Value for which P(Y<k | X) is computed.
#'@param model An object of class "mle2" produced using the function fisk_flexfit.
#'@param features A numeric vector specifying the value of covriates at which the conditional probability should be evaluated; the covariates in the vector should appear in the same order as they do in the model. Where a model does not depend on covariates the argument may be left blank.
#'@param visualise Logical. If TRUE (the default) the conditional distribution is plotted at P(Y<k | x) is shaded. 
#'@param xlim Numeric vectors of length 2, giving the coordinate range of the dependent variable.
#'@details This function uses the same parametrization of the Fisk distribution as is used in Kleiber and Kotz (2003). The probability probability density function is used is:
#'@details f(y) = αy^α-1/[β^α(1+(y/β)^α]^2
#'@details The function returns:
#'@details P(Y>k | X) = [1+(k/α)^{-β}]^{-1}
#'@details α may be a function of covariates; in which case, the cannonical log link function is used.
#'@references Kleiber, Christian, and Samuel Kotz. Statistical Size Distributions In Economics And Actuarial Sciences. pp. 107-147. John Wiley & Sons, 2003.
#'@export 

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
  