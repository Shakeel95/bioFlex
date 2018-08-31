#' Conditional mean for Exponential distribution
#' 
#' exp_flexpredict returns the conditional mean E(Y|X) of a model fitted via the function exp_flexfit; where λ has been specified to be a function of covariates the required value should be specified using the ‘features’ parameter. exp_flexpredict also allows for the correlation of estimated parameters via the Cholesky decomposition of the variance-covariance matrix.
#'@param model An object of class "mle2" produced using the function exp_flexfit.
#'@param features  A numeric vector specifying the value of covriates at which the conditional mean should be evaluated; the covariates in the vector should appear in the same order as they do in the model. Where a model does not depend on covariates the argument may be left blank. 
#'@param draws The number of random draws from multivariate random normal representing correlated parameters. If parameter correlation is not required draws should be set to zero. 
#'@details This function uses the most common parametrization of the Exponential distribution. The probability probability density function is used is:
#'@details f(y) = λexp(-λy)
#'@details The function returns: 
#'@details E(Y|X) = λ^{-1}
#'@details λ may be a function of covariates; in which case, the canonical log link function is used.  


exp_flexpredict <- function(model, features, draws = 5) {
  
  mod <- as.data.frame(tidy(model))
  
  #================================================================#
  # linear predictor has intercept and is a function of covariates # 
  #================================================================#
  
  if (isTRUE("Intercept" %in% mod[,1])) {
    params <- auto_cholesky(model = model, draws = draws)
    Intercept <- params[1]
    betas <- params[2:length(params)]
    lambda <- exp(Intercept + sum(features*betas))
    return(as.numeric(lambda^{-1}))
    
    #=====================================#
    # mu is not a function of covariates  #
    #=====================================#
    
  } else if (!isTRUE("beta1" %in% mod[,1])) {
    lambda <- -1 
    while(isTRUE(lambda<0)) {
      lambda <- auto_cholesky(model = model, draws = draws)
    }
    return(as.numeric(lambda^{-1}))
    
    #====================================================================#
    # linear predictor has no intercept but is a function of covariates  # 
    #====================================================================#
    
  } else {
    betas <- auto_cholesky(model = model, draws = draws)
    lambda <- exp(sum(betas*features))
    return(as.numeric(lambda^{-1}))
  }
}