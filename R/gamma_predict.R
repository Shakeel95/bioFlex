#' Conditional mean for Gamma distribution
#' 
#' gamma_flexpredict returns the conditional mean E(Y|X) of a model fitted via the function gamma_flexfit; where α has been specified to be a function of covariates the required value should be specified using the ‘features’ parameter. gamma_flexpredict also allows for the correlation of estimated parameters via the Cholesky decomposition of the variance-covariance matrix.
#'@param model An object of class "mle2" produced using the function gamma_flexfit.
#'@param features  A numeric vector specifying the value of covriates at which the conditional mean should be evaluated; the covariates in the vector should appear in the same order as they do in the model. Where a model does not depend on covariates the argument may be left blank. 
#'@param draws The number of random draws from multivariate random normal representing correlated parameters. If parameter correlation is not required draws should be set to zero.
#'@details This function uses the the most common parametrization of the Gamma distribution.The probability probability density function is used is:
#'@details f(y) = λ^α/Γ(α)•y^α-1 exp(-λy)
#'@details The function returns:
#'@details E(Y|X) = α/λ
#'@details α may be a function of covariates; in which case, the cannonical log link function is used.
#'@references Kempthorne. "Parameter Estimation Fitting Probability Distributions Method Of Moments." MIT 18.443 (2015)

gamma_flexpredict <- function(model, features, draws = 5) {
  
  mod <- as.data.frame(tidy(model))
  
  #================================================================#
  # linear predictor has intercept and is a function of covariates # 
  #================================================================#
  
  if (isTRUE("Intercept" %in% mod[,1])) {
    
    lambda <- -1 
    while(isTRUE(lambda<0)) {
      params <- auto_cholesky(model = model, draws = draws)
      lambda <- params[1]
      Intercept <- params[2]
      betas <- params[3:length(params)]
    }
    alpha <- exp(Intercept + sum(features*betas))
    return(as.numeric(alpha/lambda))
    
    #=====================================#
    # mu is not a function of covariates  #
    #=====================================# 
    
  } else if (!isTRUE("beta1" %in% mod[,1])) {
    
    alpha <- -1
    lambda <- -1 
    while(isTRUE((lambda<0)) | (alpha<0)){
      params <- auto_cholesky(model = model, draws = draws)
      alpha <- params[1]
      lambda <- params[2]
    }
    return(as.numeric(alpha/lambda))
    
    #====================================================================#
    # linear predictor has no intercept but is a function of covariates  # 
    #====================================================================#
    
  } else {
    
    lambda <- -1 
    while(isTRUE(lambda<0)) {
      params <- auto_cholesky(model = model, draws = draws)
      lambda <- params[1]
      betas <- params[2:length(params)]
    }
    alpha <- exp(sum(features*betas))
    return(as.numeric(alpha/lambda))
    
  }
}
  