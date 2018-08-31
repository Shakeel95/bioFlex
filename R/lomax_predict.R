#' Conditional mean for Lomax distribution
#' 
#' lomax_flexpredict returns the conditional mean E(Y|X) of a model fitted via the function lomax_flexfit; where λ has been specified to be a function of covariates the required value should be specified using the ‘features’ parameter. lomax_flexpredict also allows for the correlation of estimated parameters via the Cholesky decomposition of the variance-covariance matrix.
#'@param model An object of class "mle2" produced using the function lomax_flexfit.
#'@param features  A numeric vector specifying the value of covriates at which the conditional mean should be evaluated; the covariates in the vector should appear in the same order as they do in the model. Where a model does not depend on covariates the argument may be left blank. 
#'@param draws The number of random draws from multivariate random normal representing correlated parameters. If parameter correlation is not required draws should be set to zero.
#'@details This function uses the same parametrization of the Lomax distribution as is used in Kleiber and Kotz (2003). The probability probability density function is used is:
#'@details f(y) = (α/λ) [1 + (x/λ)]^-(α+1)
#'@details The function returns:
#'@details E(Y|X) = λ/(α-1)
#'@details λ may be a function of covariates; in which case, the cannonical log link function is used.
#'@references Kleiber, Christian, and Samuel Kotz. Statistical Size Distributions In Economics And Actuarial Sciences. pp. 107-147. John Wiley & Sons, 2003. Print.

lomax_flexpredict <- function(model, features, draws = 5) {
  
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
    return(as.numeric(lambda/(alpha -1)))
    
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
    return(as.numeric(lambda/(alpha -1)))
    
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
    return(as.numeric(lambda/(alpha -1)))
  }
}