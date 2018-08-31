#' Conditional mean for Fisk distribution
#' 
#' fisk_flexpredict returns the conditional mean E(Y|X) of a model fitted via the function fisk_flexfit; where α has been specified to be a function of covariates the required value should be specified using the ‘features’ parameter. fisk_flexpredict also allows for the correlation of estimated parameters via the Cholesky decomposition of the variance-covariance matrix.
#'@param model An object of class "mle2" produced using the function fisk_flexfit.
#'@param features  A numeric vector specifying the value of covriates at which the conditional mean should be evaluated; the covariates in the vector should appear in the same order as they do in the model. Where a model does not depend on covariates the argument may be left blank. 
#'@param draws The number of random draws from multivariate random normal representing correlated parameters. If parameter correlation is not required draws should be set to zero.
#'@details This function uses the same parametrization of the Fisk distribution as is used in Kleiber and Kotz (2003). The probability probability density function is used is:
#'@details f(y) = αy^α-1/[β^α(1+(y/β)^α]^2
#'@details The function returns:
#'@details E(Y|X) = απ/sin(π/β)
#'@details α may be a function of covariates; in which case, the cannonical log link function is used.
#'@references Kleiber, Christian, and Samuel Kotz. Statistical Size Distributions In Economics And Actuarial Sciences. pp. 107-147. John Wiley & Sons, 2003.

fisk_flexpredict <- function(model, features, draws = 5) {
  
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
    return(as.numeric((alpha*pi/beta)/(sin(pi/beta))))
    
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
    return(as.numeric((alpha*pi/beta)/(sin(pi/beta))))
    
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
    return(as.numeric((alpha*pi/beta)/(sin(pi/beta))))
  }
}