#' Conditional mean for Generalized Gamma distribution
#' 
#' ggamma_flexpredict returns the conditional mean E(Y|X) of a model fitted via the function ggamma_flexfit; where μ has been specified to be a function of covariates the required value should be specified using the ‘features’ parameter. ggamma_flexpredict also allows for the correlation of estimated parameters via the Cholesky decomposition of the variance-covariance matrix.
#'@param model An object of class "mle2" produced using the function ggamma_flexfit.
#'@param features  A numeric vector specifying the value of covriates at which the conditional mean should be evaluated; the covariates in the vector should appear in the same order as they do in the model. Where a model does not depend on covariates the argument may be left blank. 
#'@param draws The number of random draws from multivariate random normal representing correlated parameters. If parameter correlation is not required draws should be set to zero.
#'@details This function uses the same parametrization of the Fisk distribution as is used in Stacy (1962); note that this is differs significantly to the parametrization used in many common R packages. The probability probability density function is used is:
#'@details f(y) = [α/Γ(φ)] • μ^αφ y^αφ-1 exp(-(μy)^α)
#'@details The function returns:
#'@details E(Y|X) = μΓ(φ+α^{-1})/Γ(φ)
#'@details μ may be a function of covariates; in which case, the cannonical log link function is used.
#'@references Kleiber, Christian, and Samuel Kotz. Statistical Size Distributions In Economics And Actuarial Sciences. pp. 107-147. John Wiley & Sons, 2003.

ggamma_flexpredict <- function(model, features, draws = 5) {
  
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
    mu <- exp(Intercept + sum(betas*features))
    return(as.numeric(mu*gamma(phi+alpha^{-1})/gamma(phi)))
    
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
    return(as.numeric(mu*gamma(phi+alpha^{-1})/gamma(phi)))
    
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
    mu <- exp(sum(features*betas))
    return(as.numeric(mu*gamma(phi+alpha^{-1})/gamma(phi)))
  }
}