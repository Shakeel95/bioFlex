#' Conditional mean for Singh-Maddala distribution
#' 
#' sinmad_flexpredict returns the conditional mean E(Y|X) of a model fitted via the function sinmad_flexfit; where b has been specified to be a function of covariates the required value should be specified using the ‘features’ parameter. sinmad_flexpredict also allows for the correlation of estimated parameters via the Cholesky decomposition of the variance-covariance matrix.
#'@param model An object of class "mle2" produced using the function sinmad_flexfit.
#'@param features  A numeric vector specifying the value of covriates at which the conditional mean should be evaluated; the covariates in the vector should appear in the same order as they do in the model. Where a model does not depend on covariates the argument may be left blank. 
#'@param draws The number of random draws from multivariate random normal representing correlated parameters. If parameter correlation is not required draws should be set to zero.
#'@details This function uses the same parametrization of the Singh-Maddala distribution as is used in Kleiber and Kotz (2003).  The probability probability density function is used is:
#'@details f(y) = aqy^a-1/[b^a(1+(y/b)^a)^1+q]
#'@details The function returns:
#'@details E(Y|X) = bΓ(1+1/a)Γ(q-1/a)/Γ(q)
#'@details b may be a function of covariates; in which case, the cannonical log link function is used.
#'@references Kleiber, Christian, and Samuel Kotz. Statistical Size Distributions In Economics And Actuarial Sciences. pp. 107-147. John Wiley & Sons, 2003.

sinmad_flexpredict <- function(model, features, draws = 5) {
  
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
    return(as.numeric((b*gamma(1+a^{-1})*gamma(q-a^{-1}))/gamma(q)))
    
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
    return(as.numeric((b*gamma(1+a^{-1})*gamma(q-a^{-1}))/gamma(q)))
    
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
    return(as.numeric((b*gamma(1+a^{-1})*gamma(q-a^{-1}))/gamma(q)))
  }
}