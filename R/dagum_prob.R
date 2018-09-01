#' Conditional probability for Dagum distribution
#' 
#' dagum_flexprob returns the conditional probability P(Y<k | X) of a model fitted via the function dagum_flexfit; where b has been specified to be a function of covariates the required value should be specified using the ‘features’ parameter. The function includes a procedure for visualizing the conditional probability. dagum_flexprob also allows for the correlation of estimated parameters via the Cholesky decomposition of the variance-covariance matrix.
#'@param K Value for which P(Y<k | X) is computed.
#'@param model An object of class "mle2" produced using the function dagum_flexfit.
#'@param features A numeric vector specifying the value of covriates at which the conditional probability should be evaluated; the covariates in the vector should appear in the same order as they do in the model. Where a model does not depend on covariates the argument may be left blank.
#'@param visualise Logical. If TRUE (the default) the conditional distribution is plotted at P(Y<k | x) is shaded. 
#'@param xlim Numeric vectors of length 2, giving the coordinate range of the dependent variable.
#'@param draws The number of random draws from multivariate random normal representing correlated parameters. If parameter correlation is not required draws should be set to zero.
#'@details This function uses the same parametrization of the Dagum distribution as is used in Kleiber and Kotz (2003). The probability probability density function is used is:
#'@details f(y) = apy^ap-1/[b^ap(1+(y/b)^a)^p+1]
#'@details The function returns:
#'@details P(Y< k | X) = [1+(k/b)^{-a}]^{-p}
#'@details b may be a function of covariates; in which case, the cannonical log link function is used.
#'@references Kleiber, Christian, and Samuel Kotz. Statistical Size Distributions In Economics And Actuarial Sciences. pp. 107-147. John Wiley & Sons, 2003.
#'@export

dagum_flexprob <- function(K, model, features, visualise = TRUE, xlim, draws = 5) {
  
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
      return(as.numeric(pdagum(K, scale = b, shape1.a = a, shape2.p = p)))
    } else {
      return(as.numeric(pdagum(K, scale = b, shape1.a = a, shape2.p = p)))
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
      return(as.numeric(pdagum(K, scale = b, shape1.a = a, shape2.p = p)))
    } else {
      return(as.numeric(pdagum(K, scale = b, shape1.a = a, shape2.p = p)))
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
      return(as.numeric(pdagum(K, scale = b, shape1.a = a, shape2.p = p)))
    } else {
      return(as.numeric(pdagum(K, scale = b, shape1.a = a, shape2.p = p)))
    }
  }
}