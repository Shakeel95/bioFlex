#' Conditional probability for Exponential distribution
#' 
#'exp_flexprob returns the conditional probability P(Y<k | X) of a model fitted via the function exp_flexfit; where λ has been specified to be a function of covariates the required value should be specified using the ‘features’ parameter. The function includes a procedure for visualizing the conditional probability. exp_flexprob also allows for the correlation of estimated parameters via the Cholesky decomposition of the variance-covariance matrix.
#'@param K Value for which P(Y<k | X) is computed.
#'@param model An object of class "mle2" produced using the function exp_flexfit.
#'@param features A numeric vector specifying the value of covriates at which the conditional probability should be evaluated; the covariates in the vector should appear in the same order as they do in the model. Where a model does not depend on covariates the argument may be left blank.
#'@param visualise Logical. If TRUE (the default) the conditional distribution is plotted at P(Y<k | x) is shaded. 
#'@param xlim Numeric vectors of length 2, giving the coordinate range of the dependent variable.
#'@param draws The number of random draws from multivariate random normal representing correlated parameters. If parameter correlation is not required draws should be set to zero.
#'@details This function uses the most common parametrization of the Exponential distribution. The probability probability density function is used is:
#'@details f(y) = λexp(-λy)
#'@details The function returns: 
#'@details P(Y<k | X) = 1-exp(-kλ)
#'@details λ may be a function of covariates; in which case, the canonical log link function is used.
#'@export

exp_flexprob <- function(K, model, features, visualise = TRUE, xlim, draws = 5) {
  
  mod <- as.data.frame(tidy(model))
  
  #================================================================#
  # linear predictor has intercept and is a function of covariates # 
  #================================================================#
  
  if (isTRUE("Intercept" %in% mod[,1])) {
    
    params <- auto_cholesky(model = model, draws = draws)
    Intercept <- params[1]
    betas <- params[2:length(params)]
    lambda <- exp(Intercept + sum(features*betas))
    
    if(isTRUE(visualise)) {
      preview <<- function(x) {
        dexp(x, rate = lambda)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(as.numeric(pexp(K, rate = lambda)))
    } else {
      return(as.numeric(pexp(K, rate = lambda)))
    }
    
    #=====================================#
    # mu is not a function of covariates  #
    #=====================================#
    
  } else if (!isTRUE("beta1" %in% mod[,1])) {
    
    lambda <- -1 
    while(isTRUE(lambda<0)) {
      lambda <- auto_cholesky(model = model, draws = draws)
    }
    
    if(isTRUE(visualise)) {
      preview <<- function(x) {
        dexp(x, rate = lambda)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(as.numeric(pexp(K, rate = lambda)))
    } else {
      return(as.numeric(pexp(K, rate = lambda)))
    }
    
    #====================================================================#
    # linear predictor has no intercept but is a function of covariates  # 
    #====================================================================#
    
  } else {
    
    betas <- auto_cholesky(model = model, draws = draws)
    lambda <- exp(sum(betas*features))
    
    if(isTRUE(visualise)) {
      preview <<- function(x) {
        dexp(x, rate = lambda)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(as.numeric(pexp(K, rate = lambda)))
    } else {
      return(as.numeric(pexp(K, rate = lambda)))
    }
  }
}
  