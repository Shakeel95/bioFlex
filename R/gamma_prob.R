#' Conditional probability for Gamma distribution
#' 
#'gamma_flexprob returns the conditional probability P(Y<k | X) of a model fitted via the function gamma_flexfit; where λ has been specified to be a function of covariates the required value should be specified using the ‘features’ parameter. The function includes a procedure for visualizing the conditional probability. gamma_flexprob also allows for the correlation of estimated parameters via the Cholesky decomposition of the variance-covariance matrix.
#'@param K Value for which P(Y<k | X) is computed.
#'@param model An object of class "mle2" produced using the function gamma_flexfit.
#'@param features A numeric vector specifying the value of covriates at which the conditional probability should be evaluated; the covariates in the vector should appear in the same order as they do in the model. Where a model does not depend on covariates the argument may be left blank.
#'@param visualise Logical. If TRUE (the default) the conditional distribution is plotted at P(Y<k | x) is shaded. 
#'@param xlim Numeric vectors of length 2, giving the coordinate range of the dependent variable.
#'@param draws The number of random draws from multivariate random normal representing correlated parameters. If parameter correlation is not required draws should be set to zero.
#'@details This function uses the the most common parametrization of the Gamma distribution.The probability probability density function is used is:
#'@details f(y) = λ^α/Γ(α)•y^α-1 exp(-λy)
#'@details The function returns:
#'@details P(Y<k | X) = (1/Γ(α))γ(α,kλ)
#'@details α may be a function of covariates; in which case, the cannonical log link function is used.
#'@references Kempthorne. "Parameter Estimation Fitting Probability Distributions Method Of Moments." MIT 18.443 (2015)
#'@references R. V. Hogg and A. T. Craig (1978) Introduction to Mathematical Statistics, 4th edition. New York: Macmillan. (See Section 3.3.)
#'@export 

gamma_flexprob <- function(K, model, features, visualise = TRUE, xlim, draws = 5) {
  
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
    
    if(isTRUE(visualise)) {
      preview <<- function(x) {
        dgamma(x, shape = alpha, rate = lambda)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(as.numeric(pgamma(K, shape = alpha, rate = lambda)))
    } else {
      return(as.numeric(pgamma(K, shape = alpha, rate = lambda)))
    }
    
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
    
    if(isTRUE(visualise)) {
      preview <<- function(x) {
        dgamma(x, shape = alpha, rate = lambda)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(as.numeric(pgamma(K, shape = alpha, rate = lambda)))
    } else {
      return(as.numeric(pgamma(K, shape = alpha, rate = lambda)))
    }
    
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
    
    if(isTRUE(visualise)) {
      preview <<- function(x) {
        dgamma(x, shape = alpha, rate = lambda)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(as.numeric(pgamma(K, shape = alpha, rate = lambda)))
    } else {
      return(as.numeric(pgamma(K, shape = alpha, rate = lambda)))
    }
  }
}
  