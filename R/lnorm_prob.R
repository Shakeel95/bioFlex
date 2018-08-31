#'Conditional probability for Log-Normal distribution
#'
#'lnorm_flexprob returns the conditional probability P(Y<k | X) of a model fitted via the function lnorm_flexfit; where μ has been specified to be a function of covariates the required value should be specified using the ‘features’ parameter. The function includes a procedure for visualizing the conditional probability. lnorm_flexprob also allows for the correlation of estimated parameters via the Cholesky decomposition of the variance-covariance matrix.
#'@param K Value for which P(Y<k | X) is computed.
#'@param model An object of class "mle2" produced using the function lnorm_flexfit.
#'@param features A numeric vector specifying the value of covriates at which the conditional probability should be evaluated; the covariates in the vector should appear in the same order as they do in the model. Where a model does not depend on covariates the argument may be left blank.
#'@param visualise Logical. If TRUE (the default) the conditional distribution is plotted at P(Y<k | x) is shaded. 
#'@param xlim Numeric vectors of length 2, giving the coordinate range of the dependent variable.
#'@param draws The number of random draws from multivariate random normal representing correlated parameters. If parameter correlation is not required draws should be set to zero.
#'@details This function uses the two parameter parametrization of the Log-Normal distribution. The probability probability density function is used is:
#'@details f(y) = [yσ(2π)^1/2]^-1 exp(-log(y-μ)^2/(2σ^2))
#'@details The function returns:
#'@details P(Y<k | X) = 0.5+0.5erf[log(k-μ)/(√2σ)]
#'@details μ may be a function of covariates; in which case, the identity link function is used.
#'@references Faith Ginos, Brenda. "Parameter Estimation For The Lognormal Distribution." Brigham Young University Scholars Archive (2018): 1-111. Web. 10 Aug. 2018.
#'@references Kleiber, Christian, and Samuel Kotz. Statistical Size Distributions In Economics And Actuarial Sciences. pp. 107-147. John Wiley & Sons, 2003.
#'@export 

lnorm_flexprob <- function(K, model, features, visualise = TRUE, xlim, draws = 5) {
  
  mod <- as.data.frame(tidy(model))
  
  #================================================================#
  # linear predictor has intercept and is a function of covariates # 
  #================================================================#
  
  if (isTRUE("Intercept" %in% mod[,1])) {
    
    sigma <- -1 
    while(isTRUE(sigma <0)) {
      params <- auto_cholesky(model = model, draws = draws)
      sigma <- params[1]
      Intercept <- params[2]
      betas <- params[3:length(params)]
    }
    mu <- Intercept + sum(betas*features)
    
    if (isTRUE(visualise)) {
      preview <<- function(x) {
        dlnorm(x, meanlog = mu, sdlog = sigma)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(plnorm(K, meanlog = mu, sdlog = sigma))
    } else {
      return(plnorm(K, meanlog = mu, sdlog = sigma))
    }
    
    #=====================================#
    # mu is not a function of covariates  #
    #=====================================# 
    
  } else if (!isTRUE("beta1" %in% mod[,1])) {
    
    sigma <- -1
    while(isTRUE(sigma < 0)) {
      params <- auto_cholesky(model = model, draws = draws)
      sigma <- params[1]
      mu <- params[2]
    }
    
    if (isTRUE(visualise)) {
      preview <<- function(x) {
        dlnorm(x, meanlog = mu, sdlog = sigma)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(plnorm(K, meanlog = mu, sdlog = sigma))
    } else {
      return(plnorm(K, meanlog = mu, sdlog = sigma))
    }
    
    #====================================================================#
    # linear predictor has no intercept but is a function of covariates  # 
    #====================================================================#
    
  } else {
    
    sigma <- -1 
    while(isTRUE(sigma<0)) {
      params <- auto_cholesky(model = model, draws = draws)
      sigma <- params[1]
      betas <- params[2:length(params)]
    }
    mu <- sum(betas*features)
    
    if (isTRUE(visualise)) {
      preview <<- function(x) {
        dlnorm(x, meanlog = mu, sdlog = sigma)
      }
      plot(preview, xlim = xlim, ylab = "Density", xlab = "", lwd = 3)
      Shade(preview, breaks = c(0,K), xlim = xlim)
      abline(a = 0, b = 0)
      return(plnorm(K, meanlog = mu, sdlog = sigma))
    } else {
      return(plnorm(K, meanlog = mu, sdlog = sigma))
    }
    
  }
}