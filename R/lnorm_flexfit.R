#' Fitting Log-Normal distribution via maximum likelihood
#'
#'lnorm_flexfit is used to fit a Log-Normal distribution to a strictly positive response variable. The parameter mu may be specified either as a function of covariates or as a constant estimated using the response variable alone.
#'If mu is specified to be a function of covariates, the identity link function is used.
#'@param formula An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#'@param data An optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lnorm_flexfit is called.
#'@param weights An optional vector of weights to be used in the fitting process. Should be NULL or a numeric vector.
#'@param subset An optional vector specifying a subset of observations to be used in the fitting process.
#'@param ownstart An optional list containing starting values for the maximum likelihood estimation procedure. If a model with an intercept has been specified, the list must be of the form ownstart = list(sigma = , beta0 = , beta1 =, …); if a mode with no intercept has been specified the list must be of the form ownstart = list(sigma = , beta1 = , …); if the rate parameter is not a function of covariates the list must be of the form ownstart = list(sigma = , mu = ). It is important that the list have as many elements as there are parameters in the model, and that these be supplied in the order set out above.
#'@param key A logical parameter dictating whether a key is produced alongside the model’s output.
#'@param warnings A logical parameter dictating whether warnings from the maximum likelihood estimation procedure are produced alongside the model’s output.
#'@param ... Additional arguments to be passed to the function optim within the maximum likelihood estimation procedure. Useful arguments include the gradient descent algorithm to be used and bounds on parameter values; see the stats package.
#'@details This function uses the two parameter parametrization of the Log-Normal distribution.Starting values for the maximum likelihood estimation procedure are taken from Brenda Faith Ginos’ Masters dissertation on the Log-Normal distribution. The probability probability density function is used is:
#'@details f(y) = [yσ(2π)^{1/2}]^{-1} exp(-log(y-μ)^{2}/(2σ^{2}))
#'@details The parameters μ and σ are not strictly location and scale parameters. Rather, they are the scale and location parameters for the normally distributed logarithm log(y), where y follows a Log-Normal distribution. In the first stage estimation functions  μ and σ are referred to as ‘effective’ scale and location parameters. 
#'@details When the argument formula specifies a full model with an intercept, μ takes the following form and is estimated via a two step (least squares and maximum likelihood) procedure:
#'@details μ = β0 + β1x1 + …. + βkxk
#'@details When the formula argument specifies a model without an intercept, μ takes the bellow form and is estimated via a two step (least squares and maximum likelihood) procedure. Unless theory suggests that an intercept should not be used, users are advised to use a model with an intercept as the maximum likelihood estimation procedure is more stable.
#'@details μ = β1x1 + …. + βkxk
#'@details When a null model is specified (formula = y ~ 0) μ is not estimated as a function of covariates. The starting value for the maximum likelihood estimation procedure is obtained by calling the function lnorm_mean. 
#'@references Faith Ginos, Brenda. "Parameter Estimation For The Lognormal Distribution." Brigham Young University Scholars Archive (2018): 1-111. Web. 10 Aug. 2018.
#'@references  Kleiber, Christian, and Samuel Kotz. Statistical Size Distributions In Economics And Actuarial Sciences. pp. 107-147. John Wiley & Sons, 2003. Print.
#'@export


lnorm_flexfit <- function(formula, data = NULL, weights = NULL, subset = NULL, ownstart = NULL, key = TRUE, warnings = FALSE, ...) {

  #===================#
  # first stage code  #
  #===================#

  # construct data frame equivalen to formula
  data <- lm(formula = formula, data = data, method = "model.frame")
  # check if more than 30 covariates are used
  if(ncol(data) > 31)stop("At this time lnorm_fit() accepts at most 30 covariates")

  # construct first stage regression
  dat <- data
  dat[,1] <- log(dat[,1]) - lnorm_var(dat[,1])/2
  lm_start <- lm(formula = formula, data = dat)

  if (isTRUE("(Intercept)" %in% names(coefficients(lm_start)))) {

    #================================================================#
    # linear predictor has intercept and is a function of covariates #
    #================================================================#

    # build list of starting values for mle2
    if (is.null(ownstart)) {
      list_start <- list(sigma = sqrt(lnorm_var(data[,1])), Intercept = coefficients(lm_start)[1])
      for (i in 2:length(coefficients(lm_start))) {
        list_start[i+1] = coefficients(lm_start)[i]
      }
      for (i in 2:length(coefficients(lm_start))) {
        names(list_start)[i+1] <- paste("beta", i-1, sep = "")
      }
    } else {
      list_start <- ownstart
    }

    # define log-likelihood function
    LL <- function(sigma, Intercept, beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, beta6 = 0, beta7 = 0, beta8 = 0, beta9 = 0, beta10 = 0, beta11 = 0, beta12 = 0, beta13 = 0, beta14 = 0, beta15 = 0, beta16 = 0, beta17 = 0, beta18 = 0, beta19 = 0, beta20 = 0, beta21 = 0, beta22 = 0, beta23 = 0, beta24 = 0, beta25 = 0, beta26 = 0, beta27 = 0, beta28 = 0, beta29 = 0, beta30 = 0) {
      mu = 0
      for (i in 2:ncol(data)) {
        mu = mu + data[,i]*get(paste0("beta",i-1))
      }
      mu = Intercept + mu
      R = dlnorm(data[,1], meanlog = mu, sdlog = sigma, log = TRUE)
      -sum(R)
    }

    # maximise log-likelihood
    if(!isTRUE(warnings)) {
      MLE <- suppressWarnings(mle2(LL, start = list_start, ...))
    } else {
      MLE <- mle2(LL, start = list_start, ...)
    }

    # generate key and return
    if(isTRUE(key)) {
      Parameters <- c("sigma", "Intercept")
      Equals <- c("->","->")
      Interpretation <- c("sd of response variable's logarithm", "intercept in linear predictor")
      for (i in 2:length(coefficients(lm_start))) {
        Parameters <- c(Parameters, names(list_start)[i+1])
        Equals <- c(Equals, "->")
        Interpretation <- c(Interpretation, paste(names(coefficients(lm_start))[i], " coefficient in linear predictor"))
      }
      key <- data.frame(Parameters, Equals, Interpretation)
      print("Key:")
      print(key)
      return(MLE)
    } else {
      return(MLE)
    }

    #=====================================#
    # mu is not a function of covariates  #
    #=====================================#

  } else if (is.null(names(coefficients(lm_start)))) {

    # build list of starting values for mle2
    if (is.null(ownstart)) {
      list_start <- list(sigma = sqrt(lnorm_var(data[,1])), mu = lnorm_mean(data[,1]))
    } else {
      list_start <- ownstart
    }

    # define log-likelihood function
    LL <- function(sigma, mu) {
      R = dlnorm(data[,1], meanlog = mu, sdlog = sigma, log = TRUE)
      -sum(R)
    }

    # maximise log-likelihood
    if(!isTRUE(warnings)) {
      MLE <- suppressWarnings(mle2(LL, start = list_start, ...))
    } else {
      MLE <- mle2(LL, start = list_start, ...)
    }

    # generate key and return
    if(isTRUE(key)) {
      Parameters <- c("mu","sigma")
      Equals <- c("->", "->")
      Interpretation <- c("mean of the response variable's logarithm","sd of response variable's logarithm")
      key <- data.frame(Parameters, Equals, Interpretation)
      print("Key:")
      print(key)
      return(MLE)
    } else {
      return(MLE)
    }

    #====================================================================#
    # linear predictor has no intercept but is a function of covariates  #
    #====================================================================#

  } else {

    # build list of starting values for mle2
    if (is.null(ownstart)) {
      list_start <- list(sigma = sqrt(lnorm_var(data[,1])))
      for (i in 1:length(coefficients(lm_start))) {
        list_start[i+1] = coefficients(lm_start)[i]
      }
      for (i in 1:length(coefficients(lm_start))) {
        names(list_start)[i+1] <- paste("beta", i, sep = "")
      }
    } else {
      list_start <- ownstart
    }


    # define log-likelihood
    LL <- function(sigma, beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, beta6 = 0, beta7 = 0, beta8 = 0, beta9 = 0, beta10 = 0, beta11 = 0, beta12 = 0, beta13 = 0, beta14 = 0, beta15 = 0, beta16 = 0, beta17 = 0, beta18 = 0, beta19 = 0, beta20 = 0, beta21 = 0, beta22 = 0, beta23 = 0, beta24 = 0, beta25 = 0, beta26 = 0, beta27 = 0, beta28 = 0, beta29 = 0, beta30 = 0) {
      mu = 0
      for (i in 2:ncol(data)) {
        mu = mu + data[,i]*get(paste0("beta",i-1))
      }
      R = dlnorm(data[,1], meanlog = mu, sdlog = sigma, log = TRUE)
      -sum(R)
    }

    # maximise log-likelihood
    if(!isTRUE(warnings)) {
      MLE <- suppressWarnings(mle2(LL, start = list_start, ...))
    } else {
      MLE <- mle2(LL, start = list_start, ...)
    }

    # generate key and return
    if(isTRUE(key)) {
      Parameters <- c("sigma")
      Equals <- c("->")
      Interpretation <- c("sd of response variable's logarithm")
      for (i in 1:length(coefficients(lm_start))) {
        Parameters <- c(Parameters, names(list_start)[i+1])
        Equals <- c(Equals, "->")
        Interpretation <- c(Interpretation, paste(names(coefficients(lm_start))[i], " coefficient in linear predictor"))
      }
      key <- data.frame(Parameters, Equals, Interpretation)
      print("Key:")
      print(key)
      return(MLE)
    } else {
      return(MLE)
    }
  }
}
