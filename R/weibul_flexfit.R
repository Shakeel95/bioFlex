#' Fitting Weibull distribution via maximum likelihood
#'
#'weibull_flexfit is used to fit a Weibull distribution to a strictly positive response variable. The scale parameter may be specified either as a function of covariates or as a constant estimated using the response variable alone.
#'@param formula An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#'@param data An optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which weibull_flexfit is called.
#'@param weights An optional vector of weights to be used in the fitting process. Should be NULL or a numeric vector.
#'@param subset An optional vector specifying a subset of observations to be used in the fitting process.
#'@param ownstart An optional list containing starting values for the maximum likelihood estimation procedure. If a model with an intercept has been specified, the list must be of the form ownstart = list(kappa = , beta0 = , beta1 =, …); if a mode with no intercept has been specified the list must be of the form ownstart = list(kappa = , beta1 = , …); if the rate parameter is not a function of covariates the list must be of the form ownstart = list(lambda = , kappa = ). It is important that the list have as many elements as there are parameters in the model, and that these be supplied in the order set out above.
#'@param key A logical parameter dictating whether a key is produced alongside the model’s output.
#'@param warnings A logical parameter dictating whether warnings from the maximum likelihood estimation procedure are produced alongside the model’s output.
#'@param ... Additional arguments to be passed to the function optim within the maximum likelihood estimation procedure. Useful arguments include the gradient descent algorithm to be used and bounds on parameter values; see the stats package.
#'@details This function uses the same parametrization of the Weibull distribution as is used in Kleiber and Kotz (2003). Starting values for the maximum likelihood estimation procedure take from Niola et al. (2006), which uses a slightly different parametrization. The probability probability density function is used is: 
#'@details f(y) = (κ/λ) (y/λ)^{κ-1} exp(-y/λ)^{κ}
#'@details κ is the shape parameter and λ is the scale parameter. 
#'@details When the argument formula specifies a full model with an intercept, λ takes the following form and is estimated via a two step (least squares and maximum likelihood) procedure:
#'@details λ = β0 + β1x1 + …. + βkxk
#'@details When the formula argument specifies a model without an intercept, λ takes the bellow form and is estimated via a two step (least squares and maximum likelihood) procedure. Unless theory suggests that an intercept should not be used, users are advised to use a model with an intercept as the maximum likelihood estimation procedure is more stable.
#'@details λ = β1x1 + …. + βkxk
#'@details When a null model is specified (formula = y ~ 0) λ is not estimated as a function of covariates. The starting value for the maximum likelihood estimation procedure is obtained by calling the function weibull_scale.
#'@references Kleiber, Christian, and Samuel Kotz. Statistical Size Distributions In Economics And Actuarial Sciences. pp. 107-147. John Wiley & Sons, 2003. Print.
#'@references NIOLA, VINCENZO, ROSARIO OLIVIERO, and GIUSEPPE QUAREMBA. "A Method Of Moments For The Estimation Of Weibull Pdf Parameters." Proceedings of the 8th WSEAS Int. Conference on Automatic Control, Modeling and Simulation 8 (2014)
#'@export


weibull_flexfit <- function(formula, data = NULL, weights = NULL, subset = NULL, ownstart = NULL, key = TRUE, warnings = FALSE, ...) {

  #===================#
  # first stage code  #
  #===================#

  # construct data frame equivalen to formula
  data <- lm(formula = formula, data = data, method = "model.frame")
  # check if more than 30 covariates are used
  if(ncol(data) > 31)stop("At this time weibull_fit() accepts at most 30 covariates")


  # construct first stage regression
  dat <- data
  dat[,1] <- log(data[,1]/gamma(1+weibull_shape(data[,1])^{-1}))
  lm_start <- lm(formula = formula, data = dat)

  if (isTRUE("(Intercept)" %in% names(coefficients(lm_start)))) {

    #================================================================#
    # linear predictor has intercept and is a function of covariates #
    #================================================================#

    # build list of starting values for mle2
    if(is.null(ownstart)) {
      list_start <- list(kappa = weibull_shape(data[,1]), Intercept = coefficients(lm_start)[1])
      for (i in 2:ncol(data)) {
        list_start[i+1] = coefficients(lm_start)[i]
      }
      for (i in 2:length(coefficients(lm_start))) {
        names(list_start)[i+1] <- paste("beta", i-1, sep = "")
      }
    } else {
      list_start <- ownstart
    }

    # define log-likelihood
    LL <- function(kappa, Intercept, beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, beta6 = 0, beta7 = 0, beta8 = 0, beta9 = 0, beta10 = 0, beta11 = 0, beta12 = 0, beta13 = 0, beta14 = 0, beta15 = 0, beta16 = 0, beta17 = 0, beta18 = 0, beta19 = 0, beta20 = 0, beta21 = 0, beta22 = 0, beta23 = 0, beta24 = 0, beta25 = 0, beta26 = 0, beta27 = 0, beta28 = 0, beta29 = 0, beta30 = 0) {
      mu = 0
      for (i in 2:ncol(data)) {
        mu = mu + data[,i]*get(paste0("beta",i-1))
      }
      mu <- mu + Intercept
      R = dweibull(data[,1], shape = kappa, scale = exp(mu), log = TRUE)
      -sum(R)
    }

    # maximise log-likelihood
    if(!isTRUE(warnings)) {
      MLE <- suppressWarnings(mle2(LL, start = list_start, ...))
    } else {
      MLE <- mle2(LL, start = list_start, ...)
    }

    # generate key and return
    if (isTRUE(key)) {
      Parameters <- c("kappa", "Intercept")
      Equals <- c("->","->")
      Interpretation <- c("shape parameter", "intercept in linear predictor")
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

    #=========================================#
    # lambda is not a function of covariates  #
    #=========================================#

  } else if (is.null(names(coefficients(lm_start)))) {

    # build list of starting values for mle2
    if (is.null(ownstart)) {
      list_start <- list(lambda = weibull_scale(data[,1]), kappa = weibull_shape(data[,1]))
    } else {
      list_start <- ownstart
    }

    # define log-likelihood
    LL <- function(lambda, kappa) {
      R = dweibull(data[,1], shape = kappa, scale = lambda, log = TRUE)
      -sum(R)
    }

    # maximise log-likelihood
    if(!isTRUE(warnings)) {
      MLE <- suppressWarnings(mle2(LL, start = list_start, ...))
    } else {
      MLE <- mle2(LL, start = list_start, ...)
    }

    # generate key and return
    if (isTRUE(key)) {
      Parameters <- c("lambda", "kappa")
      Equals <- c("->","->")
      Interpretation <- c("scale parameter", "shape parameter")
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
      list_start <- list(kappa = weibull_shape(data[,1]))
      for (i in 1:length(coefficients(lm_start))) {
        list_start[i+1] <- coefficients(lm_start)[i]
      }
      for (i in 1:length(coefficients(lm_start))) {
        names(list_start)[i+1] <- paste("beta",i, sep = "")
      }
    } else {
      list_start <- ownstart
    }

    # define log-likelihood
    LL <- function(kappa, beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, beta6 = 0, beta7 = 0, beta8 = 0, beta9 = 0, beta10 = 0, beta11 = 0, beta12 = 0, beta13 = 0, beta14 = 0, beta15 = 0, beta16 = 0, beta17 = 0, beta18 = 0, beta19 = 0, beta20 = 0, beta21 = 0, beta22 = 0, beta23 = 0, beta24 = 0, beta25 = 0, beta26 = 0, beta27 = 0, beta28 = 0, beta29 = 0, beta30 = 0) {
      mu = 0
      for (i in 2:ncol(data)) {
        mu = mu + data[,i]*get(paste0("beta",i-1))
      }
      R = dweibull(data[,1], shape = kappa, scale = exp(mu), log = TRUE)
      -sum(R)
    }

    # maximise log-likelihood
    if(!isTRUE(warnings)) {
      MLE <- suppressWarnings(mle2(LL, start = list_start, ...))
    } else {
      MLE <- mle2(LL, start = list_start, ...)
    }

    # generate key and return
    if (isTRUE(key)) {
      Parameters <- c("kappa")
      Equals <- c("->")
      Interpretation <- c("shape parameter")
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
