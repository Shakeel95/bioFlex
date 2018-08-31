#' Fitting Generalized Beta of the Second Kind via maximum likelihood
#'
#'gb2_flexfit is used to fit a Generalized Beta distribution of the Second Kind distribution to a strictly positive response variable. The scale parameter may be specified either as a function of covariates or as a constant estimated using the response variable alone.
#'If the scale parameter is specified to be a function of covariates, the canonical log link function is used.
#'@param formula An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#'@param data An optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which gamma_flexfit is called.
#'@param weights An optional vector of weights to be used in the fitting process. Should be NULL or a numeric vector.
#'@param subset An optional vector specifying a subset of observations to be used in the fitting process.
#'@param ownstart An optional list containing starting values for the maximum likelihood estimation procedure. If a model with an intercept has been specified, the list must be of the form ownstart = list(a = , p = , q = , beta0 = , beta1 =, …); if a mode with no intercept has been specified the list must be of the form ownstart = list(a = , p = , q = , beta1 = , …); if the rate parameter is not a function of covariates the list must be of the form ownstart = list(b = , a = , p = , q =). It is important that the list have as many elements as there are parameters in the model, and that these be supplied in the order set out above.
#'@param key A logical parameter dictating whether a key is produced alongside the model’s output.
#'@param warnings A logical parameter dictating whether warnings from the maximum likelihood estimation procedure are produced alongside the model’s output.
#'@param ... Additional arguments to be passed to the function optim within the maximum likelihood estimation procedure. Useful arguments include the gradient descent algorithm to be used and bounds on parameter values; see the stats package.
#'@details This function uses the same parametrization of the Generalized Beta distribution of the Second Kind (GB2) as is used in Kleiber and Kotz (2003). Starting values for the maximum likelihood estimation procedure are obtained by maximizing the likelihood function of the GB2 with starting values p = 1 and q = 1, and the remaining starting values taken from the method of moments estimates of the parameters in the Fisk distribution. The probability probability density function is used is:
#'@details f(y) = ay^{ap-1}/[b^{ap}Β(p,q)(1+(y/b)^{a})^{p+q}]
#'@details When the argument formula specifies a full model with an intercept, b takes the following form and is estimated via a three step (maximum likelihood, least squares, maximum likelihood) procedure:
#'@details b = exp(β0 + β1x1 + …. + βkxk)
#'@details When the formula argument specifies a model without an intercept, b takes the bellow form. Estimating a model without an intercept is not advised as the estimates are extremely unstable; in simulations and test with medical data, the maximum likelihood procedure often often did not converge. Users wanting to estimate a model without an intercept are advised to supply their own starting values. 
#'@details b = exp(β1x1 + …. + βkxk)
#'@details When a null model is specified (formula = y ~ 0) b is not estimated as a function of covariates. The starting value for the maximum likelihood estimation procedure is obtained by calling the function weibull_scale.
#'@references Kleiber, Christian, and Samuel Kotz. Statistical Size Distributions In Economics And Actuarial Sciences. John Wiley & Sons, 2003. Print.
#'@export

gb2_flexfit <- function(formula, data = NULL, weights = NULL, subset = NULL, ownstart = NULL, key = TRUE, warnings = FALSE, warnings1 = FALSE, ...) {

  #===================#
  # first stage code  #
  #===================#

  # construct data frame equivalen to formula
  data <- lm(formula = formula, data = data, method = "model.frame")
  # check if more than 30 covariates are used
  if(ncol(data) > 31)stop("At this time weibull_fit() accepts at most 30 covariates")

  #=========================================================================#
  # perform first stage maximum likelihood estimation to obtain parameters  #
  #=========================================================================#

  # build list of starting values for mle2
  list_start1 <- list(a = fisk_shape(data[,1]), b = fisk_scale(data[,1]), p = 1, q = 1)

  # define log-likelihood
  LL1 <- function(b, a, p, q) {
    R = dgb2(data[,1], scale = b, shape1 = a, shape2 = p, shape3 = q)
    -sum(log(R))
  }

  # maximise log-likelihood and save parameter values
  if (!isTRUE(warnings1)) {
    MLE1 <- suppressWarnings(mle2(LL1, start = list_start1))
  } else {
    MLE1 <- mle2(LL1, start = list_start1)
  }
  b1 <- as.numeric(tidy(MLE1)[1,2])
  a1 <- as.numeric(tidy(MLE1)[2,2])
  p1 <- as.numeric(tidy(MLE1)[3,2])
  q1 <- as.numeric(tidy(MLE1)[4,2])

  # construct first stage regression
  dat <- data
  cons <- (gamma(p1 + a1^{-1})*gamma(q1 - a1^{-1}))/(gamma(p1)*gamma(q1))
  dat[,1] <- log(dat[,1]/cons)
  lm_start <- lm(formula = formula, data = dat)

  if (isTRUE("(Intercept)" %in% names(coefficients(lm_start)))) {

    #================================================================#
    # linear predictor has intercept and is a function of covariates #
    #================================================================#

    # build list of starting values for mle2
    if (is.null(ownstart)) {
      list_start <- list(a = a1, p = p1, q = q1, Intercept = coefficients(lm_start)[1])
      for (i in 2:ncol(data)) {
        list_start[i+3] = coefficients(lm_start)[i]
      }
      for (i in 2:ncol(data)) {
        names(list_start)[i+3] <- paste("beta", i-1, sep = "")
      }
    } else {
      list_start <- ownstart
    }

    # define log-likelihood
    LL <- function(a, p, q, Intercept, beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, beta6 = 0, beta7 = 0, beta8 = 0, beta9 = 0, beta10 = 0, beta11 = 0, beta12 = 0, beta13 = 0, beta14 = 0, beta15 = 0, beta16 = 0, beta17 = 0, beta18 = 0, beta19 = 0, beta20 = 0, beta21 = 0, beta22 = 0, beta23 = 0, beta24 = 0, beta25 = 0, beta26 = 0, beta27 = 0, beta28 = 0, beta29 = 0, beta30 = 0) {
      mu = 0
      for (i in 2:ncol(data)) {
        mu = mu + data[,i]*get(paste0("beta",i-1))
      }
      mu <- mu + Intercept
      R = dgb2(data[,1], scale = exp(mu), shape1 = a, shape2 = p, shape3 = q)
      -sum(log(R))
    }

    # maximise log-likelihood
    if(!isTRUE(warnings)) {
      MLE <- suppressWarnings(mle2(LL, start = list_start, ...))
    } else {
      MLE <- mle2(LL, start = list_start, ...)
    }

    # generate key and return
    if (isTRUE(key)) {
      Parameters <- c("a", "p", "q", "Intercept")
      Equals <- c("->","->", "->","->")
      Interpretation <- c("first shape parameter; controls heaviness of tails", "second shape parameter; controls skewness along with q", "third shape parameter; controls skewness along with p", "intercept in linear predictor")
      for (i in 2:length(coefficients(lm_start))) {
        Parameters <- c(Parameters, names(list_start)[i+3])
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

    #====================================#
    # b is not a function of covariates  #
    #====================================#

  } else if (is.null(names(coefficients(lm_start)))) {

    # build list of starting values for mle2
    if(is.null(ownstart)) {
      list_start <- list(b = b1, a = a1, p = p1, q = q1)
    } else {
      list_start <- ownstart
    }

    # define log likelihood
    LL <- function(b, a, p, q) {
      R = dgb2(data[,1], scale = b, shape1 = a, shape2 = p, shape3 = q)
      -sum(log(R))
    }

    # maximise log-likelihood
    if(!isTRUE(warnings)) {
      MLE <- suppressWarnings(mle2(LL, start = list_start, ...))
    } else {
      MLE <- mle2(LL, start = list_start, ...)
    }

    # generate key and return
    if (isTRUE(key)) {
      Parameters <- c("a", "b" ,"p", "q")
      Equals <- c("->","->", "->","->")
      Interpretation <- c("first shape parameter; controls heaviness of tails", "scale parameter" ,"second shape parameter; controls skewness along with q", "third shape parameter; controls skewness along with p")
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
      list_start <- list(a = a1, p = p1, q = q1)
      for (i in 1:length(coefficients(lm_start))) {
        list_start[i+3] = coefficients(lm_start)[i]
      }
      for (i in 1:length(coefficients(lm_start))) {
        names(list_start)[i+3] <- paste("beta", i, sep = "")
      }
    } else {
      list_start <- ownstart
    }

    # define log-likelihood
    LL <- function(a, p, q, beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, beta6 = 0, beta7 = 0, beta8 = 0, beta9 = 0, beta10 = 0, beta11 = 0, beta12 = 0, beta13 = 0, beta14 = 0, beta15 = 0, beta16 = 0, beta17 = 0, beta18 = 0, beta19 = 0, beta20 = 0, beta21 = 0, beta22 = 0, beta23 = 0, beta24 = 0, beta25 = 0, beta26 = 0, beta27 = 0, beta28 = 0, beta29 = 0, beta30 = 0) {
      mu = 0
      for (i in 2:ncol(data)) {
        mu = mu + data[,i]*get(paste0("beta",i-1))
      }
      R = dgb2(data[,1], scale = exp(mu), shape1 = a, shape2 = p, shape3 = q)
      -sum(log(R))
    }

    # maximise log-likelihood
    if(!isTRUE(warnings)) {
      MLE <- suppressWarnings(mle2(LL, start = list_start, ...))
    } else {
      MLE <- mle2(LL, start = list_start, ...)
    }

    # generate key and return
    if (isTRUE(key)) {
      Parameters <- c("a", "p", "q")
      Equals <- c("->","->", "->")
      Interpretation <- c("first shape parameter; controls heaviness of tails", "second shape parameter; controls skewness along with q", "third shape parameter; controls skewness along with p")
      for (i in 1:length(coefficients(lm_start))) {
        Parameters <- c(Parameters, names(list_start)[i+3])
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
