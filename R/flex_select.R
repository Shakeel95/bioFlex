#' Select best fitting distribution
#'
#' flex_select is used to rank the overall distributional fit of a set of candidate distributions to a dataset. By default 11 distributions nested by  the Generalized Beta distribution of the Second Kind are considered. The ranking is based on an IC selected by the user. 
#'@param formula An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#'@param data An optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which flex_select is called.
#'@param IC String variable indicating which information criteria should be used to rank the proposed distributions. Currently, the following are supported: “AIC”, “adjBIC”, “BIC”, “CAIC”, “AIC3”, and “AICc”.
#'@param distribution A vector of string variables indicating which distributions should be ranked. Currently, the following are supported: "exp", "lnorm", "gamma", "weibull", "lomax", "fisk", "ggamma", "betapr", "sinmad", "dagum", and "gb2". By default, all distributions are ranked.
#'@details Since all distributions estimated are special cases of the Generalized Beta distribution of the Second Kind (GB2), comparing the information criteron of a distribution with a second distribution is equivalent to performing a likelihood ratio test as long as the second is a special case of the first. Non-nested distributions can also be compared in the same fashion by first comparing both to a distribution which nests them; usually the GB2.
#'@details It is important to understand that different ICs imply different α levels in the likelihood ratio test; different ICs therefore imply different preferences for over distributions with high and low numbers of parameters. 
#'@details The AIC criteration is good for choosing models with predictive power. However, the criterion is not consistent and is biased towards overparameterized models. The criteration takes the form: 
#'@details AIC = -2l(θ) + 2p
#'@details The adjusted BIC (adjBIC) criteration may penalize highly parametric models more than more than the BIC or less than the AIC, depending on the sample size; the likely model selection mistake therefore depends on the sample size. The criteration takes the form:
#'@details adjBIC = -2l(θ) +log([n+2]/24)p
#'@details The BIC criteration favors the selection of a parsimonious model. Although BIC is consistent, it often selects underrepresented models. The criteration takes the form:
#'@details BIC = -2l(θ) +log(n)p
#'@details The corrected AIC (CAIC) criteration likewise favors the selection of a parsimonious model and is most likely to select an underparameterized model. The CAIC is in fact more likely to select less parameterized models than the BIC. Dziak et al. (2017) note that additional parameter in the CAIC was selected somewhat arbitrarily, and the criterion therefore has no clear advantage over the BIC. The criteration takes the form:
#'@details CAIC = -2l(θ) + [log(n)+1]p
#'@details The AICc criteration was originally developed for rime series data, and applies a slightly heavier penalty to model complexity than the AIC. Model selection results will be similar to those under the AIC, as long as the sample size is not too large relative to the number of parameters. The criteration takes the form:
#'@details CAIC = AIC + (k+1)(k+2)/(n-k-2)
#'@details (k is the number of included regression coefficients, including an intercept)
#'@details The AIC3 criteration applies a heavier penalty to additional parameters than the AIC. Dziak et al. (2017) note that despite good performance in simulations, there is the AIC3 has little theoretical basis. The criteration takes the form:
#'@details AIC3 = -2l(θ) + 3p
#'@references Dziak, John et al. "Sensitivity And Specificity Of Information Criteria." Technical Report Series #12-119 (2017)
#'@export

flex_select <- function(formula, data = NULL, IC = "AIC", distributions = c("exp", "lnorm", "gamma", "weibull", "lomax", "fisk", "ggamma", "betapr", "sinmad", "dagum", "gb2")) {

  #===================#
  # first stage code  #
  #===================#

  # construct data frame equivalen to formula
  lm_test <- lm(formula = formula, data = data)
  data <- lm(formula = formula, data = data, method = "model.frame")
  # check if more than 30 covariates are used
  if(ncol(data) > 31)stop("At this time weibull_fit() accepts at most 30 covariates")

  CDF <- ecdf(data[,1])

  #================================================================#
  # linear predictor has intercept and is a function of covariates #
  #================================================================#

  if (isTRUE("(Intercept)" %in% names(coefficients(lm_test)))) {

    # exponential distribution
    if (isTRUE("exp" %in% distributions)) {
      exp <- tryCatch(exp_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (!suppressWarnings(is.na(exp))) {
        lambda <- 0
        for (i in 2:ncol(data)) {
          lambda <- lambda + as.numeric(tidy(exp)[i,2])*mean(data[,i])
        }
        lambda <- exp(lambda + as.numeric(tidy(exp)[1,2]))
        preview <- function(x) {
          dexp(x, rate = lambda)
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "exponential distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        exp <- which_IC(data = data, model = exp, IC = IC)
        print("the exponential distribution has been estimated")
      } else {
        print("the exponential distribution could not be estimated")
      }
    } else {
      exp <- NA
    }

    # log-normal distribution
    if (isTRUE("lnorm" %in% distributions)) {
      lnorm <- tryCatch(lnorm_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (!suppressWarnings(is.na(lnorm))) {
        params <- as.data.frame(tidy(lnorm))
        mu <- 0
        for (i in 2:ncol(data)) {
          mu <- mu + params[i+1,2]*mean(data[,i])
        }
        mu <- mu + params[2,2]
        preview1 <- function(x) {
          dlnorm(x, meanlog = mu, sdlog = params[1,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "log-normal distributional fit")
        par(new = TRUE)
        plot(preview1, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        lnorm <- which_IC(data = data, model = lnorm, IC = IC)
        print("the log-normal distribution has been estimated")
      } else {
        print("the log-normal distribution could not be estimated")
      }
    } else {
      lnorm <- NA
    }

    # gamma distribution
    if (isTRUE("gamma" %in% distributions)) {
      gamma <- tryCatch(gamma_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.na(gamma))) {
        params <- as.data.frame(tidy(gamma))
        alpha <- 0
        for (i in 2:ncol(data)) {
          alpha <- alpha + params[i+1,2]*mean(data[,i])
        }
        alpha <- exp(alpha + params[2,2])
        preview <- function(x) {
          dgamma(x, shape = alpha, rate = params[1,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "gamma distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        gamma <- which_IC(data = data, model = gamma, IC = IC)
        print("the gamma distribution has been estimated")
      } else {
        print("the gamma distribution could not be estimated")
      }
    } else {
      gamma <- NA
    }

    # weibull distribution
    if (isTRUE("weibull" %in% distributions)) {
      weibull <- tryCatch(weibull_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.na(weibull))) {
        params <- as.data.frame(tidy(weibull))
        lambda <- 0
        for (i in 2:ncol(data)) {
          lambda <- lambda + params[i+1,2]*mean(data[,i])
        }
        lambda <- exp(lambda + params[2,2])
        preview <- function(x) {
          dweibull(x, shape = params[1,2], scale = lambda)
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "weibull distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        weibull <- which_IC(data = data, model = weibull, IC = IC)
        print("the weibull distribution has been estimated")
      } else {
        print("the weibull distribution could not be estimated")
      }
    } else {
      weibull <- NA
    }

    # lomax distribution
    if (isTRUE("lomax" %in% distributions)) {
      lomax <- tryCatch(lomax_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.na(lomax))) {
        params <- as.data.frame(tidy(lomax))
        lambda <- 0
        for (i in 2:ncol(data)) {
          lambda <- lambda + params[i+1,2]*mean(data[,i])
        }
        lambda <- exp(lambda + params[2,2])
        preview <- function(x) {
          dLomax(x, shape = lambda, scale = params[1,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "lomax distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        lomax <- which_IC(data = data, model = lomax, IC = IC)
        print("the lomax distribution has been estimated")
      } else {
        print("the lomax distribution has not been estimated")
      }
    } else {
      lomax <- NA
    }

    # fisk distribution
    if (isTRUE("fisk" %in% distributions)) {
      fisk <- tryCatch(fisk_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.na(fisk))) {
        params <- as.data.frame(tidy(fisk))
        alpha <- 0
        for (i in 2:ncol(data)) {
          alpha <- alpha + params[i+1,2]*mean(data[,i])
        }
        alpha <- exp(alpha + params[2,2])
        preview <- function(x) {
          dllogis(x, shape = params[1,2], scale = alpha)
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "fisk distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        fisk <- which_IC(data = data, model = fisk, IC = IC)
        print("the fisk distribution has been estimated")
      } else {
        print("the fisk distribution has not been estimated")
      }
    } else {
      fisk <- NA
    }

    # generalised gamma distribution
    if (isTRUE("ggamma" %in% distributions)) {
      ggamma <- tryCatch(ggamma_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.na(ggamma))) {
        params <- as.data.frame(tidy(ggamma))
        mu <- 0
        for (i in 2:ncol(data)) {
          mu <- mu + params[i+2,2]*mean(data[,i])
        }
        mu <- exp(mu + params[3,2])^{-1}
        preview <- function(x) {
          dStacy_gamma(x, mu = mu, alpha = params[1,2], phi = params[2,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "generalised gamma distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        ggamma <- which_IC(data = data, model = ggamma, IC = IC)
        print("wowee, the generalised gamma distribution has been estimated")
      } else
        print("the generalised gamma distribution has not been estimated")
    } else {
      ggamma <- NA
    }

    # beta prime distribution
    if (isTRUE("betapr" %in% distributions)) {
      betapr <- tryCatch(betapr_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.na(betapr))) {
        params <- as.data.frame(tidy(betapr))
        mu <- 0
        for (i in 2:ncol(data)) {
          mu <- mu + params[i+1,2]*mean(data[,i])
        }
        mu <- exp(mu + params[2,2])
        preview <- function(x) {
          dbetapr(x, shape1 = mu, shape2 = params[1,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "beta prime distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        betapr <- which_IC(data = data, model = betapr, IC = IC)
        print("the beta prime distribution has been estimated")
      } else {
        print("the beta prime distribution could not be estimated")
      }
    } else {
      betapr <- NA
    }

    # singh-maddala distribution
    if (isTRUE("sinmad" %in% distributions)) {
      sinmad <- tryCatch(sinmad_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.na(sinmad))) {
        params <- as.data.frame(tidy(sinmad))
        mu <- 0
        for (i in 2:ncol(data)) {
          mu <- mu + params[i+2,2]*mean(data[,i])
        }
        mu <- exp(mu + params[3,2])
        preview <- function(x) {
          dsinmad(x, scale = mu, shape1.a = params[1,2], shape3.q = params[2,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "singh-maddala distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        sinmad <- which_IC(data = data, model = sinmad, IC = IC)
        print("the singh-maddala distribution has been estimated")
      } else {
        print("odd, the singh-maddala distribution could not be estimated")
      }
    } else {
      sinmad <- NA
    }

    # dagum distribution
    if (isTRUE("dagum" %in% distributions)) {
      dagum <- tryCatch(dagum_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.na(dagum))) {
        params <- as.data.frame(tidy(dagum))
        mu <- 0
        for (i in 2:ncol(data)) {
          mu <- params[i+2,2]*mean(data[,i])
        }
        mu <- exp(mu + params[3,2])
        preview <- function(x) {
          ddagum(x, scale = mu, shape1.a = params[1,2], shape2.p = params[2,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "dagum distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        dagum <- which_IC(data = data, model = dagum, IC = IC)
        print("the dagum distribution has been estimated")
      } else {
        print("the dagum distribution has not been estimated")
      }
    } else {
      dagum <- NA
    }

    # generalised beta of the second kind distribution
    if (isTRUE("gb2" %in% distributions)){
      gb2 <- tryCatch(gb2_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.null(gb2))) {
        params <- as.data.frame(tidy(gb2))
        mu <- 0
        for (i in 2:ncol(data)) {
          mu <- mu + params[i+3,2]*mean(data[,i])
        }
        mu <- exp(mu + params[4,2])
        preview <- function(x) {
          dgb2(x, shape1 = params[1,2], scale = mu, shape2 = params[2,2], shape3 = params[3,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "generalised beta of the second kind distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        gb2 <- which_IC(data = data, model = gb2, IC = IC)
        print("wowzers gb2 has been estimated, comparing information criteria is now equivalent to performing a likelihood ratio test!")
      } else {
        print("oh no, the generalised beta of the second kind could not be estimated!")
      }
    } else {
      gb2 <- NA
    }

    #=====================================================#
    # location parameter is not a function of covariates  #
    #=====================================================#

  } else if (is.null(names(coefficients(lm_test)))) {

    # exponential distribution
    if (isTRUE("exp" %in% distributions)) {
      exp <- tryCatch(exp_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (!suppressWarnings(is.na(exp))) {
        lambda <- exp(as.numeric(tidy(exp)[1,2]))
        preview <- function(x) {
          dexp(x, rate = lambda)
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "exponential distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        exp <- which_IC(data = data, model = exp, IC = IC)
        print("the exponential distribution has been estimated")
      } else {
        print("the exponential distribution could not be estimated")
      }
    } else {
      exp <- NA
    }

    # log-normal distribution
    if (isTRUE("lnorm" %in% distributions)) {
      lnorm <- tryCatch(lnorm_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (!suppressWarnings(is.na(lnorm))) {
        params <- as.data.frame(tidy(lnorm))
        preview <- function(x) {
          dlnorm(x, meanlog = params[2,2], sdlog = params[1,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "log-normal distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        lnorm <- which_IC(data = data, model = lnorm, IC = IC)
        print("the log-normal distribution has been estimated")
      } else {
        print("the log-normal distribution could not be estimated")
      }
    } else {
      lnorm <- NA
    }

    # gamma distribution
    if (isTRUE("gamma" %in% distributions)) {
      gamma <- tryCatch(gamma_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.na(gamma))) {
        params <- as.data.frame(tidy(gamma))
        preview <- function(x) {
          dgamma(x, shape = params[1,2], rate = params[2,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "gamma distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        gamma <- which_IC(data = data, model = gamma, IC = IC)
        print("the gamma distribution has been estimated")
      } else {
        print("the gamma distribution could not be estimated")
      }
    } else {
      gamma <- NA
    }

    # weibull distribution
    if (isTRUE("weibull" %in% distributions)) {
      weibull <- tryCatch(weibull_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.na(weibull))) {
        params <- as.data.frame(tidy(weibull))
        preview <- function(x) {
          dweibull(x, shape = params[2,2], scale = params[1,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "weibull distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        weibull <- which_IC(data = data, model = weibull, IC = IC)
        print("the weibull distribution has been estimated")
      } else {
        print("sorry, the weibull distribution could not be estimated")
      }
    } else {
      weibull <- NA
    }

    # lomax distribution
    if (isTRUE("lomax" %in% distributions)) {
      lomax <- tryCatch(lomax_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.na(lomax))) {
        params <- as.data.frame(tidy(lomax))
        preview <- function(x) {
          dLomax(x, shape = params[1,2], scale = params[2,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "lomax distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        lomax <- which_IC(data = data, model = lomax, IC = IC)
        print("the lomax distribution has been estimated")
      } else {
        print("the lomax distribution has not been estimated")
      }
    } else {
      lomax <- NA
    }

    # fisk distribution
    if (isTRUE("fisk" %in% distributions)) {
      fisk <- tryCatch(fisk_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.na(fisk))) {
        params <- as.data.frame(tidy(fisk))
        preview <- function(x) {
          dllogis(x, scale = params[1,2], shape = params[2,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "fisk distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        fisk <- which_IC(data = data, model = fisk, IC = IC)
        print("hooray, the fisk distribution has been estimated!")
      } else {
        print("oh no, the fisk distribution could not be estimated!")
      }
    } else {
      fisk <- NA
    }

    # generalised gamma
    if (isTRUE("ggamma" %in% distributions)) {
      ggamma <- tryCatch(ggamma_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.na(ggamma))) {
        params <- as.data.frame(tidy(ggamma))
        preview <- function(x) {
          dStacy_gamma(x, alpha = params[1,2], mu = params[2,2], phi = params[3,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "generalised gamma distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        ggamma <- which_IC(data = data, model = ggamma, IC = IC)
        print("wow, the generalised gamma distribution has been estimated")
      } else {
        print("the generalised gamma distribution has not been estimated")
      }
    } else {
      ggamma <- NA
    }

    # beta prime distribution
    if (isTRUE("betapr" %in% distributions)) {
      betapr <- tryCatch(betapr_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.na(betapr))) {
        params <- as.data.frame(tidy(betapr))
        preview <- function(x) {
          dbetapr(x, shape1 = params[1,2], shape2 = params[2,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "beta prime distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        betapr <- which_IC(data = data, model = betapr, IC = IC)
        print("hot dang, the beta prime distribution has been estimated")
      } else {
        print("the beta prime dietribution could not be estimated")
      }
    } else {
      betapr <- NA
    }

    # singh-maddala distribution
    if (isTRUE("sinmad" %in% distributions)) {
      sinmad <- tryCatch(sinmad_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.na(sinmad))) {
        params <- as.data.frame(tidy(sinmad))
        preview <- function(x) {
          dsinmad(x, shape1.a = params[1,2], scale = params[2,2], shape3.q = params[3,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "singh-maddala distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        sinmad <- which_IC(data = data, model = sinmad, IC = IC)
        print("the singh-maddala distribution has been estimated")
      } else {
        print("the singh-maddala distribution not has been estimated")
      }
    } else {
      sinmad <- NA
    }

    # dagum distribution
    if (isTRUE("dagum" %in% distributions)) {
      dagum <- tryCatch(dagum_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.na(dagum))) {
        params <- as.data.frame(tidy(dagum))
        preview <- function(x) {
          ddagum(x, scale = params[2,2], shape1.a = params[1,2], shape2.p = params[3,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "dagum distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        dagum <- which_IC(data = data, model = dagum, IC = IC)
        print("the dagum distribution has been estimated")
      } else {
        print("the dagum distribution could not be estimated")
      }
    } else {
      dagum <- NA
    }

    # generalised beta of the second kind
    if (isTRUE("gb2" %in% distributions)) {
      gb2 <- tryCatch(gb2_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.na(gb2))) {
        params <- as.data.frame(tidy(gb2))
        preview <- function(x) {
          dgb2(x, shape1 = params[2,2], scale = params[1,2], shape2 = params[3,2], shape3 = params[4,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "generalised beta of the second kind distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        gb2 <- which_IC(data = data, model = gb2, IC = IC)
        print("wowzers gb2 has been estimated, comparing information criteria is now equivalent to performing a likelihood ratio test!")
      } else {
        print("oh no, the generalised beta of the second kind could not be estimated!")
      }
    } else {
      gb2 <- NA
    }


    #====================================================================#
    # linear predictor has no intercept but is a function of covariates  #
    #====================================================================#

  } else {

    # exponential distribution
    if (isTRUE("exp" %in% distributions)) {
      exp <- tryCatch(exp_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (!suppressWarnings(is.na(exp))) {
        lambda <- 0
        for (i in 2:ncol(data)) {
          lambda <- lambda + as.numeric(tidy(exp)[i-1,2])*mean(data[,i])
        }
        lambda <- exp(lambda)
        preview <- function(x) {
          dexp(x, rate = lambda)
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "exponential distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        exp <- which_IC(data = data, model = exp, IC = IC)
        print("the exponential distribution has been estimated")
      } else {
        print("the exponential distribution could not be estimated")
      }
    } else {
      exp <- NA
    }

    # log-normal distribution
    if (isTRUE("lnorm" %in% distributions)) {
      lnorm <- tryCatch(lnorm_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (!suppressWarnings(is.na(lnorm))) {
        params <- as.data.frame(tidy(lnorm))
        mu <- 0
        for (i in 2:ncol(data)) {
          mu <- mu + mean(data[,i])*params[i,2]
        }
        preview <- function(x) {
          dlnorm(x, meanlog = mu, sdlog = params[1,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "log-normal distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        lnorm <- which_IC(data = data, model = lnorm, IC = IC)
        print("wow, the log-normal distribution has been estimated")
      } else {
        print("the log-normal distribution could not be estimated")
      }
    } else {
      lnorm <- NA
    }

    # gamma distribution
    if (isTRUE("gamma" %in% distributions)) {
      gamma <- tryCatch(gamma_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.na(gamma))) {
        params <- as.data.frame(tidy(gamma))
        alpha <- 0
        for (i in 2:ncol(data)){
          alpha <- alpha + params[i,2]*mean(data[,i])
        }
        alpha <- exp(alpha)
        preview <- function(x) {
          dgamma(x, shape = alpha, rate = params[1,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "gamma distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        gamma <- which_IC(data = data, model = gamma, IC = IC)
        print("the gamma distribution has been estimated")
      } else {
        print("the gamma distribution could not be estimated")
      }
    } else {
      gamma <- NA
    }

    # weibull distribution
    if (isTRUE("weibull" %in% distributions)) {
      weibull <- tryCatch(weibull_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.na(weibull))) {
        params <- as.data.frame(tidy(weibull))
        lambda <- 0
        for (i in 2:ncol(data)) {
          lambda <- lambda + params[i,2]*mean(data[,i])
        }
        lambda <- exp(lambda)
        preview <- function(x) {
          dweibull(x, shape = params[1,2], scale = lambda)
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "weibull distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        weibull <- which_IC(data = data, model = weibull, IC = IC)
        print("the weibull distribution has been estimated")
      } else {
        print("the weibull distribution has not been estimated, plenty more distributions...")
      }
    } else {
      weibull <- NA
    }

    # lomax distribution
    if (isTRUE("lomax" %in% distributions)) {
      lomax <- tryCatch(lomax_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.na(lomax))) {
        params <- as.data.frame(tidy(lomax))
        lambda <- 0
        for (i in 2:ncol(data)) {
          lambda <- lambda + params[i,2]*mean(data[,i])
        }
        lambda <- exp(lambda)
        preview <- function(x) {
          dLomax(x, scale = lambda, shape = params[1,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "lomax distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        lomax <- which_IC(data = data, model = lomax, IC = IC)
        print("the lomax distribution has been estimated")
      } else {
        print("the lomax distribution could not be estimated")
      }
    } else {
      lomax <- NA
    }

    # fisk distribution
    if (isTRUE("fisk" %in% distributions)) {
      fisk <- tryCatch(fisk_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.na(fisk))) {
        params <- as.data.frame(tidy(fisk))
        alpha <- 0
        for (i in 2:ncol(data)) {
          alpha <- alpha + params[i,2]*mean(data[,i])
        }
        alpha <- exp(alpha)
        preview <- function(x) {
          dllogis(x, scale = alpha, shape = params[1,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "fisk distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        fisk <- which_IC(data = data, model = fisk, IC = IC)
        print("hooray, the fisk distribution has been estimated!")
      } else {
        print("oh no, the fisk distribution has not been estimated!")
      }
    } else {
      firsk <- NA
    }

    # generalised gamma distribution
    if (isTRUE("ggamma" %in% distributions)) {
      ggamma <- tryCatch(ggamma_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.na(ggamma))) {
        params <- as.data.frame(tidy(ggamma))
        mu <- 0
        for (i in 2:ncol(data)){
          mu <- mu + params[i+1,2]*mean(data[,i])
        }
        mu <- exp(mu)
        preview <- function(x) {
          dStacy_gamma(x, mu = mu, alpha = params[1,2], phi = params[2,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "generalised gamma distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        ggamma <- which_IC(data = data, model = ggamma, IC = IC)
        print("the generalised gamma distribution has been estimated")
      } else {
        print("the generalised gamma distribution could not be estimated")
      }
    } else {
      ggamma <- NA
    }

    # beta prime distribution
    if (isTRUE("betapr" %in% distributions)) {
      betapr <- tryCatch(betapr_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.na(betapr))) {
        params <- as.data.frame(tidy(betapr))
        mu <- 0
        for (i in 2:ncol(data)) {
          mu <- mu + params[i,2]*mean(data[,2])
        }
        mu <- exp(mu)
        preview <- function(x) {
          dbetapr(x, shape1 = mu, shape2 = params[1,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "beta prime distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        betapr <- which_IC(data = data, model = betapr, IC = IC)
        print("not long left, the beta prime distribution has been estimated")
      } else {
        print("not long left, however the beta prime distribution could not be estimated")
      }
    } else {
      betapr <- NA
    }

    # singh-maddala distribution
    if (isTRUE("sinmad" %in% distributions)) {
      sinmad <- tryCatch(sinmad_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.na(sinmad))) {
        params <- as.data.frame(tidy(sinmad))
        mu <- 0
        for (i in 2:ncol(data)) {
          mu <- mu + params[i+1,2]*mean(data[,i])
        }
        mu <- exp(mu)
        preview <- function(x) {
          dsinmad(x, scale = mu, shape1.a = params[1,2], shape3.q = params[2,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "singh-maddala distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        sinmad <- which_IC(data = data, model = sinmad, IC = IC)
        print("good job, the singh-maddala distribution has been estimated")
      } else {
        print("odd, the singh-maddala distribution could not be estimated")
      }
    } else {
      sinmad <- NA
    }

    # dagum distribution
    if (isTRUE("dagum" %in% distributions)) {
      dagum <- tryCatch(dagum_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.na(dagum))) {
        params <- as.data.frame(tidy(dagum))
        mu <- 0
        for (i in 2:ncol(data)) {
          mu <- mu + params[i+1,2]*mean(data[,i])
        }
        preview <- function(x) {
          ddagum(x, scale = mu, shape1.a = params[1,2], shape2.p = params[2,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "dagum distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        dagum <- which_IC(data = data, model = dagum, IC = IC)
        print("top model choice, the singh-maddala distribution has been estimated, this usually doesn't happen")
      } else {
        print("the dagum distribution has not been estimated, consider including an intercept in the linear predictor")
      }
    } else {
      dagum <- NA
    }

    # generalised beta of the second kind
    if (isTRUE("gb2" %in% distributions)) {
      gb2 <- tryCatch(gb2_flexfit(formula = formula, data = data, key = FALSE), error = function(e) {return(NA)})
      if (suppressWarnings(!is.na(gb2))) {
        params <- as.data.frame(tidy(gb2))
        mu <- 0
        for (i in 2:ncol(data)){
          mu <- mu + params[i+2,2]*mean(data[,2])
        }
        mu <- exp(mu)
        preview <- function(x) {
          dgb2(x, shape1 = params[1,2], scale = mu, shape2 = params[2,2], shape3 = params[3,2])
        }
        hist(data[,1], probability = TRUE, col = "grey", xlab = names(data)[1], main = "generalised beta of the second kind distributional fit")
        par(new = TRUE)
        plot(preview, xlim = c(min(data[,1]), max(data[,1])), lwd = 3, col = "red", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        gb2 <- which_IC(data = data, model = gb2, IC = IC)
        print("wowzers gb2 has been estimated, comparing information criteria is now equivalent to performing a likelihood ratio test!")
      } else {
        print("oh no, the generalised beta of the second kind could not be estimated!")
      }
    } else {
      gb2 <- NA
    }

  }

  ranking <- t(sort(data.frame(exp, lnorm, gamma, weibull, lomax, fisk, ggamma, betapr, sinmad, dagum,gb2)))
  colnames(ranking) <- IC
  return(ranking)
}
