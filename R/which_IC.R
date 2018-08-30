#'General Purpose Information Criterion Function
#'
#'General purpose function for computing one of six preset information criteria.
#'@param data
#'@param model A fitted model object for which there exists a logLik method to extract the corresponding log-likelihood, or an object inheriting from class logLik.
#'@param IC ...

which_IC <- function(data, model, IC) {

  L <- as.numeric(logLik(model))
  n <- nrow(data)
  k <- as.numeric(nrow(tidy(model)))

  if (isTRUE("AIC" %in% IC)) {
    return(-2*L + 2*k)
  } else if (isTRUE("adjBIC" %in% IC)) {
    return(-2*L + log((n+2)/24)*k)
  } else if (isTRUE("BIC" %in% IC)) {
    return(-2*L + log(n)*k)
  } else if (isTRUE("CAIC" %in% IC)) {
    return(-2*L + (log(n) + 1)*k)
  } else if (isTRUE("AIC3" %in% IC)) {
    return(-2*L + 3*k)
  } else if (isTRUE("AICc" %in% IC)) {
    return(-2*L + 2*k + 2*((k+1)*(k+2))/(n-k-2))
  }

}
