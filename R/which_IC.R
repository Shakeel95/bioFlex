#'General Purpose Information Criterion Function
#'
#'General purpose function for computing one of six preset information criteria. This function is called by flex_select , and is used to rank distributions in order of best overall fit. 
#'@param data An mandatory data frame or maxtrix containing the variables in the model.
#'@param model An object of class "mle2".
#'@param IC String variable indicating which information criteria should be outputted. Currently, the following are supported: “AIC”, “adjBIC”, “BIC”, “CAIC”, “AIC3”, and “AICc”.
#'@references Dziak, John et al. "Sensitivity And Specificity Of Information Criteria." Technical Report Series #12-119 (2017)
#'@export

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
