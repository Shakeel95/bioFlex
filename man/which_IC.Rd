% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/which_IC.R
\name{which_IC}
\alias{which_IC}
\title{General Purpose Information Criterion Function}
\usage{
which_IC(data, model, IC)
}
\arguments{
\item{data}{An mandatory data frame or maxtrix containing the variables in the model.}

\item{model}{An object of class "mle2".}

\item{IC}{String variable indicating which information criteria should be outputted. Currently, the following are supported: “AIC”, “adjBIC”, “BIC”, “CAIC”, “AIC3”, and “AICc”.}
}
\description{
General purpose function for computing one of six preset information criteria. This function is called by flex_select , and is used to rank distributions in order of best overall fit.
}
\references{
Dziak, John et al. "Sensitivity And Specificity Of Information Criteria." Technical Report Series #12-119 (2017)
}
