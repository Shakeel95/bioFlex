% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gamma_prob.R
\name{gamma_flexprob}
\alias{gamma_flexprob}
\title{Conditional probability for Gamma distribution}
\usage{
gamma_flexprob(K, model, features, visualise = TRUE, xlim, draws = 5)
}
\arguments{
\item{K}{Value for which P(Y<k | X) is computed.}

\item{model}{An object of class "mle2" produced using the function gamma_flexfit.}

\item{features}{A numeric vector specifying the value of covriates at which the conditional probability should be evaluated; the covariates in the vector should appear in the same order as they do in the model. Where a model does not depend on covariates the argument may be left blank.}

\item{visualise}{Logical. If TRUE (the default) the conditional distribution is plotted at P(Y<k | x) is shaded.}

\item{xlim}{Numeric vectors of length 2, giving the coordinate range of the dependent variable.}

\item{draws}{The number of random draws from multivariate random normal representing correlated parameters. If parameter correlation is not required draws should be set to zero.}
}
\description{
gamma_flexprob returns the conditional probability P(Y<k | X) of a model fitted via the function gamma_flexfit; where λ has been specified to be a function of covariates the required value should be specified using the ‘features’ parameter. The function includes a procedure for visualizing the conditional probability. gamma_flexprob also allows for the correlation of estimated parameters via the Cholesky decomposition of the variance-covariance matrix.
}
\details{
This function uses the the most common parametrization of the Gamma distribution.The probability probability density function is used is:

f(y) = λ^α/Γ(α)•y^α-1 exp(-λy)

The function returns:

P(Y<k | X) = (1/Γ(α))γ(α,kλ)

α may be a function of covariates; in which case, the cannonical log link function is used.
}
\references{
Kempthorne. "Parameter Estimation Fitting Probability Distributions Method Of Moments." MIT 18.443 (2015)

R. V. Hogg and A. T. Craig (1978) Introduction to Mathematical Statistics, 4th edition. New York: Macmillan. (See Section 3.3.)
}
