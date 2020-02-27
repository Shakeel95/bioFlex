# bioFlex

The package bioFlex implements the framework developed in [Jones et al. (2014)](https://onlinelibrary.wiley.com/doi/full/10.1002/hec.3178) for the estimation of a variable’s full conditional distribution. bioFlex was developed to estimate data generating processes in Health Economics but can be used with any strictly positive response variable; the package uses eleven distributions nested by the Generalized Beta distribution of the Second Kind to capture higher level features such as conditional skewness and conditional kurtosis.The package’s main use is the estimation of the conditional distribution, however once this has been obtained bioFlex can be used to calculate useful features of the distribution such as E(Y|X) and P(Y<k|X).
