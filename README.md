# bioFlex

# 1. Fitting models to real data with 'bioFlex' 

In Dr John Schorling’s Diabetes dataset, glycosolated haemoglobin (`glyhb`) above 7.0 is taken as a positive diagnosis of diabetes. A researcher may be interested in estimating a statistical model for `glyhb`, given a vector of covariates, in order to identify patients at risk of developing diabetes. 

In this example, the covariates considered are:  a dummy variable for the subject’s location     `Buckingham `, a dummy variable for the subject’s gender ` male `, height in inches ` height `,  weight in pounds `weight`, age `age`, waist size in inches `waist`, hip size in inches `hip`, and a series of dummy variables for the subjects build `small` and `large`.

### Data 

The `Diabetes` dataset consist of 17 variables on 373 subjects from 1046 subjects who were interviewed in a study to understand the prevalence of obesity, diabetes, and other cardiovascular risk factors in central Virginia for African Americans. 

```{r, echo=FALSE, results='asis'}
knitr::kable(head(Diabetes))
```

The response variable, `glyhb`,  is strictly positive with a heavy right tail, so modelling the data with a flexible parametric distribution is appropriate. It is good practice to take note of the skewness and the kurtosis, as these are good benchmarks for the estimated model. 

```{r}
# Load dataset 
attach(Diabetes)
# Calculate the skewness 
Skew(glyhb)
# Calculate the kurtosis 
Kurt(glyhb)
```

```{r, echo = FALSE}
plot(density(glyhb), lwd = 3, main = "Density Plot")
hist(glyhb,  main = "Histogram", xlab = "Glycosolated Hemoglobin")
```
