# LVHML
A Latent Variable Approach to Learning High-dimensional Multivariate longitudinal Data, as proposed by Lee at al. (2025+).

To install this package in R, run the following commands:  

```R
library(devtools) 
install_github("Arthurlee51/LVHML")
```

## Usage 
lvhml_est(Y, R, X ) 


## Brief description
The LVHML package implements the estimator developed by Lee et al. (2024+), which analyses High-dimensional multivariate longitudinal data via a latent variable framework. The main function lvhml_est() performs parameter estimation, latent dimension determination, and, optionally, estimation of the asymptotic variance and standard errors of covariates' regression coefficients in a binary response setting, as in simulation and data analysis. For detailed parameter descriptions, usage instructions, and examples, execute the following after installation:

```R
?lvhml_est
```

## Reference 
Sze Ming Lee, Yunxiao Chen, and Tony Sit. **Pairwise Comparisons without Stochastic Transitivity: Model, Theory and Applications**. 2025+. (manuscript)
