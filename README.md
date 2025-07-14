# LVHML
A Latent Variable Approach to Learning High-dimensional Multivariate longitudinal Data, as proposed by Lee at al. (2025+).

To install this package in R, run the following commands after downloading LVHML_0.1.0.tar.gz from the repository:

```R
library(devtools) 
install.packages("LVHML_0.1.0.tar", repos = NULL, type = "source")
```

## Usage 
lvhml_est(Y, R, X ) 


## Brief description
The LVHML package implements the estimator developed by Lee et al. (2025+), which analyses High-dimensional multivariate longitudinal data via a latent variable framework. The main function lvhml_est() performs parameter estimation, latent dimension determination, and, optionally, estimation of the asymptotic variance and standard errors of covariates' regression coefficients in a binary response setting, as in simulation and data analysis. For detailed parameter descriptions, usage instructions, and examples, execute the following after installation:

```R
?lvhml_est
```

## Reference 
Sze Ming Lee, Yunxiao Chen, and Tony Sit. **A Latent Variable Approach to Learning
High-dimensional Multivariate longitudinal Data**. 2025+. (manuscript)
