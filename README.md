# LVHDR
A Latent Variable Approach to Learning Discrete-time High-dimensional Recurrent Event Data, as proposed by Lee at al. (2024+).

To install this package in R, run the following commands:  

```R
library(devtools) 
install_github("Arthurlee51/LVHDR")
```

## Usage 
lvhdr_est(Y, R, X ) 


## Brief description
The LVHDR package implements the estimator developed by Lee et al. (2024+), which analyses Discrete-time High-dimensional Recurrent Event Data via a latent variable framework. The main function lvhdr_est() performs parameter estimation, latent dimension determination, and, optionally, estimation of the asymptotic variance and standard errors of covariates' regression coefficients. For detailed parameter descriptions, usage instructions, and examples, execute the following after installation:

```R
?lvhdr_est
```

## Reference 
Sze Ming Lee, Yunxiao Chen, and Tony Sit. **A Latent Variable Approach to Learning Discrete-time High-dimensional Recurrent Event Data**. 2024+. (manuscript)
