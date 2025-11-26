# LVHML
A Latent Variable Approach to Learning High-dimensional Multivariate longitudinal Data, as proposed by Lee at al. (2025).

To install this package in R, run the following commands:

```R
library(devtools)
install_github("Arthurlee51/LVHML")
```

## Usage 
lvhml_est(Y, R, X ) 


## Brief description
The LVHML package implements the estimator developed by Lee et al. (2025+), which analyses High-dimensional multivariate longitudinal data via a latent variable framework. The main function lvhml_est() performs parameter estimation, latent dimension determination, and, optionally, estimation of the asymptotic variance and standard errors of covariates' regression coefficients in a binary response setting, as in simulation and data analysis. For detailed parameter descriptions, usage instructions, and examples, execute the following after installation:

```R
?lvhml_est
```
## Workflow
The Workflow folder contains the codes as well as the description of the workflow to obtain the numerical results in the article. Please refer to the README file in the folder for details.

## Reference 
Sze Ming Lee, Yunxiao Chen, and Tony Sit. **A Latent Variable Approach to Learning High-dimensional Multivariate Longitudinal Data**. *arXiv preprint* arXiv:2405.15053v2, 2025.
