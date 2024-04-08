# LVHDR
A Latent Variable Approach to Learning Discrete-time High-dimensional Recurrent Event Data

To install this package in R, run the following commands:  

```R
library(devtools) 
install_github("Arthurlee51/LVHDR")
```

## Usage 
lvhdr_est(Y, R, X , Kset = 1:10, par = FALSE, n.cores = 1, Asymp = FALSE, Silent = FALSE, Z = array(0, dim = c(0, 0, 0))) 


## Brief description
This package implements the estimator proposed in Lee et al(2024+) for Learning Discrete-time High-dimensional Recurrent Event Data using a latent variable approach. See 
lvhdr_est.Rd in the man folder for description of paramters and example code. 

## Reference 
Sze Ming Lee, Yunxiao Chen, and Tony Sit. **A Latent Variable Approach to Learning Discrete-time High-dimensional Recurrent Event Data**. 2024+. (manuscript)
