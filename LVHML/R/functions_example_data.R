#' Function to generate example data.
#' @param J Number of event types.
#' @param nset sets the size of N relative to J, "s": N = 5J, "L":N = 10J
#' @param Kstar True latent dimension
#' @param ext Logical. Indicates whether the extension (section 2.4) allowing time-dependent loadings and coefficients on X is used. Default is FALSE.
#' @param gamma_fix Logical. Indicates whether constraint (section 2.5)  gamma_{jt} = t gamma_j  for j = 1,...,J is imposed. Default is FALSE.
#' @param gap Logical. Indicates whether the the setting (Section E.1) with sufficient eigengap for identification of latent factors is applied. Default is FALSE. 
#' @return A list that can be used as input for `lvhml_est` and true values of parameters for comparison.
#'         \item{Ustar}{A matrix storing the true values of U.}
#'         \item{Thetastar}{A matrix storing the true values of Theta.}
#'         \item{N,J,Tp}{Return values of N, J and Tp for referece.}
#' @seealso \code{\link{lvhml_est}} for the description of the observed data Y, R and X.
#'@export
lvhml_ex_data <- function(J=100,nset="s", Kstar =3, ext = FALSE, gamma_fix = FALSE, gap = FALSE){
  # Determine sample size
  N <- switch(nset,s = 5 * J,L = 10 * J)

  #Tp is the number of time-points, equivalent to T in the manucsript.
  Tp = 4
  #Generate parameters according to setting
  parameters <- genpar.func(Kstar, N, J, Tp, ext = ext, gamma_fix = gamma_fix, gap = gap)
  
  #Generate data based on paramters
  dataset = gendata.func(parameters$Lambda ,Kstar, N, J, Tp)

  return(list(Y = dataset$Y, R = dataset$R, X = parameters$X, Ustar = parameters$Ustar, Thetastar = parameters$Thetastar, Kstar = Kstar, N=N, J=J, Tp = Tp  ))
}


#---------------------------------------------------------------------------------------------------------------
#Supporting internal functions:
#' @keywords internal
#Function to generate parameters. 
#withZ: Logical. Indicates whether the time-dependent covariate Z_t is included. Default is FALSE.
genpar.func <- function(Kstar,N,J,Tp, withZ = FALSE, ext = FALSE, gamma_fix = FALSE,  gap = FALSE){
  #Generate Gamma, tranpose of A(tA),Theta, X and Beta. 
  g_len <- if(gamma_fix) 1 else Tp
  if(!gamma_fix){
    Gamma = matrix(runif(J * g_len, -1, 1), nrow = J, ncol = g_len)
  }else{
    Gamma = matrix(runif(J * g_len, -0.25, 0.25), nrow = J, ncol = g_len)
  }
  A_len <- if(ext) Kstar*Tp else Kstar
  tA = matrix(truncnorm::rtruncnorm(J * A_len, a = -3, b = 3), nrow = A_len, ncol = J)
  Theta = matrix(truncnorm::rtruncnorm(N * Kstar, a = -3, b = 3), N, Kstar)
  if(gap==TRUE){
    Theta = matrix(truncnorm::rtruncnorm(N * Kstar, a = -1, b = 1), N, Kstar)#
    Theta<- Theta*matrix( (1:Kstar)/2, nrow = N, ncol = Kstar, byrow = TRUE )
  }
  
    px = 5
    Xprep1 <- rbinom(N, 2, 0.5) + 1  #First dummy variable 
    Xprep2 <- rbinom(N, 2, 0.5) + 1  #Second dummy variable
    X = matrix(0, nrow = N, ncol = px)
    for (l in 1:2) {
      X[, l] <- (Xprep1 == l)
      X[,l+2] <- (Xprep2 == l)
    }
    X[, 5] = runif(N, -1, 1)

  b_len <- if(ext) px*Tp else px
  tBeta = matrix(runif(b_len * J, 0.5, 1), b_len, J)
  
  #Initialise Z and finalised transpose of V
  if(!withZ){
    pz = 0
    Z = array(0, c(N, pz, Tp))
    tVstar = matrix(0, pz, J)
  }else{
    pz <- 2
    Z = array(runif(N*pz*Tp, -1, 1),c(N, pz, Tp) )
    tVstar = matrix(runif(pz * J, -1, 1), pz, J)
  }
 
 
  # Normalization function
  Stars = lvhml_norm(J, N,  px, tA, Theta, Gamma, tBeta,X, Tp,ext, gamma_fix)
  tAstar = Stars$tAstar
  Thetastar = Stars$Thetastar
  tBetastar = Stars$tBetastar
  #Update tBetastar for easier comparison in FDR in later stage. Set half of coefficients to zero.
  #b_prep_len: Number of covariates(not dummy variables).
    b_prep_len <- if(ext) (px-2)*Tp else px-2
    matrix_prep <- matrix(sample(c(0,1),b_prep_len  * J , replace = TRUE) , b_prep_len, J)
    matrix_ind <- matrix(sample(c(0,1),b_len  * J , replace = TRUE) , b_len, J)
    if(!ext){
      matrix_ind[1,] <-  matrix_prep[1,]
      matrix_ind[2,] <- matrix_prep[1,]
      matrix_ind[3,]<- matrix_prep[2,]
      matrix_ind[4,]<- matrix_prep[2,]
      matrix_ind[5,]<- matrix_prep[3,]
    }else{
      for( t in 1:Tp){
        matrix_ind[(t-1)*px+1,] <-  matrix_prep[(t-1)*(px-2) + 1,]
        matrix_ind[(t-1)*px+2,] <- matrix_prep[(t-1)*(px-2) + 1,]
        matrix_ind[(t-1)*px+3,] <-  matrix_prep[(t-1)*(px-2) + 2,]
        matrix_ind[(t-1)*px+4,] <- matrix_prep[(t-1)*(px-2) + 2,]
        matrix_ind[(t-1)*px+5,] <- matrix_prep[(t-1)*(px-2) + 3,]
      }
    }
   tBetastar <-tBetastar*matrix_ind
   Gammastar = Stars$Gammastar

  # Initialization of Lambda and other components
  IN = matrix(1, nrow = N, ncol = 1)
  Lambda <- array(0, dim = c(N, J, Tp))

  #Generate lambba_ijt
  #Generate Gamma component based on the status of gamma_fix.
  for (t in 1:Tp) {
    if(!gamma_fix){
      gamma_part <- IN %*% Gammastar[, t]
    }else{
      gamma_part <- t*IN %*%t(Gammastar)
    }
  if(!ext){
      Lambda[, , t] = xi.func(Thetastar %*% tAstar + X %*% tBetastar + Z[,,t]%*%tVstar +   gamma_part)
  }else{
      Lambda[, , t] = xi.func(Thetastar %*% tAstar[((t-1)*Kstar+1):((t-1)*Kstar+Kstar),] + X %*% tBetastar[((t-1)*px+1):((t-1)*px+px),] + Z[,,t]%*%tVstar +   gamma_part )
  }
  }

  Ustar <- GetU.func(Gammastar, t(tBetastar), t(tAstar), t(tVstar))

  # Returning the output list
  out <- list(
    "Ustar" = Ustar,
    "Thetastar" = Thetastar,
    "Z" = Z,
    "X" = X,
    "Lambda" = Lambda,
    "px" = px,
    "pz" = pz,
    "Tp" = Tp,
    "N" = N
  )
  return(out)
}


#Function to generate data according to the parameters.
#' @keywords internal
# Function to generate data
gendata.func <- function(Lambda ,Kstar,N, J, Tp) {
  # Initialize an array to store y_{rij} values
  Y = array(NA, dim = c(N, J, Tp))

  # R is the matrix for GM indicator, ensuring at least one timepoint is observed
  # Overall GM rate is around 30 to 40%
  # Generate all possible combinations of GM indicators
  all_combinations <- expand.grid(rep(list(0:1), Tp))
  all_combinations <- as.matrix(all_combinations[-1, ])
  colnames(all_combinations) <- NULL
  rownames(all_combinations) <- NULL

  # Sample indices and create R matrix
  sampled_indices <- sample(nrow(all_combinations), N, replace = TRUE)
  R <- all_combinations[sampled_indices, ]

  # Generate Failure time data
  for (t in 1:Tp) {
    # FT: Failure at time t, accounting for the influence of C
    FT = matrix(0, nrow = N, ncol = J)
    Failed = (matrix(runif(N * J, 0, 1), nrow = N, ncol = J) <= Lambda[, , t])
    FT[Failed] = 1
    FT[R[, t] == 0, ] = NA  # Setting NA for unobserved timepoints
    Y[, , t] = FT
  }

  # Returning the output list with generated data
  out = list("Y" = Y, "R" = R)
  return(out)
}




