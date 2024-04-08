#' Function to generate example data.
#To be updated. Add parameters and output descriptions as well as making a better than that connects with the package name
#' @param J Number of event types.
#' @param nset sets the size of N relative to J, "s": N = 5J, "L":N = 10J
#' @param Kstar True latent dimension
#' @return A list structure that can be used as input for `lvhdr_est` and true values of parameters for comparison.
#'         \item{Ustar}{A matrix storing the true values of U.}
#'         \item{Thetastar}{A matrix storing the true values of Theta.}
#'         \item{N,J,Tp}{Return values of N, J and Tp for referece.}
#' @seealso \code{\link{lvhdr_est}} for the description of the observed data Y, R and X.
#'@export
lvhdr_ex_data <- function(J=100,nset="s", Kstar =3){
  if(nset =="s"){
    N = 5*J
  }
  if(nset =="L"){
    N = 10*J
  }

  #Tp is the number of time-points, equivalent to T in the manucsript.
  Tp = 4
  #Generate parameters according to setting
  parameters <- genpar.func(Kstar, N, J, Tp)

  #Generate data based on paramters
  dataset = gendata.func(parameters$Lambda ,Kstar, N, J, Tp)

  return(list(Y = dataset$Y, R = dataset$R, X = parameters$X, Ustar = parameters$Ustar, Thetastar = parameters$Thetastar, Kstar = Kstar, N=N, J=J, Tp = Tp  ))
}


#---------------------------------------------------------------------------------------------------------------
#Supporting internal functions:
#' @keywords internal
#Function to generate parameters
genpar.func <- function(Kstar,N,J,Tp){
  Gamma = matrix(runif(J * Tp, -1, 1), nrow = J, ncol = Tp)
  tA = matrix(truncnorm::rtruncnorm(J * Kstar, a = -3, b = 3), nrow = Kstar, ncol = J)
  Theta = matrix(truncnorm::rtruncnorm(N * Kstar, a = -3, b = 3), N, Kstar)
  px = 4
  pz = 0
  cov_len = px + pz
  Xprep = rbinom(N, 3, 0.5) + 1
  X = matrix(0, nrow = N, ncol = px)
  for (l in 1:3) {
    X[, l] <- (Xprep == l)
  }
  X[, 4] = runif(N, -1, 1)
  Z = array(0, c(N, pz, Tp))
  tBeta = matrix(runif(cov_len * J, -1, 1), cov_len, J)

  # Normalization function
  Stars = lvhdr_norm(J, N, cov_len, px, tA, Theta, Gamma, tBeta,X)
  tAstar = Stars$tAstar
  Thetastar = Stars$Thetastar
  tBetastar = Stars$tBetastar
  Gammastar = Stars$Gammastar

  # Initialization of Lambda and other components
  IN = matrix(1, nrow = N, ncol = 1)
  Lambda <- array(0, dim = c(N, J, Tp))
  XZ = CombXZ.func(X, Z, px, pz, Tp, cov_len, N)

  #Generate lambba_ijt
  for (t in 1:Tp) {
      Lambda[, , t] = xi.func(Thetastar %*% tAstar + XZ[, , t] %*% tBetastar + IN %*% Gammastar[, t])
  }

  Ustar <- GetU.func(Gammastar, t(tBetastar), t(tAstar))

  # Returning the output list
  out <- list(
    "Ustar" = Ustar,
    "Thetastar" = Thetastar,
    "Z" = Z,
    "X" = X,
    "Lambda" = Lambda,
    "cov_len" = cov_len,
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




