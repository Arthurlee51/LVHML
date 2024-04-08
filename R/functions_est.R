#' @useDynLib LVHDR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' @title Estimation Function for LVHDR Model
#' @description Main function to implement the proposed estimator. This includes choosing the number of factors (K) and calculating asymptotic variances if required. The function is capable of handling time-dependent covariates.
#'
#' @param Y A N x J x Tp array storing occurrences of events (1 for occurred, 0 for no occurrence). Entries for unobserved data should be marked as NA.
#' @param R A N x Tp matrix representing the missing indicator (1 for observed, 0 for missing).
#' @param X A N x px matrix storing static covariates.
#' @param Kset A vector containing the candidate numbers of factors. The default value is 1:10.
#' @param par Logical. Indicates whether to use parallel computing. Default is FALSE.
#' @param n.cores The number of cores to use for parallel computing. Used in conjunction with par = TRUE. Default is 1.
#' @param Asymp Logical. Indicates whether to estimate the asymptotic variance. Default is FALSE.
#' @param Silent Logical. Option to suppress output showing progress. Default is FALSE.
#' @param Z A N x J x Tp array for time-dependent covariates. Defaults as empty array if not provided.
#'
#' @return A list containing the results of the estimation, including:
#'         \item{Thetahat}{Estimated Theta}
#'         \item{Gammahat}{Estimated Gamma}
#'         \item{Betahat}{Estimated Beta}
#'         \item{Ahat}{Estimated A}
#'         \item{Khat}{Estimated K}
#'         \item{ICs}{Information criteria for values of K in Kset}
#'         \item{Sigma_E}{Estimated asymptotic variance matrices for sqrt(N)*betahat_j (if Asymp = TRUE)}
#'         \item{SE}{The standard errors of Betahat (if Asymp = TRUE)}
#'
#' @examples
#' # Generate example dataset
#' set.seed(123)
#' dataset <- lvhdr_ex_data(J = 100, nset = "s", Kstar = 3)
#' attach(dataset)
#' out <- lvhdr_est(Y, R, X, Kset = 1:10, par = FALSE, n.cores = 1, Asymp = FALSE, Silent = FALSE, Z = array(0, dim = c(0, 0, 0)))
#'
#' @export
lvhdr_est <- function(Y, R, X , Kset = 1:10, par = FALSE, n.cores = 1, Asymp = FALSE, Silent = FALSE, Z = array(0, dim = c(0, 0, 0))) {
  # Initialize variables and lists for the estimation process
  all_output <- list()  # List storing output across all values of K
  lk <- length(Kset)  # Number of K values to be tested
  ICs <- rep(0, lk)
  lastobjs <- rep(0, lk)

  # Get dimensions of input matrices
  px <- ncol(X)
  pz <- dim(Z)[2]
  cov_len <- px + pz
  N <- dim(Y)[1]
  J <- dim(Y)[2]
  Tp <- dim(Y)[3]

  # Create XZ, combining static and time-dependent covariates
  XZ <- CombXZ.func(X, Z, px, pz, Tp, cov_len, N)
  XZmat <- GetXZmat.func(XZ, cov_len, N, Tp)

  # Convert Y to logical matrix to reduce memory usage
  Ymatlog <- GetYmat.func(Y, N, J, Tp) == 1

  # Indices of time points where r_{it} = 1
  indices_R <- unlist(apply(R, 1, function(row) which(row == 1)))
  lengths_R <- rowSums(R)  # Lengths of the time points

  # Loop through each K value in Kset for estimation
  for (i in seq_along(Kset)) {
    K <- Kset[i]
    if (!Silent) {
      print(sprintf("K = %d", K))
    }

    # Compute initial values using SVD-based approach
    init_list <- SVDinit.func(Y, R, XZmat, indices_R, lengths_R, Tp, K, cov_len, N, J, par, n.cores)

    # Estimation of parameters
    out <- lvhdr_est_fit_K(Y, Ymatlog, R, XZ, X, XZmat, indices_R, lengths_R, Tp, K, cov_len, px, N, J, init_list$Tinit, init_list$Uinit, par, n.cores, Silent)
    all_output <- append(all_output, list(out))  # Store results for each K

    # Calculate and store IC
    ICs[i] <- -2 * out$lastobj * N * J + K * max(N, J) * log(J*sum(R) / max(N, J))
    lastobjs[i] <- out$lastobj
  }

  # Select Khat that minimises IC
  Khat <- if (lk == 1) Kset[1] else Kset[which.min(ICs)]

  # Get the data for the chosen Khat and return relevant parameters
  chosen <- all_output[[which.min(ICs)]]
  Thetahat <- chosen$Thetahat
  Uhat <- chosen$Uhat
  Gammahat <- GetGammafromU.func(Uhat, Tp)
  Betahat <- GetBetafromU.func(Uhat, Tp, cov_len)
  Ahat <- GetAfromU.func(Uhat, Tp, cov_len, Khat)

  if (!Asymp) {
    return(list(Thetahat = Thetahat, Gammahat = Gammahat, Betahat = Betahat, Ahat = Ahat, Khat = Khat, ICs = ICs))
  } else {
    AsymVarhatstore <- calasympvar.func(X, Z, Thetahat, Uhat, R, J, N, Tp)
    Sigma_E <- AsymVarhatstore[(Tp + 1):(Tp + cov_len), (Tp + 1):(Tp + cov_len), ]
    SE <- t(sapply(1:J, function(j) sqrt(diag(Sigma_E[,,j])))) / sqrt(N)
    return(list(Thetahat = Thetahat, Gammahat = Gammahat, Betahat = Betahat, Ahat = Ahat, Khat = Khat, ICs = ICs, Sigma_E = Sigma_E, SE = SE))
}
}

#-----------------------------------------------------------------------------------------------------------
#' @title Normalize Model Parameters
#' @description  This function normalizes the model parameters (A, Theta, Gamma, and Beta) to satisfy identifiability conditions.
#'
#' @param J Number of events or outcome variables.
#' @param N Number of subjects.
#' @param cov_len Total number of covariates (static and time-dependent).
#' @param px Number of static covariates.
#' @param tA Transposed A matrix.
#' @param Theta Theta matrix.
#' @param Gamma Gamma matrix.
#' @param tBeta Transposed Beta matrix representing covariate(static and time-dependent) effects.
#' @param X Matrix of static covariates.
#'
#' @return A list containing normalized versions of tA (tAstar), Theta (Thetastar), Gamma (Gammastar), and tBeta (tBetastar).
#'
#'@export
lvhdr_norm = function(J,N,cov_len,px,tA,Theta,Gamma,tBeta,X){
  # Initial checks for NA values in critical parameters
  if (any(is.na(tA)) || any(is.na(Theta))) {
    message("NA detected in tA or Theta, normalization cannot proceed.")
    return(list(tAstar = NaN, Thetastar = NaN, Gammastar = NaN, tBetastar = NaN))
  }

    IN = matrix(1, nrow =N, ncol =1)
    INZi <- if (cov_len > 0) cbind(IN, X) else IN

    # Calculate Tp from Gamma and adjust Theta
    Tp =ncol(Gamma)
    TILTH = Theta - INZi %*% solve(t(INZi)%*%INZi) %*% t(INZi) %*% Theta

    # Projection of Theta and calculation of adjusted Gamma and Beta
    Q = as.matrix(solve(t(INZi)%*%INZi)%*%t(INZi)%*%Theta%*%tA)
    Gammastar = Gamma + matrix(Q[1,], nrow=J, ncol=Tp)
    if(cov_len>0){
      tBetastar = tBeta
      if(px>0){
        tBetastar[1:px,] = tBetastar[1:px,] +Q[2:(px+1),]
      }
    }else{
      tBetastar=matrix(0, nrow=0, ncol=0)
    }

    # Normalizing A and Theta
    SigmaJA = tA%*%t(tA)/J
    SigmaNT = t(TILTH)%*%TILTH/N
    object = eigen(expm::sqrtm(SigmaJA)%*%SigmaNT%*%expm::sqrtm(SigmaJA))
    HNJ = t(solve(expm::sqrtm(SigmaJA))%*%object$vectors)
    tAstar = HNJ%*%tA #tAstar%*%t(tAstar)/J Gives Identity matrix
    Thetastar= TILTH%*%solve(HNJ)#t(Thetastar)%*%Thetastar/N gives diagonal matrix
  return(list( "tAstar" = tAstar, "Thetastar"= Thetastar, "Gammastar"=Gammastar, "tBetastar"= tBetastar))
}

#-------------------------------------------------------------------------------------------------------------



