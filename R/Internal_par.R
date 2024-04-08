#R script containing internal functions essential for estimation of parameters
#------------------------------------------------------------------------------------------------------------
#' Internal Estimation Function for LVHDR Model
#'
#' This function estimates A, Beta, and Theta for a given value of K. It's designed
#' to be used internally within the package and supports parallel computing, handling
#' of time-dependent covariates, and convergence checking.
#'
#' @param Ymatlog Logical matrix derived from Y to reduce memory usage.
#' @param XZ Combined static and time-dependent covariates.
#' @param XZmat Matrix form of combined X and Z covariates for internal use.
#' @param indices_R Indices of observed entries in R.
#' @param lengths_R Number of observed time points for each subject.
#' @param Tp Number of time points.
#' @param K Number of factors to estimate.
#' @param cov_len Length of combined covariates.
#' @param px Number of static covariates.
#' @param N Number of subjects.
#' @param J Number of events.
#' @param Tinit Initial values for Theta.
#' @param Uinit Initial values for U (Gamma, A, Beta).
#'
#' @seealso \code{\link{lvhdr_est}} for the description of \code{Y}, \code{R}, \code{X}, \code{par}, \code{n.cores} and \code{Silent}.
#' @keywords internal
#'
#' @return A list containing the estimation results, including the success status, estimates of Theta and U, and the number of iterations.
lvhdr_est_fit_K <- function(Y, Ymatlog, R, XZ, X, XZmat, indices_R, lengths_R, Tp, K, cov_len, px, N, J, Tinit, Uinit, par, n.cores, Silent) {
  # Initialization
  success <- 1  # To store whether convergence is successful
  maxiter <- 300  # Maximum number of iterations
  Uit <- array(0, dim = c(J, Tp + cov_len + K, maxiter))  # Estimates of Gamma, A, and Beta
  Thetait <- array(0, dim = c(N, K, maxiter))  # Theta estimates
  niter <- 1  # Iteration counter
  tol <- 1e-6  # Tolerance for convergence
  reit <- 0  # Re-iteration counter

  # Initial estimates
  Thetait[, , 1] <- Tinit
  Uit[, , 1] <- Uinit

  # Objective function initialization
  obj.func <- function(tU, Theta) rcpp_objprep(Ymatlog, R, tU, Theta, XZmat, N, J, cov_len, K, Tp)
  Improvementsobj <- rep(0, maxiter - 1)
  objstore <- rep(0, maxiter)
  objstore[1] <- obj.func(t(Uit[, , 1]), as.matrix(Thetait[, , 1]))

  # Iteration loop
  while (niter == 1 || Improvementsobj[niter - 1] > tol) {
    if (!Silent) {
      print(niter)
    }
    niter <- niter + 1

    # Max iteration check
    if (niter > maxiter) {
      if (reit == 0) {
        reit <- 1
        print("Max iteration reached, retrying")
        niter <- 2  # Reset to continue the loop properly
        Thetait[, , 1] <- matrix(rtruncnorm(N * K, a = -1, b = 1, mean = 0, sd = 1), N, K)
        Uit[, , 1] <- matrix(rtruncnorm(J * (Tp + K + cov_len), a = -1, b = 1, mean = 0, sd = 1), J, Tp + K + cov_len)
        objstore[1] <- obj.func(t(Uit[, , 1]), as.matrix(Thetait[, , 1]))
      } else {
        print("Max iteration reached, does not converge")
        success <- 0
        break
      }
    }

    # Estimation update
    Uit[, , niter] <- t(Uest.func(Y, R, Uit[, , niter - 1], as.matrix(Thetait[, , niter - 1]), XZmat, indices_R, lengths_R, Tp, K, cov_len, N, J, par, n.cores))
    Thetait[, , niter] <- t(Thetaest.func(Y, R, as.matrix(Thetait[, , niter - 1]), Uit[, , niter], XZ, Tp, K, cov_len, N,J,  par, n.cores))

    # Finite check for Theta
    if (any(!is.finite(Thetait[, , niter]))) {
      print(sprintf("K = %d, Theta not finite, stopping", K))
      success <- 0
      break
    }

    # Objective value update
    objstore[niter] <- obj.func(t(Uit[, , niter]), as.matrix(Thetait[, , niter]))
    Improvementsobj[niter - 1] <- objstore[niter] - objstore[niter - 1]

    # Handle NA in objective improvement
    if (is.na(Improvementsobj[niter - 1])) {
      print(sprintf("K = %d, improvement NA, does not converge", K))
      success <- 0
      break
    }
  }

  # Finalize estimates
  final_iter <- ifelse(niter > maxiter, maxiter, niter)
  tBeta <- t(GetBetafromU.func(Uit[, , final_iter], Tp, cov_len))
  hats <- lvhdr_norm(J, N, cov_len, px, t(GetAfromU.func(Uit[, , final_iter], Tp, cov_len, K)), Thetait[, , final_iter], GetGammafromU.func(Uit[, , final_iter], Tp), tBeta, X)
  Thetahat <- hats$Thetastar
  Uhat <- if (cov_len > 0) {
    GetU.func(hats$Gammastar, t(hats$tBetastar), t(hats$tAstar))
  } else {
    GetU.func(hats$Gammastar, matrix(0, nrow = 0, ncol = 0), t(hats$tAstar))
  }

  # Return results
  return(list(
    success = success,
    Thetahat = Thetahat,
    Uhat = Uhat,
    Uit = Uit,
    Thetait = Thetait,
    Improvementsobj = Improvementsobj,
    niter = final_iter,
    K = K,
    objstore = objstore,
    lastobj = objstore[final_iter]
  ))
}
#----------------------------------------------------------------------------------------------------------------
#' Compute Initial Values Using SVD-based Algorithm
#'
#' This function calculates initial values for parameters using a Singular Value Decomposition (SVD)-based algorithm.
#' It's intended for internal use and supports parallel computing for efficient operation.
#'
#'@seealso \code{lvhdr_est_fit_K} for the description of parameters.
#' @keywords internal
#'
#' @return A list containing initial values for U (Gamma, A, Beta) and T (Theta).
SVDinit.func <- function(Y, R, XZmat, indices_R, lengths_R, Tp, K, cov_len, N, J, par, n.cores) {
  epsilon <- 0.01  # Constant used for bounding values

  #Initialization
  tildeM <- matrix(0, nrow = N, ncol = J)
  Gammainit <- matrix(0, nrow = J, ncol = Tp)

  # Compute proportion of observed responses across time periods
  phat <- colSums(R) / N

  # Transform Y into -1, 0, 1 based on value and missingness
  lijts <- 2 * Y - 1
  lijts[is.na(lijts)] <- 0

  # Apply SVD to transformed Y for each time point
  for (t in 1:Tp) {
    svdresult <- svd(lijts[, , t], nu = J, nv = J)
    tildeK <- max(K + 1, which(svdresult$d > 2 * sqrt(N * phat[t])))
    indices <- 1:tildeK
    tildelijts <- as.matrix(svdresult$u[, indices]) %*% diag(svdresult$d[indices], nrow = length(indices)) %*% t(svdresult$v[, indices])

    # Compute intermediate M matrix based on tildelijts
    M <- matrix(0, nrow = N, ncol = J)
    indices <- which((tildelijts >= -1 + epsilon) & (tildelijts <= 1 - epsilon))
    M[indices] <- invxi.func(0.5 * (tildelijts[indices] + 1))
    M[(tildelijts < -1 + epsilon)] <- invxi.func(epsilon)
    M[(tildelijts > 1 - epsilon)] <- invxi.func(1 - epsilon)

    Gammainit[, t] <- colSums(M) / N
    tildeM <- tildeM + M - matrix(Gammainit[, t], byrow = TRUE, nrow = N, ncol = J)
  }

  tildeM <- tildeM / Tp  # Normalize tildeM by the number of timepoints

  # SVD of tildeM to obtain initial Theta and A
  svdresult <- svd(tildeM, nu = J, nv = J)
  Tinit <- matrix(sqrt(N) * svdresult$u[, 1:K], nrow = N, ncol = K)

  Ainit <- matrix(0, nrow = J, ncol = K)
  for (k in 1:K) {
    Ainit[, k] <- svdresult$d[k] * svdresult$v[, k]
  }
  Ainit <- Ainit / sqrt(N)  # Normalize Ainit

  # Obtain initial Beta values through GLM, if applicable
  if (cov_len > 0) {
    Betainit <- t(Betaestglm.func(Y, R, Gammainit, Ainit, Tinit, XZmat, indices_R, lengths_R, Tp, K, cov_len, N, J, par, n.cores))
    # Transpose for the special case of cov_len = 1
    if (cov_len == 1) {
      Betainit <- t(Betainit)
    }
  } else {
    Betainit <- matrix(0, nrow = 0, ncol = 0)
  }

  Uinit <- GetU.func(Gammainit, Betainit, Ainit)

  return(list(Uinit = Uinit, Tinit = Tinit))
}
#----------------------------------------------------------------------------------------------------------------
#' Beta Estimation for Initial Values
#'
#' This function estimates the initial values for Beta parameters using Generalized Linear Models (GLMs).
#' It supports parallel processing and is designed for internal package use.
#'@seealso \code{lvhdr_est_fit_K} for the description of parameters.
#' @keywords internal
#'
#' @return A matrix containing the estimated Beta coefficients for each event.
Betaestglm.func <- function(Y, R, Gamma, A, Theta, XZmat, indices_R, lengths_R, Tp, K, cov_len, N, J, par, n.cores) {
  if (par) {
    out <- parallel::mcmapply(
      function(j) as.vector(Betaparglm.func(t(Y[, j, ]), indices_R, lengths_R, Gamma[j, ], A[j, ], Theta, XZmat, N, cov_len, K, Tp)),
      1:J, mc.cores = n.cores
    )
  } else {
    out <- sapply(
      1:J, function(j) as.vector(Betaparglm.func(t(Y[, j, ]), indices_R, lengths_R, Gamma[j, ], A[j, ], Theta, XZmat, N, cov_len, K, Tp))
    )
  }
  return(matrix(out, nrow = length(out) / J, ncol = J))
}

#' Beta Parameter Estimation with GLM
#'
#' This internal function performs the estimation of Beta parameters for a single event using GLM.
#' It's utilized by \code{Betaestglm.func} for either sequential or parallel processing.
#'
#' @param yris Transposed Y matrix for the jth event.
#' @param gammaj Gamma coefficients for the jth event.
#' @param aj A coefficients for the jth event.
#' @seealso \code{lvhdr_est_fit_K} for the description of other parameters.
#'
#' @keywords internal
#'
#' @return A vector of coefficients estimated by the GLM.
Betaparglm.func <- function(yris, indices_R, lengths_R, gammaj, aj, Theta, XZmat, N, cov_len, K, Tp) {
  yj <- yris[!is.na(yris)]

  Thetaaj <- Theta %*% aj
  id <- rep(1:N, times = lengths_R)

  Eforj <- rcpp_BetagetEforj(indices_R - 1, lengths_R, XZmat, N, Tp, cov_len)$Eforj
  glmfit <- glm.fit(Eforj, yj, family = binomial(), offset = Thetaaj[id] + gammaj[indices_R])

  coef <- glmfit$coefficients

  return(coef)
}

#--------------------------------------------------------------------------------------------------------------
#' Update the Values of U(Gamma, A, and Beta)
#'
#' This function updates the values of U based on the current estimates.
#' It supports parallel processing for efficiency.
#'
#' @param U Current estimates of Gamma, A, and Beta parameters.
#' @param Theta Current estimates Theta.
#' @seealso \code{lvhdr_est_fit_K} for the description of other parameters.
#'
#' @keywords internal
#'
#' @return A matrix containing the updated value of U.
Uest.func <- function(Y, R, U, Theta, XZmat, indices_R, lengths_R, Tp, K, cov_len, N, J, par, n.cores) {
  if (par) {
    out <- parallel::mcmapply(
      function(j) as.vector(Upar.func(t(Y[, j, ]), indices_R, lengths_R, U[j, ], Theta, XZmat, N, cov_len, K, Tp)),
      1:J, mc.cores = n.cores
    )
  } else {
    out <- sapply(
      1:J, function(j) as.vector(Upar.func(t(Y[, j, ]), indices_R, lengths_R, U[j, ], Theta, XZmat, N, cov_len, K, Tp))
    )
  }
  return(matrix(out, nrow = length(out) / J, ncol = J))
}

#' Parallel Function for U Parameter Update
#'
#' Performs the update of U parameters for a single event. It calculates gradients and
#' Hessians to find the optimal update direction and magnitude for U parameters.
#'
#' @param yris Transposed Y matrix for the jth event.
#' @param uj Current estimates of the jth row of U for the event.
#' @param Theta Current estimates Theta.
#'
#'@seealso \code{lvhdr_est_fit_K} for the description of other parameters.
#'
#' @keywords internal
#'
#' @return A vector of updated U parameters for the event.
Upar.func <- function(yris, indices_R, lengths_R, uj, Theta, XZmat, N, cov_len, K, Tp) {
  yj <- yris[!is.na(yris)]
  uj_len <- length(uj)

  # Adjust indices for Rcpp and calculate necessary matrices for optimization
  touse <- rcpp_getEforjEujXiM0(uj, indices_R - 1, lengths_R, Theta, XZmat, N, Tp, K, cov_len, uj_len)
  Eforj <- touse$Eforj
  Euj <- touse$Euj
  XiM0 <- touse$XiM0

  # Calculate gradient and Hessian
  grad <- rcpp_cal_grad(yj, XiM0, Eforj, N, uj_len)
  Hes <- rcpp_cal_Hes(XiM0, Eforj, uj, N, length(XiM0))

  # Determine descent direction
  drt <- get_drt.func(Hes, grad, uj_len)

  # Update jth row of U.
  out <- rcpp_getnewpar_func(XiM0, yj, N, grad, drt, uj, Eforj, rep(0, length(XiM0)))
  return(out)
}

#-----------------------------------------------------------------------------------------------------------------
#' Update Theta Parameters
#'
#' This function updates the Theta parameters based on the current estimates of U, Theta, and other model parameters.
#' It supports parallel processing for efficient computation.
#'
#' @param U Current estimates of Gamma, A, and Beta parameters.
#' @param Theta Current estimates Theta.
#'
#'@seealso \code{lvhdr_est_fit_K} for the description of other parameters.
#'
#' @keywords internal
#'
#' @return A matrix of updated Theta parameters for each subject.
Thetaest.func = function(Y,R,Theta,U,XZ,Tp,K,cov_len,N,J,par,n.cores){
  Gamma = GetGammafromU.func(U,Tp)
  A = as.matrix(GetAfromU.func(U,Tp,cov_len,K) )
  Beta = as.matrix(GetBetafromU.func(U,Tp,cov_len) )

  if(par){
    out =parallel::mcmapply(function(i) as.vector(Thetapar.func(t(Y[i,,]), R[i,],Theta[i,],A,Gamma,Beta,XZ,K,cov_len,i, J)) ,1:N,mc.cores = n.cores )
  }else{
    out= sapply(1:N, function(i) as.vector(Thetapar.func(t(Y[i,,]), R[i,],Theta[i,],A,Gamma,Beta,XZ,K,cov_len,i, J)) )
  }

  return(out)
}


#' Parallel Function for Theta Parameter Update
#'
#' Computes the update for Theta parameters for a single subject using current model estimates.
#'
#' @param yrjs Transposed Y matrix for the ith subject.
#' @param Ri Missing indicator vector for the ith subject.
#' @param thetai Current Theta estimate for the ith subject.
#' @param A Current estimate of A.
#' @param Gamma Current estimate of Gamma.
#' @param Beta Current estimate of Beta.
#' @param i Subject index.
#'
#'#'@seealso \code{lvhdr_est_fit_K} for the description of other parameters.
#'
#' @keywords internal
#'
#' @return A vector of updated Theta parameters for the subject.
Thetapar.func = function(yrjs, Ri,thetai,A,Gamma,Beta,XZ,K,cov_len,i, J){
  yi = yrjs[!is.na(yrjs)]

  id = rep(1:J, each=sum(Ri))
  # Prepare matrices and offsets for optimization
  tAfori= t(A[id,])
  ind = yrjs
  ind[!is.na(ind)]=1
  tGammaprep = t(Gamma)*ind
  if (cov_len>0){
    tBetaXZiprep = t(Beta%*%XZ[i,,])*ind
    tooff = tGammaprep[!is.na(tGammaprep)] +tBetaXZiprep[!is.na(tBetaXZiprep)]
  }else{
    tooff = tGammaprep[!is.na(tGammaprep)]
  }

  #Compute gradient and Hessian
  eUi  = as.vector(thetai%*%tAfori)+tooff
  XiM0=xi.func(eUi)
  grad<- rcpp_cal_grad(yi,XiM0, t(tAfori), J,K)
  Hes<-rcpp_cal_Hes(as.vector(XiM0), t(tAfori), thetai, J, length(XiM0))
  drt <- get_drt.func(Hes,grad,K)

  #Update parameter by line search
  out = rcpp_getnewpar_func(XiM0,yi,J,grad,drt, thetai,t(tAfori ),tooff)
  return(out)
}

#-----------------------------------------------------------------------------------------------------------------
#Simple supporting functions for parameters estimation
#projection function
#'@keywords internal
proj.func = function(out,K){
  constant=5*sqrt(K)
  value = sqrt(sum(out^2))
  if ( value >constant){
    out=out*constant/value
  }
  return(out)
}

#Function to get XZmat
#'@keywords internal
GetXZmat.func = function(XZ,cov_len,N,Tp){
  if(cov_len>0){
    # Initialize an empty matrix with dimensions p * (N*Tp), such that XZmat[, ((i-1)*Tp+1):(i*Tp)] refers to the p \times Tp matrix for the ith individual
    XZmat <- matrix(nrow = cov_len, ncol = N * Tp)
    # Fill the matrix using a loop over N
    for (i in 1:N){
      XZmat[, ((i-1)*Tp+1):(i*Tp)] <- XZ[i,,]
    }
  }else{
    XZmat = matrix(0, nrow=0, ncol=0)
  }
  return(XZmat)
}

#Function to create Ymat
#'@keywords internal
GetYmat.func = function(Y,N,J,Tp){
  Ymat <- matrix(nrow = N, ncol = J * Tp)
  # Fill the matrix using a loop over J
  for (j in 1:J){
    Ymat[, ((j-1)*Tp+1):(j*Tp)] <- Y[,j,]
  }
  return(Ymat)
}


#Function to get descend direction from Hessian matrix and gradient
#'@keywords internal
get_drt.func = function(Hes,grad,par_len){
  if (det(Hes)!=0){
    drt = as.vector(-solve(Hes)%*%grad)
  }else{
    drt = as.vector(-solve(Hes - 0.1*diag(par_len))%*%grad )
  }
  return(drt)
}

#inverse function of the logistic function
#'@keywords internal
invxi.func = function(p){
  out = log(p/(1-p))
  return(out)
}


#This script contains the function to generate data for simulation
#Function to combine Gamma,Beta, A to form U
#'@keywords internal
GetU.func = function(Gamma,Beta,A){
  if(nrow(Beta)>0){
    return(cbind(Gamma,Beta,A))
  }else{
    return(cbind(Gamma,A))
  }

}

#Function to extract Gamma given U
#'@keywords internal
GetGammafromU.func = function(U,Tp){
  return(U[,1:Tp])
}


#Function to extract Beta given U
#'@keywords internal
GetBetafromU.func = function(U,Tp,cov_len ){
  if(cov_len >0){
    return(U[,(Tp+1):(Tp+cov_len)])
  }else{
    return(matrix(0, nrow=0, ncol=0))
  }
}

#Function to extract A given U
#'@keywords internal
GetAfromU.func = function(U,Tp,p,Kstar){
  return(U[,(Tp+p+1):(Tp+p+Kstar)])
}

#Function to combine X and Z
#'@keywords internal
CombXZ.func = function(X,Z,px,pz,Tp,cov_len,N){
  if(cov_len>0){
    XZ = array(0,c(N,cov_len,Tp) )
    if(px>0){
      XZ[,1:px,] = X
    }
    if(pz>0){
      XZ[,(px+1):cov_len,] = Z
    }
  }else{
    XZ = array(0, dim = c(0,0,0))
  }

  return(XZ)
}

#logit function
#'@keywords internal
xi.func = function(x){
  expx = exp(x)
  out=expx/(1 + expx)
  #Make sure it does not explode
  out[x>700]=1
  out[x< -800]=0
  return(out)
}

