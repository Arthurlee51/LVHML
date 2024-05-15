#R script containing internal functions essential for estimation of asymptotic variance
#----------------------------------------------------------------------------------------------------------------
#' Compute Phi_j
#'
#' Computes the Phi_j matrix for the calculation of asymptotic variance, based on EEstarstore and ujE values.
#'
#' @param EEstarstore A 4 dimensional array storing pre-computed values of e_{ir}e_{ir}^{t} for all i and r.
#' @param ujE uj times e_{ir} values for a specific j.
#' @param ujlen Length of the uj.
#'
#'@seealso \code{lvhml_est_fit_K} for the description of other parameters.
#'
#'#' @keywords internal
#'
#' @return The computed Phi_j matrix.
calPhi_j.func <- function(EEstarstore, ujE, R, N, J, Tp, ujlen) {
  rhoipp <- -exp(ujE) / (1 + exp(ujE))^2
  Evarrhoipp <- matrix(0, nrow = N, ncol = Tp)

  for (t in 1:Tp) {
    Evarrhoipp[, t] <- R[, t] * rhoipp[, t]
  }

  Phi_j <- matrix(0, nrow = ujlen, ncol = ujlen)
  for (t in 1:Tp) {
    for (i in 1:N) {
      Phi_j <- Phi_j + Evarrhoipp[i, t] * EEstarstore[,, i, t]
    }
  }

  Phi_j <- Phi_j / N
  return(Phi_j)
}


#' Compute Asymptotic Variance of uj for all j.
#'
#' Calculates the asymptotic variance of the U parameter for each event (j).
#'
#'
#' @keywords internal
#'
#'#'@seealso \code{lvhml_est_fit_K} for the description of parameters.
#'
#' @return An array containing the estimated asymptotic variance for each event.
calasympvar.func <- function(X, Z, Thetahat, Uhat, R, J, N, Tp) {
  ujlen <- ncol(Uhat)
  Ehat<-array(0,dim=c(N,ujlen,Tp))
  px <- ncol(X)
  pz <- dim(Z)[2]
  cov_len <- px + pz
  XZ <- CombXZ.func(X, Z, px, pz, Tp, cov_len, N)

  UERhat <- UER.func(Uhat, Thetahat, XZ, Tp, cov_len, J, N)

  for ( t in 1:Tp){
    Dprep = matrix(0, nrow=N,ncol=Tp)
    Dprep[,t]=1
    if(cov_len>0){
      Ehat[,,t]<- cbind(Dprep,XZ[,,t],Thetahat )
    }else{
      Ehat[,,t]<- cbind(Dprep, Thetahat)
    }
  }
  EEhatstore <- calEEstore.func(Ehat, ujlen, N, Tp)

  AsymVarhatstore <- array(0, dim = c(ujlen, ujlen, J))

  for (j in 1:J) {
    ujEhat <- UERhat[, j, ]
    Phi_j_hat <- calPhi_j.func(EEhatstore, ujEhat, R, N, J, Tp, ujlen)
    AsymVarhat <- solve(-Phi_j_hat)
    AsymVarhatstore[,, j] <- AsymVarhat
  }

  return(AsymVarhatstore)
}
#-------------------------------------------------------------------------------------------------------------
#Simple supporting functions related to the computation of asymptotic variance
#Return a array that gives uj^t e_{ir} for every i and j
#'@keywords internal
UER.func = function(U, Theta, XZ,Tp,cov_len,J,N){
  out <- array(0, dim = c(N,J,Tp))
  IN = rep(1, N)
  for ( t in 1:Tp){
    if(cov_len>0){
      out[,,t]<-   cbind(IN,XZ[,,t],Theta)%*%t(U[,c(t,(Tp+1):ncol(U))] )
    }else{
      out[,,t]<-   cbind(IN,Theta)%*%t(U[,c(t,(Tp+1):ncol(U))] )
    }
  }
  return(out)
}

#Function to calculate e_{ir}e_{ir}^{\top} for all combinations of i and r.
calEEstore.func = function(E, ujlen, N,Tp){
  EEstore = array(0,dim = c(ujlen, ujlen,N,Tp))
  for ( i in 1:N){
    for ( r in 1:Tp){
      EEstore[,,i,r]<- E[i,,r]%*%t(E[i,,r])
    }
  }
  return(EEstore)
}
