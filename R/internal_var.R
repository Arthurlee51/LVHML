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
calasympvar.func <- function(X, Z, Thetahat, Uhat, R, J, N, Tp, ext,gamma_fix) {
  ujlen <- ncol(Uhat)
  #Ehat<-array(0,dim=c(N,ujlen,Tp))
  px <- ncol(X)
  pz <- dim(Z)[2]
  K <- ncol(Thetahat)

  UERhat <- UER.func(Uhat, Thetahat, X,Z, Tp, px,pz, J, N,ext,gamma_fix)

  Ehat <- compute_Ehat( Thetahat, X, Z,  ujlen,K,Tp, px,pz, N,ext, gamma_fix)
  
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
#Return a array that gives e_{it} u_{j}^{top} for every i, j and t.
#'@keywords internal
UER.func = function(U, Theta, X, Z,Tp,px,pz,J,N,ext,gamma_fix){
  out <- array(0, dim = c(N,J,Tp))
  IN = rep(1, N)
  K <- ncol(Theta)
  g_len <- if (gamma_fix) 1 else Tp
  b_len <- if (px > 0) if (ext) px * Tp else px else 0
  
  for ( t in 1:Tp){
    if(!gamma_fix){
      out[,,t]<- IN%*%t(U[,t])
    }else{
      out[,,t]<- t*IN%*%t(U[,1])
    }
    
    if(px>0){
      if(!ext){
        out[,,t] <- out[,,t] + X%*%t(U[,(g_len+1):(g_len+px)])
      }else{
        out[,,t] <- out[,,t] + X%*%t(U[,(g_len+(t-1)*px+1 ):(g_len+t*px)])
      }
    }
    
    if(pz>0){
      out[,,t] <- out[,,t] + Z[,,t]%*%t(U[,(g_len+b_len+1):(g_len+b_len+pz)])
    }
    
    if(!ext){
      out[,,t] <- out[,,t] + Theta%*%t(U[,(g_len+b_len+pz+1):ncol(U)])
    }else{
      out[,,t] <- out[,,t] + Theta%*%t(U[,(g_len+b_len+pz+(t-1)*K+1):(g_len+b_len+pz+t*K)])
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


#Function to compute Ehat
compute_Ehat <- function( Thetahat, X, Z,  ujlen,K,Tp, px,pz, N,ext, gamma_fix) {
  Ehat<-array(0,dim=c(N,ujlen,Tp))
  
  for ( t in 1:Tp){
    #Compute Ehat based on t.i.e. e_{it}
    if(!gamma_fix){
      Dprep = matrix(0, nrow=N,ncol=Tp)
      Dprep[,t]=1
    }else{
      Dprep <- matrix(t, nrow=N,ncol=1)
    }
    
    tobind <- Dprep
    if(px>0){
      if(!ext){
        tobind <- cbind(tobind ,X)
      }else{
        DXprep = matrix(0, nrow = N, ncol = Tp*px)
        DXprep[,((t-1)*px+1):(t*px)] <- X
        tobind <- cbind(tobind , DXprep)
      }
    }
    
    if(pz>0){
      tobind <- cbind(tobind ,Z[,,t])
    }
    
    if(!ext){
      Ehat[,,t]<- cbind(tobind, Thetahat)
    }else{
      DTheta <- matrix(0, nrow = N, ncol = Tp*K)
      DTheta[,((t-1)*K+1):(t*K)] <- Thetahat
      Ehat[,,t]<- cbind(tobind, DTheta )
    }
  }
  return(Ehat)
}
