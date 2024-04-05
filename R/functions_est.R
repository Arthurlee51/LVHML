#est.func: Main function to implement the proposed estimator, including choosing number of factors K and calculating asymptotic variances.
# Function to estimate A, Theta, and K
#Z is for the extension handling time-dependent covariates.
lvhdr_est <- function(Y, R,  X = matrix(0, nrow = 0, ncol = 0), Tp, Kset=1:10, par = FALSE, n.cores = 1, Asymp = FALSE, Z = array(0, dim = c(0, 0, 0))) {
  # Initialize variables and lists for the estimation process
  all_output = list()#list storing output across all values of K.
  lk = length(Kset)  # Number of K values to be tested
  ICs = rep(0, nrow = lk)
  lastobjs = rep(0, lk)

  # Get dimensions of input matrices
  px = ncol(X)
  pz = dim(Z)[2]
  cov_len = px + pz
  N = dim(Y)[1]
  J = dim(Y)[2]

  # Create XZ, combining static and time-dependent covariates
  XZ = CombXZ.func(X, Z, px, pz, Tp, cov_len, N)
  XZmat <- GetXZmat.func(XZ, cov_len, N, Tp)

  # Convert Y to logical matrix to reduce memory usage
  Ymatlog <- GetYmat.func(Y, N, J, Tp) == 1

  # Count of non-NA Y elements
  sumY <- length(Y[!is.na(Y)])

  #
  indices_R <- unlist(apply(R, 1, function(row) which(row == 1)))  # Indices of timepoints where r_{it} = 1
  lengths_R <- rowSums(R)  # Lengths of the timepoints

  # Loop through each K value in Kset for estimation
  for (i in 1:length(Kset)) {
    K = Kset[i]
    print(sprintf("K = %d", K))
    #Get initial values
    # #Tinit: Initial values for Theta,
    # #Unint: Initial values for U.
    val = 0.5

    #Compute initial values using SVD-based approach
    init_list <- SVDinit.func(Y, R, XZmat,indices_R,lengths_R, Tp, K,cov_len, N, J, par, n.cores )

    #Estimation of parameters.
    out = lvhdr_est_fit_K(Y, Ymatlog, R, XZ, X, XZmat,indices_R, lengths_R, Tp, K, cov_len, px,N, J, init_list$Tinit, init_list$Uinit ,par, n.cores)
    all_output <- append(all_output, list(out))  # Store results for each K

    # Calculate IC for each K
    ICs[i] <- -2 * out$lastobj * N * J + K * max(N, J) * log(sumY / max(N, J))
    lastobjs[i] <- out$lastobj
  }

  # Select winning K and return output
  if (lk == 1) {
    Khat = Kset[1]
  } else {
    Khat <- Kset[which.min(ICs)]
  }


  #Get the data for the chosen Khat and return relevant parameters.
  chosen <- all_output[[which.min(ICs)]]
  Thetahat <- chosen$Thetahat
  Uhat <- chosen$Uhat
  Gammahat <- GetGammafromU.func(Uhat,Tp)
  Betahat <- GetBetafromU.func(Uhat,Tp,cov_len)
  Ahat <- GetAfromU.func(Uhat,Tp,cov_len,Khat)
  #Return output
  if(!Asymp){
    return(list(Thetahat = Thetahat, Gammahat = Gammahat, Betahat = Betahat, Ahat = Ahat, Khat = Khat, all_output = all_output,ICs = ICs))
  }else{
    AsymVarhatstore <- calasympvar.func(X,Z, Thetahat , Uhat,R ,J,N,Tp)
    #Get Sigma_E, the asymptotic sqrt(N)*betaj, j = 1,...,J
    Sigma_E <- AsymVarhatstore[(Tp+1):(Tp+cov_len),(Tp+1):(Tp+cov_len), ]
    #Compute standard error for individual betas
    #SE: standard error
    SE <- t(sapply(1:J, function(j) sqrt(diag(Sigma_E[,,j]))))/sqrt(N)
    return(list(Thetahat = Thetahat, Gammahat = Gammahat, Betahat = Betahat, Ahat = Ahat, Khat = Khat, all_output = all_output,ICs = ICs,Sigma_E=Sigma_E, SE = SE  ))
  }
}


# estfit_K.func: Function to estimate A, Beta, and Theta for a given K
lvhdr_est_fit_K <- function(Y, Ymatlog, R, XZ, X, XZmat, indices_R, lengths_R, Tp, K, cov_len, px, N, J,Tinit, Uinit , par, n.cores) {
#  est.func <- function(Y, Ymatlog, R, XZ, X, XZmat, Tp, K, N, J, par, n.cores) {
  # Initialization of variables
  success = 1  # To store whether convergence is successful
  maxiter = 300  # Maximum number of iterations
  Uit = array(0, dim = c(J, Tp + cov_len + K, maxiter))  # Store estimates of Gamma, tA and tBeta
  Thetait = array(0, dim = c(N, K, maxiter))
  niter = 1  # Iteration counter

  # Apply initial values
  Thetait[, , 1] = Tinit
  Uit[, , 1] = Uinit

  # Set tolerance for convergence
  tol = 10^(-6)
  reit = 0  # Re-iteration counter

  # Define the objective function
  obj.func <- function(tU, Theta) rcpp_objprep(Ymatlog, R, tU, Theta, XZmat, N, J, cov_len, K, Tp)

  # Store the improvement and objective values
  Improvementsobj = rep(0, maxiter - 1)
  objstore = rep(0, maxiter)
  objstore[1] = obj.func(t(Uit[, , 1]), as.matrix(Thetait[, , 1]))

  # Iteration part
  while (niter == 1 || Improvementsobj[niter - 1] > tol) {
    print(niter)
    niter = niter + 1

    # Check for max iteration and reiterate if necessary
    if (niter == (maxiter + 1)) {
      if (reit == 0) {
        reit = reit + 1
        print("Max iteration reached, retrying")
        niter = 1
        Thetait[, , 1] = matrix(rtruncnorm(N * K, a = -1, b = 1, N, K))
        Uit[, , 1] = matrix(rtruncnorm(J * (Tp + K + cov_len), a = -1, b = 1), J, Tp + K + cov_len)
        objstore[1] = obj.func(t(Uit[, , 1]), as.matrix(Thetait[, , 1]))
        next
      } else {
        print("Max iteration reached")
        niter = niter - 1
        success = 0
        break
      }
    }

    # Estimate U and Theta in each iteration
    Uit[, , niter] = t(Uest.func(Y, R, Uit[, , niter - 1], as.matrix(Thetait[, , niter - 1]), XZmat, indices_R, lengths_R, Tp, K, cov_len, N, J, par, n.cores))
    Thetait[, , niter] = t(Thetaest.func(Y, R, as.matrix(Thetait[, , niter - 1]), Uit[, , niter], XZ, Tp, K, cov_len, N, par, n.cores))

    # Break the loop if Theta is not finite
    if (any(!is.finite(Thetait[, , niter]))) {
      print(sprintf("K=%d, Theta not finite", K))
      break
    }

    # Calculate and store the objective value
    objstore[niter] = obj.func(t(Uit[, , niter]), as.matrix(Thetait[, , niter]))
    Improvementsobj[niter - 1] = objstore[niter] - objstore[niter - 1]

    # Handle NA in improvement of objective value
    if (is.na(Improvementsobj[niter - 1])) {
      if (reit == 0) {
        reit = reit + 1
        print(sprintf("K=%d, improvement NA, retrying", K))
        niter = 1
        Thetait[, , 1] = matrix(rtruncnorm(N * K, a = -1, b = 1), N, K)
        Uit[, , 1] = matrix(rtruncnorm(J * (Tp + K + cov_len), a = -1, b = 1), J, Tp + K + cov_len)
        objstore[1] = obj.func(t(Uit[, , 1]), Thetait[, , 1])
        next
      } else {
        print(sprintf("K=%d, improvement NA, does not converge", K))
        success = 0
        break
      }
    }

  }

  # Normalize and retrieve final estimates of A, Beta, Gamma, and Theta
  tBeta = t(GetBetafromU.func(Uit[, , niter], Tp, cov_len))
  hats = norm.func(J, N, cov_len, px, t(GetAfromU.func(Uit[, , niter], Tp, cov_len, K)), Thetait[, , niter], GetGammafromU.func(Uit[, , niter], Tp), tBeta, X)
  Thetahat = hats$Thetastar
  Uhat = if (cov_len > 0) {
    GetU.func(hats$Gammastar, t(hats$tBetastar), t(hats$tAstar))
  } else {
    GetU.func(hats$Gammastar, matrix(0, nrow = 0, ncol = 0), t(hats$tAstar))
  }

  # Return the list of results including success status, Theta estimates, and objective values
  return(list(
    success = success,
    Thetahat = Thetahat,
    Uhat = Uhat,
    Uit = Uit,
    Thetait = Thetait,
    Improvementsobj = Improvementsobj,
    niter = niter,
    K = K,
    objstore = objstore,
    lastobj = objstore[niter]
  ))
}

#-----------------------------------------------------------------------------------------------------------
#Functions for computing initial values using the SVD-based algorithm
#Method 3:SVD-based algorithm + glm
 SVDinit.func = function(Y, R, XZmat,indices_R,lengths_R, Tp, K,cov_len, N, J, par, n.cores ){

   #Epsilon is a constant.
   epsilon <- 0.01

   #Compute phat: proportion of observed responses across time periods
   phat = colSums(R)/N

   #Compute lijts:Transform y_ijt into -1,0,1 based on value and missingness
   lijts <- 2*Y - 1
   lijts[is.na(lijts)] <- 0

   #initiate tildeM, Gammainit
   tildeM = matrix(0, nrow=N, ncol=J)
   Gammainit = matrix(0, nrow=J, ncol = Tp)

   #Apply SVD to lijts
   for ( t in 1:Tp){
     svdresult <- svd(lijts[,,t], nu = J, nv = J)
     #Compute tilde lijts for each t
     tildeK = max(K+1, which(svdresult$d > 2*sqrt(N*phat[t]) ))
     indices <- 1:tildeK
     tildelijts = as.matrix(svdresult$u[,indices])%*%diag(svdresult$d[indices], nrow=length(indices))%*%t(svdresult$v[,indices])


     #Now compute indices for M
     M = matrix(0, nrow=N, ncol=J)
     indices <- which( (tildelijts >= -1+epsilon) & (tildelijts <= 1 - epsilon) )
     M[indices]= invxi.func( 0.5*(tildelijts[indices]+1))
     M[(tildelijts < -1 + epsilon)] = invxi.func( epsilon )
     M[(tildelijts > 1 - epsilon) ] = invxi.func(1 - epsilon )

     Gammainit[,t]<- colSums(M)/N

     tildeM = tildeM  + M - matrix(Gammainit[,t], byrow = TRUE, nrow=N, ncol=J)
   }

   #Divide by number of timepoints
   tildeM = tildeM/Tp

   #Now get Theteinit and Ainit through svd
   svdresult <- svd(tildeM,  nu = J, nv = J)
   Tinit <- matrix(sqrt(N)*svdresult$u[,1:K], nrow=N,ncol=K)

   Ainit = matrix(0, nrow=J, ncol=K)
   for( k in 1:K){
     Ainit[,k]<- svdresult$d[k]*svdresult$v[,k]
   }
   Ainit =Ainit/sqrt(N)

   if(cov_len>0){
     #Now get initial values for beta(and v) through glm calculations
     Betainit <- t(Betaestglm.func(Y,R,Gammainit,Ainit,Tinit,XZmat,indices_R,lengths_R,Tp,K,cov_len,N,J,par,n.cores))
     #Tranpose for the special case cov_len=1
     if(cov_len==1){
       Betainit = t(Betainit)
     }
   }else{
     Betainit <- matrix(0,nrow=0, ncol=0)
   }

   Uinit <- GetU.func(Gammainit,Betainit,Ainit)
   return(list(Uinit = Uinit, Tinit = Tinit))
 }

#Beta est function(For generating initial values only)
Betaestglm.func = function(Y,R,Gamma,A,Theta,XZmat,indices_R,lengths_R,Tp,K,cov_len,N,J,par,n.cores){
  if(par){
    out =parallel::mcmapply(function(j) as.vector(Betaparglm.func(t(Y[,j,]),indices_R,lengths_R,Gamma[j,],A[j,],Theta,XZmat,N,cov_len,K,Tp)),1:J,mc.cores = n.cores )
  }else{
    out= sapply(1:J, function(j) as.vector(Betaparglm.func(t(Y[,j,]),indices_R,lengths_R,Gamma[j,],A[j,],Theta,XZmat,N,cov_len,K,Tp)) )
  }
  return(out)
}

#yris = t(Y[,j,]), gammaj = Gamma[j,], aj = A[j,]
Betaparglm.func = function(yris,indices_R,lengths_R,gammaj,aj,Theta,XZmat,N,cov_len,K,Tp){
  yj = yris[!is.na(yris)]

  Thetaaj<- Theta%*%aj

  id = rep(1:N, lengths_R)

  Eforj =rcpp_BetagetEforj( indices_R-1,lengths_R ,  XZmat, N, Tp,  cov_len)$Eforj
  glmfit <- glm.fit(Eforj,yj, family = binomial(), offset =  Thetaaj[id] + gammaj[indices_R] )

  coef = glmfit$coefficients

  return(coef)
}


#-------------------------------------------------------------------------------------------------------
#Other supporting functions for parameters estimation
#Update the values of Gamma, A and Beta(U)
Uest.func = function(Y,R,U,Theta,XZmat,indices_R,lengths_R,Tp,K,cov_len,N,J,par,n.cores){
    if(par){
      out =parallel::mcmapply(function(j) as.vector(Upar.func(t(Y[,j,]),indices_R,lengths_R,U[j,],Theta,XZmat,N,cov_len,K,Tp)),1:J,mc.cores = n.cores )
    }else{
      out= sapply(1:J, function(j) as.vector(Upar.func(t(Y[,j,]),indices_R,lengths_R,U[j,],Theta,XZmat,N,cov_len,K,Tp)) )
    }
  return(out)
}

#yris = t(Y[,j,]), uj = U[j,]
Upar.func = function(yris,indices_R,lengths_R,uj,Theta,XZmat,N,cov_len,K,Tp){
  yj = yris[!is.na(yris)]
  uj_len = length(uj)

  #indices_R-1 used to make sure that indices match Rcpp setting
  touse <-rcpp_getEforjEujXiM0( uj, indices_R-1,lengths_R, Theta, XZmat, N, Tp,K,cov_len,uj_len)
  Eforj <- touse$Eforj
  Euj <- touse$Euj
  XiM0 <- touse$XiM0

  #Calculate gradient
  grad<- rcpp_cal_grad(yj,XiM0,Eforj,N,uj_len )

  Hes <- rcpp_cal_Hes(XiM0, Eforj, uj, N, length(XiM0))

  #Get descend direction
  drt <- get_drt.func(Hes,grad,uj_len )
  #--------------------------------------------------------------------------------------------
  out = rcpp_getnewpar_func(XiM0,yj,N,grad,drt,uj,Eforj,rep(0,length(XiM0) ) )
  return(out)
}


Thetaest.func = function(Y,R,Theta,U,XZ,Tp,K,cov_len,N,par,n.cores){
  out = matrix(0, nrow= N, ncol =K)
  Gamma = GetGammafromU.func(U,Tp)
  A = as.matrix(GetAfromU.func(U,Tp,cov_len,K) )
  Beta = as.matrix(GetBetafromU.func(U,Tp,cov_len) )
    if(par){
      out =parallel::mcmapply(function(i) as.vector(Thetapar.func(t(Y[i,,]), R[i,],Theta[i,],A,Gamma,Beta,XZ,K,cov_len,i)) ,1:N,mc.cores = n.cores )
    }else{
      out= sapply(1:N, function(i) as.vector(Thetapar.func(t(Y[i,,]), R[i,],Theta[i,],A,Gamma,Beta,XZ,K,cov_len,i)) )
    }
  return(out)
}



#yrjs = t(Y[i,,]),Ri = R[i,]
Thetapar.func = function(yrjs, Ri,thetai,A,Gamma,Beta,XZ,K,cov_len,i){
  #Create dataframe
  yi = yrjs[!is.na(yrjs)]
  #id, that are to repeat in matrix
  id = rep(1:J, each=sum(Ri))
  #Matrix for the repeated covariates
  #Aj after taking repetitions into account.
  tAfori= t(A[id,])
  #Offset to be considered.
  #ind: Indicator for both Gamma and Beta XZ.
  ind = yrjs
  ind[!is.na(ind)]=1
  tGammaprep = t(Gamma)*ind
  #Get Beta%*%XZ for ALL r first. Then try to generalise it.
  if (cov_len>0){
    tBetaXZiprep = t(Beta%*%XZ[i,,])*ind
    tooff = tGammaprep[!is.na(tGammaprep)] +tBetaXZiprep[!is.na(tBetaXZiprep)]
  }else{
    tooff = tGammaprep[!is.na(tGammaprep)]
  }
  #calculation of gradient, M0 refers as.vector(thetai%*%Afori) + tooff
  eUi  = as.vector(thetai%*%tAfori)+tooff
  XiM0=xi.func(eUi)
  #grad =  rowSums(matrix(yi- XiM0,ncol = length(yi),nrow=length(thetai) ,byrow =TRUE)*Afori)/J

  grad<- rcpp_cal_grad(yi,XiM0, t(tAfori), J,K)
  #Hessian
  Hes<-rcpp_cal_Hes(as.vector(XiM0), t(tAfori), thetai, J, length(XiM0))
  drt <- get_drt.func(Hes,grad,K)
  #--------------------------------------------------------------------------------------------
  #Update parameter by line search
  out = rcpp_getnewpar_func(XiM0,yi,J,grad,drt, thetai,t(tAfori ),tooff)
  return(out)
}

#Function to normalise parameters according to conditions
norm.func = function(J,N,cov_len,px,tA,Theta,Gamma,tBeta,X){
  # Check for NA values
  if(any(is.na(tA)) || any(is.na(Theta))){
    print("NA detected, cannot do normalisation")
    tAstar = matrix(NaN, nrow= length(tA[,1]), ncol = length(tA[1,]))
    Thetastar = matrix(NaN, nrow= length(Theta[,1]), ncol = length(Theta[1,]))
  }else{
    IN = matrix(1, nrow =N, ncol =1)
    if(cov_len>0){
      INZi = cbind(IN,X)
    }else{
      INZi = IN
    }

    #Number of observed time points
    Tp =ncol(Gamma)
    #TILTH: Tilde THETA
    TILTH = Theta - INZi%*%solve(t(INZi)%*%INZi)%*%t(INZi)%*%Theta
    #Projection of Theta in INZ times tA
    Q = as.matrix(solve(t(INZi)%*%INZi)%*%t(INZi)%*%Theta%*%tA)
    #Get Betastar and Cstar
    Gammastar = Gamma + matrix(Q[1,], nrow=J, ncol=Tp)
    if(cov_len>0){
      tBetastar = tBeta
      if(px>0){
        tBetastar[1:px,] = tBetastar[1:px,] +Q[2:(px+1),]
      }
    }else{
      tBetastar=matrix(0, nrow=0, ncol=0)
    }

    #Get Astar and Thetastar
    SigmaJA = tA%*%t(tA)/J
    SigmaNT = t(TILTH)%*%TILTH/N
    object = eigen(expm::sqrtm(SigmaJA)%*%SigmaNT%*%expm::sqrtm(SigmaJA))
    HNJ = t(solve(expm::sqrtm(SigmaJA))%*%object$vectors)
    tAstar = HNJ%*%tA #tAstar%*%t(tAstar)/J Gives Identity matrix
    Thetastar= TILTH%*%solve(HNJ)#t(Thetastar)%*%Thetastar/N gives diagonal matrix
  }
  return(list( "tAstar" = tAstar, "Thetastar"= Thetastar, "Gammastar"=Gammastar, "tBetastar"= tBetastar))
}


#projection function
proj.func = function(out,K){
  constant=5*sqrt(K)
  value = sqrt(sum(out^2))
  if ( value >constant){
    out=out*constant/value
  }
  return(out)
}

#Function to get XZmat
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
GetYmat.func = function(Y,N,J,Tp){
  Ymat <- matrix(nrow = N, ncol = J * Tp)
  # Fill the matrix using a loop over J
  for (j in 1:J){
    Ymat[, ((j-1)*Tp+1):(j*Tp)] <- Y[,j,]
  }
  return(Ymat)
}


#Function to get descend direction from Hessian matrix and gradient
get_drt.func = function(Hes,grad,par_len){
  if (det(Hes)!=0){
    drt = as.vector(-solve(Hes)%*%grad)
  }else{
    drt = as.vector(-solve(Hes - 0.1*diag(par_len))%*%grad )
  }
  return(drt)
}

#inverse function of the logistic function
invxi.func = function(p){
  out = log(p/(1-p))
  return(out)
}


#This script contains the function to generate data for simulation
#Function to combine Gamma,Beta, A to form U
GetU.func = function(Gamma,Beta,A){
  if(nrow(Beta)>0){
    return(cbind(Gamma,Beta,A))
  }else{
    return(cbind(Gamma,A))
  }

}

#Function to extract Gamma given U
GetGammafromU.func = function(U,Tp){
  return(U[,1:Tp])
}


#Function to extract Beta given U
GetBetafromU.func = function(U,Tp,cov_len ){
  if(cov_len >0){
    return(U[,(Tp+1):(Tp+cov_len)])
  }else{
    return(matrix(0, nrow=0, ncol=0))
  }
}

#Function to extract A given U
GetAfromU.func = function(U,Tp,p,Kstar){
  return(U[,(Tp+p+1):(Tp+p+Kstar)])
}

#Function to combine X and Z
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
xi.func = function(x){
  expx = exp(x)
  out=expx/(1 + expx)
  #Make sure it does not explode
  out[x>700]=1
  out[x< -800]=0
  return(out)
}

#-----------------------------------------------------------------------------------------------------------------
#Functions for calculating asymptotic variance
#Return a array that gives uj^t e_{ir} for every i and j
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

#Function for computing Phi_j
calPhi_j.func = function(EEstarstore,ujE,R,N,J,Tp,ujlen){
  rhoipp<- -exp(ujE)/(1+exp(ujE) )^2
  #Expectation of varrho'' evaluated using ujE values
  Evarrhoipp <- matrix(0, nrow = N, ncol =Tp)
  for ( t in 1:Tp){
    Evarrhoipp[,t]<-R[,t]* rhoipp[,t]
  }

  #Compute Phi_i. Tried vectorising but seems looping is the fastest. May consider using Rcpp if too slow.
  Phi_j <- matrix(0, nrow=ujlen, ncol=ujlen)
  for ( t in 1:Tp){
    for ( i in 1:N){
      #To be updated
      Phi_j <- Phi_j + Evarrhoipp[i,t]*EEstarstore[,,i,t]
    }
  }

  Phi_j <- Phi_j/N
  return(Phi_j)
}

calEEstore.func = function(E, ujlen, N,Tp){
  EEstore = array(0,dim = c(ujlen, ujlen,N,Tp))
  for ( i in 1:N){
    for ( r in 1:Tp){
      EEstore[,,i,r]<- E[i,,r]%*%t(E[i,,r])
    }
  }
  return(EEstore)
}




#Function to compute asymptotic variance of u_j.
calasympvar.func = function(X,Z, Thetahat, Uhat,R,J,N,Tp){
  ujlen= ncol(Uhat)
  #Recover E from ujlen, R, Thetahat, Z.
  Ehat<-array(0,dim=c(N,ujlen,Tp))
  px = ncol(X)
  pz =dim(Z)[2]
  cov_len = px+pz
  XZ = CombXZ.func(X,Z,px,pz,Tp,cov_len,N)
  for ( t in 1:Tp){
    Dprep = matrix(0, nrow=N,ncol=Tp)
    Dprep[,t]=1
    if(cov_len>0){
      Ehat[,,t]<- cbind(Dprep,XZ[,,t],Thetahat )
    }else{
      Ehat[,,t]<- cbind(Dprep, Thetahat)
    }
  }
  UERhat<- UER.func(Uhat, Thetahat,XZ,Tp,cov_len,J,N )

  EEhatstore = calEEstore.func(Ehat,ujlen,N,Tp)
  #Now really calculate the estimate of asymptotic variance for each item
  AsymVarhatstore = array(0, dim =c(ujlen,ujlen,J))
  CILOstore = matrix(0, nrow=ujlen, ncol =J)
  CIUPstore = matrix(0, nrow=ujlen, ncol =J)
  for ( j in 1: J){
    ujEhat <- UERhat[,j,]
    #Compute Phi_i_hat
    Phi_j_hat <-calPhi_j.func(EEhatstore,ujEhat,R,N,J,Tp,ujlen)
    #---------------------------------------------------------
    #Now compute asymptotic variance
    AsymVarhat<- solve(-Phi_j_hat)
    AsymVarhatstore[,,j]<-AsymVarhat
  }
  return(AsymVarhatstore)
}





