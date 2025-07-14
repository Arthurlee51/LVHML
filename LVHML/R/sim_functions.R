#Functions used for simulations:
#Function to calculate Frobenius norm
Fnormsum.func = function(a,b){
  return(sqrt(sum((a - b)^2 )))
}

#Functions to evaluate performance of the estimator.
eval.func =function(Thetahat,Uhat,Thetastar,Ustar,X,Z,N,J,Tp,px,pz,K, ext, gamma_fix)
{
  #Precompute
  IN = matrix(1, nrow=N,ncol=1)
  g_len <- if (gamma_fix) 1 else Tp
  b_len <- if(ext) Tp*px else px
  Totallossstore = rep(0,Tp)
  for ( t in 1:Tp){
    if(!gamma_fix){
      tosum = IN%*%t(Uhat[,t])  - IN%*%t(Ustar[,t])
    }else{
      tosum = t*IN%*%t(Uhat[,1])  - t*IN%*%t(Ustar[,1])
    }

    if(!ext){
     tosum <- tosum + Thetahat%*%t(Uhat[, (g_len+b_len+pz+1):ncol(Uhat)]) - Thetastar%*%t(Ustar[,(g_len+b_len+pz+1):ncol(Uhat) ])
    }else{
      tosum <- tosum + Thetahat%*%t(Uhat[, (g_len+b_len+pz+(t-1)*K+1):(g_len+b_len+pz+t*K)]) - Thetastar%*%t(Ustar[,(g_len+b_len+pz+(t-1)*K+1):(g_len+b_len+pz+t*K) ])
    }

    if(px>0){
      if(!ext){
        tosum <- tosum + X%*%( t(Uhat[,(g_len+1):(g_len+px)]) - t(Ustar[,(g_len+1):(g_len+px)]) )
      }else{
        tosum <- tosum + X%*%( t(Uhat[,(g_len+(t-1)*px+1):(g_len+t*px)]) - t(Ustar[,(g_len+(t-1)*px+1):(g_len+t*px)]) )
      }
    }

    if(pz>0){
      tosum <- tosum + Z[,,t]%*%( t(Uhat[,(g_len+b_len+1):(g_len+b_len+pz)]) - t(Ustar[,(g_len+b_len+1):(g_len+b_len+pz)]) )
    }
    Totallossstore[t] = sqrt(sum(tosum^2)/(N*J))
  }
  Totalloss = max(Totallossstore)
  A_range <- (g_len+b_len+pz+1):ncol(Uhat)
  Ahat = as.matrix(Uhat[,A_range])
  Astar = as.matrix(Ustar[,A_range])

  #Compute sign matrix
  if(!ext){
    if(K==1){
      ShatA=sign(t(Astar)%*%Ahat/J)
    }else{
      ShatA=sign(diag(diag(t(Astar)%*%Ahat/J)))
    }

  }else{
    if(K==1){
      ShatA=sign(t(Astar[,1:K])%*%Ahat[,1:K]/J)
    }else{
      ShatA=sign(diag(diag(t(Astar[,1:K])%*%Ahat[,1:K]/J)))
    }

  }


  Thetaloss = Fnormsum.func(Thetahat,Thetastar%*%ShatA)/sqrt(N)
  Ulossstore = rep(0,Tp)
  for ( r in 1:Tp){
    g_pos <- if(gamma_fix) 1 else r

    if(!ext){
      Ut_range <- c(g_pos, (g_len+1):ncol(Uhat))
    }else{
      Ut_range <- c(g_pos)
      if(px>0){
        Ut_range <- c(Ut_range,(g_len + (r-1)*px+1): (g_len + r*px))
      }
      if(pz>0){
        Ut_range <- c(Ut_range,(g_len + Tp*px+1): (g_len + Tp*px+pz))
      }
      Ut_range <- c(Ut_range,(g_len + Tp*px+pz +(r-1)*K+ 1): (g_len + Tp*px+pz +r*K ) )
    }
    Urhat = Uhat[,Ut_range]
    Urstar = Ustar[,Ut_range]
    SUhat = diag(c(rep(1, 1+px+pz), diag(ShatA)) )
    Ulossstore[r]  =  Fnormsum.func(Urhat,Urstar%*%SUhat)/sqrt(J)
  }
  Uloss =max( Ulossstore)
  return(list(Ulossstore = Ulossstore,Totalloss = Totalloss,Totallossstore=Totallossstore,  Uloss = Uloss, Thetaloss = Thetaloss ))
}


#Function to compute confidence interval based on parameters estimate and data:
calCI_uj.func = function(X,Z, Thetahat, Uhat,R,J,N,Tp,ext, gamma_fix){
  ujlen= ncol(Uhat)
  #Recover Ehat from ujlen, R, Thetahat, Z.
  px = ncol(X)
  pz =dim(Z)[2]
  K <- ncol(Thetahat)
  Ehat <- compute_Ehat( Thetahat, X, Z,  ujlen,K,Tp, px,pz, N,ext, gamma_fix)
  UERhat<- UER.func(Uhat, Thetahat, X,Z, Tp, px,pz, J, N,ext,gamma_fix)

  EEhatstore = calEEstore.func(Ehat,ujlen,N,Tp)
  AsymVarhatstore = array(0, dim =c(ujlen,ujlen,J))
  CILOstore = matrix(0, nrow=ujlen, ncol =J)
  CIUPstore = matrix(0, nrow=ujlen, ncol =J)
  for ( j in 1: J){
    ujEhat <- UERhat[,j,]
    #Compute Phi_i_hat
    Phi_j_hat <-calPhi_j.func(EEhatstore,ujEhat,R,N,J,Tp,ujlen)
    #---------------------------------------------------------
    #Now compute asymptotic variance
    # Try solving the -Phi_j_hat directly
    AsymVarhat = tryCatch({
      solve(-Phi_j_hat)
    }, error = function(e) {
      NULL
    })

    # If direct inversion failed, regularize the Phi_j_hat
    if (is.null(AsymVarhat)) {
      message("Matrix inversion failed; regularizing Phi_j_hat.")
      Phi_j_hat_reg = Phi_j_hat - 0.01 * diag(ujlen)
      AsymVarhat = solve(-Phi_j_hat_reg)
    }

    AsymVarhatstore[,,j]<-AsymVarhat
    #Compute confidence intervals
    CILOstore[,j] = Uhat[j,] - 1.96*sqrt(diag(AsymVarhat))/sqrt(N)
    CIUPstore[,j] = Uhat[j,] + 1.96*sqrt(diag(AsymVarhat))/sqrt(N)
  }
  return(list("AsymVarhatstore"=AsymVarhatstore, "CILOstore"=CILOstore,"CIUPstore"=CIUPstore) )
}


#Function to compute confidence interval for theta_i based on paramters estimate and data:
calCI_thi.func = function(X,Z, Thetahat,Uhat, Ahat,R,J,N,Tp,ext, gamma_fix){
  K <- ncol(Thetahat)
  px = ncol(X)
  pz =dim(Z)[2]

  UERhat<- UER.func(Uhat, Thetahat, X,Z, Tp, px,pz, J, N,ext,gamma_fix)

  #Rescale Ahat to suit the requirement of calEEstore.
  if(!ext){
    Ahat <- array(Ahat, dim = c(J, K, Tp))
  }else{
    old_Ahat <- Ahat
    Ahat <- array(0, dim = c(J, K, Tp))
    for(t in 1:Tp){
     Ahat[,,t] <-  old_Ahat[,( (t-1)*K + 1):(t*K) ]
    }
  }

  AAhatstore = calEEstore.func(Ahat,K,J,Tp)
  AsymVarhatstore = array(0, dim =c(K,K,N))
  CILOstore = matrix(0, nrow=K, ncol =N)
  CIUPstore = matrix(0, nrow=K, ncol =N)
  for ( i in 1:N){
    Uthihat <- UERhat[i,,]
    #Compute Phi_i_hat
    Psi_i_hat <-calPsi_i.func(AAhatstore,Uthihat,R[i,],N,J,Tp,K)
    #---------------------------------------------------------
    #Compute asymptotic variance
    AsymVarhat<- solve(-Psi_i_hat)
    AsymVarhatstore[,,i]<-AsymVarhat
    #Compute confidence intervals
    CILOstore[,i] = Thetahat[i,] - 1.96*sqrt(diag(AsymVarhat))/sqrt(J)
    CIUPstore[,i] = Thetahat[i,] + 1.96*sqrt(diag(AsymVarhat))/sqrt(J)
  }
  return(list("AsymVarhatstore"=AsymVarhatstore, "CILOstore"=CILOstore,"CIUPstore"=CIUPstore) )
}



#Compute Psi_i
calPsi_i.func <- function(AAhatstore, Uthihat, ri, N, J, Tp, K) {
  rhoipp <- -exp(Uthihat) / (1 + exp(Uthihat))^2
  Evarrhoipp <- matrix(0, nrow = J, ncol = Tp)

  for (t in 1:Tp) {
    Evarrhoipp[, t] <- ri[t] * rhoipp[, t]
  }

  Psi_i <- matrix(0, nrow = K, ncol = K)
  for (t in 1:Tp) {
    for (j in 1:J) {
      Psi_i <- Psi_i + Evarrhoipp[j, t] * AAhatstore[,, j, t]
    }
  }

  Psi_i <- Psi_i / J
  return(Psi_i)
}


#Functions to run glm beta. Serve comparison in simulation only. Corresponds to LR in manuscript.
estglmbeta.func = function(Y, R, X = matrix(0, nrow = 0, ncol = 0), Tp, par = 0, n.cores = 1, ext, gamma_fix){
  # Get dimensions of input matrices
  px = ncol(X)
  N = dim(Y)[1]
  J = dim(Y)[2]

  # Convert Y to logical matrix to reduce memory usage
  Ymatlog <- GetYmat.func(Y, N, J, Tp) == 1

  indices_R <- unlist(apply(R, 1, function(row) which(row == 1)))  # Indices of timepoints where r_{it} = 1
  lengths_R <- rowSums(R)  # Lengths of the timepoints

  if(par ==0){
    Uhat= sapply(1:J, function(j) as.vector( Betaonlyparglm.func(t(Y[,j,]),indices_R,lengths_R,X,N,px,K,Tp,ext, gamma_fix) ) )
  }
  if(par ==1){
    Uhat =mcmapply(function(j) as.vector( Betaonlyparglm.func(t(Y[,j,]),indices_R,lengths_R,X,N,px,K,Tp,ext, gamma_fix) ) ,1:J,mc.cores = n.cores )
  }
  #Return output
  return(Uhat)
}

#Function to run glm function to get estimate of beta.
#yris = t(Y[,j,])
Betaonlyparglm.func = function(yris,indices_R,lengths_R,X,N,px,K,Tp, ext, gamma_fix){
  yj = yris[!is.na(yris)]

  id = rep(1:N, lengths_R)

  Eforj =rcpp_BetaonlygetEforj( indices_R-1,lengths_R ,  X, N, Tp,  px,ext, gamma_fix)$Eforj
    glmfit <- glm.fit(Eforj,yj, family = binomial())
    coef = glmfit$coefficients

  return(coef)
}


#Functions to run glmm(Generalised linear mixed models. In particular, logistic regression with random intercepts.) Corresponds to LRRI in manuscript
estglmm.func = function(Y, R, X = matrix(0, nrow = 0, ncol = 0), Tp, par = 0, n.cores = 1,ext, gamma_fix){
  # Get dimensions of input matrices
  px = ncol(X)
  N = dim(Y)[1]
  J = dim(Y)[2]

  indices_R <- unlist(apply(R, 1, function(row) which(row == 1)))  # Indices of timepoints where r_{it} = 1
  lengths_R <- rowSums(R)  # Lengths of the timepoints

  if(par ==0){
    Uhat= sapply(1:J, function(j) as.vector(parglmm.func(t(Y[,j,]),indices_R,lengths_R,X,N,px,K,Tp,ext, gamma_fix) ) )
  }
  if(par ==1){
    Uhat =mcmapply(function(j) as.vector(parglmm.func(t(Y[,j,]),indices_R,lengths_R,X,N,px,K,Tp,ext, gamma_fix) ) ,1:J,mc.cores = n.cores )
  }
  #Return output
  return(Uhat)
}

#Function to run glmm for each given j
#yris = t(Y[,j,])
parglmm.func = function(yris,indices_R,lengths_R,XZmat,N,P,K,Tp,ext, gamma_fix){
  yj = yris[!is.na(yris)]

  id = rep(1:N, lengths_R)

  Eforj =rcpp_BetaonlygetEforj( indices_R-1,lengths_R ,  X, N, Tp,  px,ext, gamma_fix)$Eforj

  #Data for the relevant part.
  data_par<- data.frame(y=yj,id=id,E = Eforj, indices_R = indices_R )
  g_len <- if(gamma_fix) 1 else Tp
  b_len <- if(ext) Tp*px else px

  colnames(data_par)[3:(2+g_len)] <-paste("gamma", 1:g_len, sep = "")

  #Set colnames for future analysis
  colnames(data_par)[(2+g_len+1):(2+g_len+b_len)]<- paste("x", 1:b_len, sep = "")

  gamma_terms <- paste("gamma", 1:g_len, sep = "")
  x_terms <- paste("x", 1:b_len, sep = "")

  # Combine into one string with the + operator between each term
  fixed_effects <- paste(c(gamma_terms, x_terms), collapse = " + ")
  # Build the entire formula string with no intercept (0 + ...) and random effect(s).
  if(!ext){
    formula_str <- paste("y ~ 0 +", fixed_effects, "+ (1 | id)")
  }else{
    formula_str <- paste("y ~ 0 +", fixed_effects, "+ (1 | id) + (1 | indices_R)")
  }

  glmfit <- glmer(formula_str, data = data_par, family = binomial)
  coef <- as.vector(fixef(  glmfit))

  return(coef)
}

