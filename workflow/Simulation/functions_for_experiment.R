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
  #Trial to make a "diagonal sign matrix"
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



#Functions to run glm beta. Serve comparison in simulation only. Corresponds to LR in manuscript.
estglmbeta.func = function(Y, R, X = matrix(0, nrow = 0, ncol = 0), Tp, par = 0, n.cores = 1, ext, gamma_fix){
  # Get dimensions of input matrices
  px = ncol(X)
  N = dim(Y)[1]
  J = dim(Y)[2]
  
  # Convert Y to logical matrix to reduce memory usage
  Ymatlog <- LVHML:::GetYmat.func(Y, N, J, Tp) == 1
  
  indices_R <- unlist(apply(R, 1, function(row) which(row == 1)))  # Indices of timepoints where r_{it} = 1
  lengths_R <- rowSums(R)  # Lengths of the timepoints
  
  #Todo: Function to perform estimation
  if(par ==0){
    Uhat= sapply(1:J, function(j) as.vector( LVHML:::Betaonlyparglm.func(t(Y[,j,]),indices_R,lengths_R,X,N,px,K,Tp,ext, gamma_fix) ) )
  }
  if(par ==1){
    Uhat =mcmapply(function(j) as.vector( LVHML:::Betaonlyparglm.func(t(Y[,j,]),indices_R,lengths_R,X,N,px,K,Tp,ext, gamma_fix) ) ,1:J,mc.cores = n.cores )
  }
  #Return output
  return(Uhat)
}

#Function to run glm function to get estimate of beta. 
#yris = t(Y[,j,])
Betaonlyparglm.func = function(yris,indices_R,lengths_R,X,N,px,K,Tp, ext, gamma_fix){
  yj = yris[!is.na(yris)]
  
  id = rep(1:N, lengths_R)
  
  Eforj =LVHML:::rcpp_BetaonlygetEforj( indices_R-1,lengths_R ,  X, N, Tp,  px,ext, gamma_fix)$Eforj
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
  
  #Todo: Function to perform estimation
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
  
  Eforj =LVHML:::rcpp_BetaonlygetEforj( indices_R-1,lengths_R ,  X, N, Tp,  px,ext, gamma_fix)$Eforj
  
  #Data for the relevant part.
  data_par<- data.frame(y=yj,id=id,E = Eforj, indices_R = indices_R )
  g_len <- if(gamma_fix) 1 else Tp
  b_len <- if(ext) Tp*px else px
  
  colnames(data_par)[3:(2+g_len)] <-paste("gamma", 1:g_len, sep = "")
  #data_par<- data.frame(y=yj,id=id,time=factor(indices_R),X =Eforj[,(Tp+1):ncol(Eforj)] )
  
  #Set colnames for future analysis
  colnames(data_par)[(2+g_len+1):(2+g_len+b_len)]<- paste("x", 1:b_len, sep = "")
  #c("x1","x2","x3","x4")
  
  # Assume g_len and b_len are defined
  gamma_terms <- paste("gamma", 1:g_len, sep = "")
  x_terms <- paste("x", 1:b_len, sep = "")
  
  # Combine into one string with the + operator between each term
  fixed_effects <- paste(c(gamma_terms, x_terms), collapse = " + ")
  # Build the entire formula string with no intercept (0 + ...) and random effect(s).
  if(!ext){
    formula_str <- paste("y ~ 0 +", fixed_effects, "+ (1 | id)")
  }else{
    print("Not applicable for this extension")
  }
  
  glmfit <- glmer(formula_str, data = data_par, family = binomial)
  coef <- as.vector(fixef(  glmfit))
  
  return(coef)
}
