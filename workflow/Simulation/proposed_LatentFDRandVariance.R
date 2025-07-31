#Script to test the performance of Asymp variance and ECP(ignore FDR) for latent variables(A and Theta)
#Test performance of plug-in estimator using simulation results, as well as performance of BH procedure in controlling FDR.
library(dplyr)
library(tidyr)
library(LVHML)
#Source function for computing asymptotic variance first. 
source("functions_for_experiment.R")

withZ <- FALSE
ext <- FALSE
gamma_fix <- FALSE
#Load data file.
Jlong = c(100,200,300,400)
for(Kstar in c(3,8)){
  #Recall Tp from simulation setting.
  Tp <- 4
  px <- 5
  pz <- 0
  if(!gamma_fix){
    g_len <- Tp 
  }else{
    g_len <- 1
  }
  
  if(!ext){
    b_len <- px
    A_len <- Kstar
  }else{
    b_len <- Tp*px
    A_len <- Tp*Kstar
  }
  u_len <- g_len + b_len + A_len
  
  for(mul in c(10)){
    Prop = matrix(0, nrow = length(Jlong), ncol = u_len )
    #summary statistics
    sum_stat = array(0, dim = c(6,u_len , length(Jlong)))
    Propstore_list = list()
    Propstore_list_T = list()
    avg_AsyVarhatstore <- list()
    mean_A_ECPs = rep(0, length(Jlong))
    #Mean of Theta's ECPS
    mean_T_ECPs = rep(0, length(Jlong))
    for( l in 1:length(Jlong)){
      J = Jlong[l]
      N = mul*J
      nsim =100
      load(sprintf("gapZ%dfixg%dext%dsimJ%dN%dKstar%dnsim%d.Rdata",withZ*1,gamma_fix*1,ext*1, J,N,Kstar,nsim) )
     
      tUstar = t(out$Ustar)
      Thetastar <- out$Thetastar
      Containstore <- array(0, dim = c(u_len, J, nsim))
      #Containstore for Theta
      Containstore_T <- array(0, dim =c(Kstar,N,nsim))
      AsymVarhatstorestore <- array(0, dim =c( u_len,u_len, J, nsim) )
      AsymVarhatstorestore_T <- array(0, dim =c( Kstar,Kstar, N, nsim) )
      signs_mat <- matrix(0, nrow = nsim, ncol = Kstar)
      for (z in 1:nsim){
        print(z)
        #Calculate sign matrices. Use A_1 if ext = TRUE and A otherwise. 
        A1ind <-   (g_len + b_len +1):(g_len + b_len +Kstar)    
        SAdiag = diag(sign(tUstar[A1ind, ]%*%out$Uhatstore[,A1ind,z]))
        signs_mat[z,] <-  SAdiag 
        SU <- diag(c( rep(1, g_len + b_len ) , rep(SAdiag, times = ifelse(ext, Tp, 1)) ) )
        Ahat <- out$Uhatstore[,(g_len + b_len +1):u_len,z]
        #Compute CI for u_j
        test <- LVHML:::calCI_uj.func(out$X,out$Z, out$Thetahatstore[,,z] , out$Uhatstore[,,z],out$Rstore[,,z] ,J,N,Tp,ext, gamma_fix)
        Thetahat <- as.matrix(out$Thetahatstore[,,z])
        test_Theta <- LVHML:::calCI_thi.func(out$X,out$Z, Thetahat,out$Uhatstore[,,z], Ahat,out$Rstore[,,z],J,N,Tp,ext, gamma_fix)
        Containstore[,,z] <- (test$CILOstore <= SU%*%tUstar)*(SU%*%tUstar <= test$CIUPstore)
        if(Kstar==1){
          SA <- diag(as.matrix(SAdiag) )
        }else{
          SA <- diag(SAdiag )
        }
        
        Containstore_T[,,z] <- (test_Theta$CILOstore <= SA%*%t(Thetastar))*(SA%*%t(Thetastar) <= test_Theta$CIUPstore)
      }
      
      Propstore_list[[l]] = apply( Containstore,c(1,2), mean )
      Propstore_list_T[[l]] = apply(Containstore_T,c(1,2),mean)
      mean_A_ECPs[l] <- mean(t(Propstore_list[[l]])[,(g_len+b_len+1):u_len])
      A_ECPs <- Propstore_list[[l]][(g_len+b_len+1):u_len,]
      mean_T_ECPs[l] <- mean(Propstore_list_T[[l]])
      T_ECPs <- Propstore_list_T[[l]]
      if(px==5){
        save(Containstore,Containstore_T,signs_mat,mean_A_ECPs,A_ECPs,mean_T_ECPs,T_ECPs,file = sprintf("Latent_Z%dfixg%dext%dECPstoreJ%dN%dKstar%d.Rdata", withZ*1, gamma_fix*1,ext*1, J,N,Kstar) )
      }
    }
  }
}

