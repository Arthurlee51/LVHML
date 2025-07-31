#Test performance of FDR control and asymptotic normality.
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
  #Get lengths of parameters.
  g_len <- if(gamma_fix) 1 else Tp
  b_len <- if(ext) Tp*px else px
  A_len <- if(ext) Tp*Kstar else Kstar
  u_len <- g_len + b_len + A_len
  
  for(mul in c(5,10)){
    Prop = matrix(0, nrow = length(Jlong), ncol = u_len )
    #summary statistics
    sum_stat = array(0, dim = c(6,u_len , length(Jlong)))
    Propstore_list = list()
    avg_AsyVarhatstore <- list()
    mean_ECPs = rep(0, length(Jlong))
for( l in 1:length(Jlong)){
J = Jlong[l]
N = mul*J
nsim =100
load(sprintf("Z%dfixg%dext%dsimJ%dN%dKstar%dnsim%d.Rdata",withZ*1,gamma_fix*1,ext*1, J,N,Kstar,nsim) )
tUstar = t(out$Ustar)
B_range <- (g_len+1):(g_len+b_len)
Betahatstore = out$Uhatstore[,B_range,]
#Also get Betastar for later FDR computation.
tBetastar = tUstar[B_range,]
#Get indices where the null hypotheses are TRUE. If TRUE =1, FALSE = 0
indices_NULL <- matrix(0,nrow =3, ncol = J)
indices_NULL[1,] <- (tBetastar[1,]==0)*(tBetastar[2,]==0)
indices_NULL[2,] <- (tBetastar[3,]==0)*(tBetastar[4,]==0)
indices_NULL[3,] <- (tBetastar[5,]==0)
#Calculate asymptotic variance and store key info in each simulation for comparison.
#Store proportion of containing true parameters
Containstore = array(0, dim = c(u_len, J, nsim))
AsymVarhatstorestore <- array(0, dim =c( u_len,u_len, J, nsim) )
FDRstore <- matrix(0, nrow=3, ncol = nsim)#Three sets of FDR for the three variables (X1(two levels), X2(two levels) and X3(continuous))
FNDRstore <-matrix(0, nrow=3, ncol = nsim) #Storing False non-discovery rate, a power measure.
powerstore <-matrix(0, nrow=3, ncol = nsim) #Storing False non-discovery rate, a power measure.
#Total number of rejected hypothesis.
Total_rej <- matrix(0, nrow=3, ncol = nsim)
for (z in 1:nsim){
  print(z)
  #Calculate sign matrices. Equivalent to sign of A_1 if ext = TRUE and A otherwise. 
  A1ind <-   (g_len + b_len +1):(g_len + b_len +Kstar)    
  SAdiag = diag(sign(tUstar[A1ind, ]%*%out$Uhatstore[,A1ind,z]))
  SU <- diag(c( rep(1, g_len + b_len ) , rep(SAdiag, times = ifelse(ext, Tp, 1)) ) )
  #Compute CI for u_j
  test <- LVHML:::calCI_uj.func(out$X,out$Z, out$Thetahatstore[,,z] , out$Uhatstore[,,z],out$Rstore[,,z] ,J,N,Tp,ext, gamma_fix)
  Containstore[,,z] <- (test$CILOstore <= SU%*%tUstar)*(SU%*%tUstar <= test$CIUPstore)
  #Propstore[,z]<- rowMeans(Contain )
  AsymVarhatstorestore[,,,z] <- test$AsymVarhatstore
  #Perform Wald's test using the asymptotic variance computed in test.
  Beta_var <- test$AsymVarhatstore[B_range ,B_range ,]
  #Only works on Wald's test for the main model. #Test against the hypotheses that the coefficients are 0.
  if(!ext & !gamma_fix){
  Test_stat_wald<- matrix(0, nrow = J, ncol = 3)
    Betahat <- Betahatstore[,,z]
    for ( j in 1:J){
      Beta_var_j <- Beta_var[,,j]
      #Test statistic for beta1 and beta2 
      Test_stat_wald[j,1]<- Betahat[j,1:2]%*%solve(Beta_var_j[1:2,1:2]/N)%*%Betahat[j,1:2]
      #Test statistic for beta3 and beta4 
      Test_stat_wald[j,2]<- Betahat[j,3:4]%*%solve(Beta_var_j[3:4,3:4]/N)%*%Betahat[j,3:4]
      #Test statistic for beta5
      Test_stat_wald[j,3]<- Betahat[j,5]%*%solve(Beta_var_j[5,5]/N)%*%Betahat[j,5]
      #Test statistics for beta_1, 2 and 3 compared with true values 
    }
    #Get p-values by comparing to the relevant chi^2 distribution, calculate them separetely. 
    p_values <- matrix(0, nrow = J, ncol = 3)
    p_values[,1]<- 1- pchisq(Test_stat_wald[,1] ,df=2)
    p_values[,2]<- 1- pchisq(Test_stat_wald[,2] ,df=2)
    p_values[,3]<- 1- pchisq(Test_stat_wald[,3] ,df=1)
    #Adjust each of the p-values using FDR, BY procedure
    #Adjust p_valus separately. 
    adjusted_p_values <- matrix(0, nrow = J, ncol = 3)
    adjusted_p_values[,1] <- p.adjust(p_values[,1], method = "BY") 
    adjusted_p_values[,2] <- p.adjust(p_values[,2], method = "BY") 
    adjusted_p_values[,3] <- p.adjust(p_values[,3], method = "BY") 
    #Get Significance matrix
    alpha <- 0.05 #alpha:FDR control threshold. 
    Sig_matrix  <- (adjusted_p_values <= alpha)
    Total_rej[1:3 ,z] <- colSums(Sig_matrix)
    #Get indicies of false discovery/number of rejections.
    FDR <- rep(0, 3)
    FNDR <- rep(0, 3)
    power <- rep(0,3)
    for(k in 1:3){
      Sig_k <- Sig_matrix[,k]
      FDR[k] <- sum(Sig_k[which(indices_NULL[k,]==1)] )/Total_rej[k,z]
      FNDR[k] <-sum(1-Sig_k[which(indices_NULL[k,]==0)] ) /(J - Total_rej[k,z])
      power[k] <- sum(Sig_k[which(indices_NULL[k,]==0)])/length(which(indices_NULL[k,]==0))
    }
    FDR[colSums(Sig_matrix)==0] <- 0
    FDRstore[1:3 ,z] <-  FDR
    FNDRstore[1:3,z] <- FNDR
    powerstore[1:3,z] <- power
  }
}
Propstore_list[[l]] = apply( Containstore,c(1,2), mean )
sum_stat[,,l] <- summary(t(Propstore_list[[l]]))
mean_ECPs[l] <- mean(t(Propstore_list[[l]])[,(g_len+1):(g_len + b_len)])
Beta_ECPs <- Propstore_list[[l]][(g_len+1):(g_len + b_len),]
#Prepare data summary to be saved
if(!ext & !gamma_fix){
  mean_FDR <- rowSums(FDRstore)/nsim
  mean_Total_rej <-  rowSums(Total_rej)/nsim
  mean_FNDR <- rowSums(FNDRstore)/nsim
  mean_power <- rowSums(powerstore )/nsim
  FDR_result <- data.frame(cbind(mean_FDR,mean_Total_rej,mean_FNDR,mean_power))
  max_FDR <- max(mean_FDR)
  save(Beta_ECPs, Tp,px,pz,  FDR_result,max_FDR,file = sprintf("Z%dfixg%dext%dECPandFDRstoreJ%dN%dKstar%d.Rdata", withZ*1, gamma_fix*1,ext*1, J,N,Kstar) )
}else{
 print("Not applicable")
}


}
}
}

