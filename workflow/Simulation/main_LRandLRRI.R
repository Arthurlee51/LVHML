# Main script to run simulations. Run both LR and LRRI results as shown in the report.
# --------------------------------------------------------------------------------

# Load required libraries
library(BB)
library(expm)
library(doParallel)
library(Rcpp)
library(truncnorm)
library(glmnet)
library(lme4)

# # Source functions from external R scripts
source("functions_for_experiment.R")
# --------------------------------------------------------------------------------
# Choosing N relative to J (options: "s" for small and "L" for large)
#nset = "s"  
Tp <- 4
withZ = FALSE#False if we do not have Z, as in simulation study.
ext = FALSE #False if we have A and Beta both time-independent, as in simulation study.
gamma_fix = FALSE#False if we allow gamma_t being time-dependent. TRUE when gamma_t = gamma for all t. Default is FALSE.
# Run simulations for different values of J
for (J in c(100,200,300,400)) {
  for(nset in c("s","L")){
    set.seed(123)
    #Generate N
    if(nset =="s"){
      N = 5*J
    }
    if(nset =="L"){
      N = 10*J
    }
    for(Kstar in c(3)){
      
      par = 0  # 0 for normal, 1 for parallel computing
      n.cores = 1  # Number of cores for parallel computing
      nsim = 100  # Number of simulations
      
      
      # Generate parameters based on the settings
      parameters <- LVHML:::genpar.func(Kstar,N,J,Tp, withZ, ext, gamma_fix)
      px <- parameters$px
      pz <- parameters$pz
      
      
      # Initialize storage arrays and variables
      Thetastar = parameters$Thetastar
      Ustar = parameters$Ustar
      Z = parameters$Z
      X = parameters$X
      #Get the length of gamma and beta based on the model assumption.
      if(!gamma_fix){
        g_len <- Tp
      }else{
        g_len <- 1
      }
      if(!ext){
        b_len <- px
      }else{
        b_len <- px*Tp
      }
      Bstar <- Ustar[,(g_len+1):(g_len+b_len)]
      Uhatstore_LR = array(0, dim = c(J, g_len + b_len , nsim))
      Uhatstore_LRRI = array(0, dim = c(J, g_len + b_len , nsim))
      Rstore = array(0, dim = c(N, Tp, nsim))
      niterstore = rep(0, nsim)
      Improvementsobjstore = rep(0, nsim)
      objstore = rep(0, nsim)
      
      # Storing loss information
      Bloss_LR = rep(0, nsim)
      Bloss_LRRI = rep(0, nsim)
      
      # Store success status of convergence (0: failure, 1: success)
      successstore = rep(0, nsim)
      
      
      for (z in 1:nsim) {
        print(sprintf("%dth simulation", z))
        
        # Data generation
        dataset = LVHML:::gendata.func(parameters$Lambda, Kstar, N, J,Tp)
        Y <- dataset$Y
        R <- dataset$R
        Rstore[,, z] <- R
        
        # Estimation process
        #Results under LR method
        coef_LR = t(estglmbeta.func(Y, R, X, Tp, par, n.cores, ext, gamma_fix))
        Uhatstore_LR[,, z] =coef_LR
        
        #Results under LRRI method 
        coef_LRRI = t(estglmm.func(Y, R, X, Tp, par , n.cores ,ext, gamma_fix))
        Uhatstore_LRRI[,,z] <- coef_LRRI
        
        # Evaluation. Compute Bloss and max MSE
        Bloss_LR[z] <- Fnormsum.func(coef_LR[,(g_len+1):(g_len + b_len)] , Bstar)/sqrt(J) 
        Bloss_LRRI[z] <- Fnormsum.func(coef_LRRI[,(g_len+1):(g_len + b_len)] , Bstar)/sqrt(J) 
      }
      
      
      #Compute MSE
      Betahatstore_LR <- Uhatstore_LR[,(g_len+1):(g_len + b_len),]
      Betahatstore_LRRI <- Uhatstore_LRRI[,(g_len+1):(g_len + b_len),]
      Difference_LR =Betahatstore_LR - array(Bstar, dim = c(J,b_len,nsim))
      Difference_LRRI <- Betahatstore_LRRI - array(Bstar, dim = c(J,b_len,nsim))
      MSE_LR = apply(Difference_LR^2,c(1,2), mean  )
      MSE_LRRI = apply(Difference_LRRI^2,c(1,2), mean  )
      max_MSE_LR <- max(MSE_LR)
      max_MSE_LRRI <- max(MSE_LRRI)
      
      # Prepare output list with all relevant information
      out = list(
        "N" = N, 
        "Tp" = Tp, 
        "Rstore" = Rstore,  
        "X" = X, 
        "Ustar" = Ustar, 
        "Uhatstore_LR" = Uhatstore_LR,
        "Uhatstore_LRRI"=Uhatstore_LRRI,
        "niterstore" = niterstore,
        "Improvementsobjstore" = Improvementsobjstore,
        "n.cores" = n.cores, 
        "Bloss_LR" = Bloss_LR,
        "Bloss_LRRI"=Bloss_LRRI,
        "max_MSE_LR"=max_MSE_LR,
        "max_MSE_LRRI"=max_MSE_LRRI,
        "n.cores" = n.cores, 
        "nsim" = nsim,
        "Kstar" = Kstar
      )
      
      # Save the output if it runs successfully until the last simulation
      if (z == nsim) {
        filename = if (par == 0) {
          sprintf("glmZ%dfixg%dext%dJ%dN%dKstar%dnsim%d.Rdata",withZ*1,gamma_fix*1,ext*1,  J, N, Kstar, nsim)
        } else {
          sprintf("glmZ%dfixg%dext%dpancore%dJ%dN%dKstar%dnsim%d.Rdata",withZ*1,gamma_fix*1,ext*1, n.cores, J, N, Kstar, nsim)
        }
        save( out, file = filename)
      }
      
    }
  }
}