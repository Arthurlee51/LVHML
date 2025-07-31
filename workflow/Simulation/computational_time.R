#Script to test computational time of the proposed model against LR and LRRI. For 5 simulations using all models under same settings. The main output is sim_time. 
# --------------------------------------------------------------------------------
# Load required libraries
library(BB)
library(expm)
library(doParallel)
library(Rcpp)
library(truncnorm)
library(glmnet)
library(lme4)
library(LVHML)


# Source functions from external R scripts
source("functions_for_experiment.R")
# --------------------------------------------------------------------------------
# Choosing N relative to J (options: "s" for small and "L" for large)
Tp <- 4#Number of time points. Set as 4 across all settings.
withZ = FALSE#False if we do not have Z, as in simulation study.
ext = FALSE #False if we have A and Beta both time-independent, as in simulation study.
gamma_fix = FALSE#False if we allow gamma_t being time-dependent. TRUE when gamma_t = gamma for all t. Default is FALSE. 
#For each nset, create a data frame storing run time. 
Jlong <- c(100,200,300,400)
sim_time <- matrix(0, nrow = length(Jlong), ncol = 4)
sim_time <- as.data.frame(sim_time)
#Names of methods compared. Proposed_Kstar refers to the proposed method assuming Kstar is known.
colnames(sim_time) <- c("proposed","proposed_Kstar","LR", "LRRI")
for(Kstar in c(3,8)){
  for (nset in c("s","L")){ 
    for (l in 1:length(Jlong)) {
      J <- Jlong[l]
      set.seed(123)
      par = FALSE # FALSE for normal, TRUE for parallel computing
      n.cores = 1 # Number of cores for parallel computing
      nsim = 5 # Number of simulations. 5 for comparing computational time.
      Kset = 1:10 # Range of K values to choose from. Should be 1:10 for simulation purpose.
      
      N <- switch(nset, s = 5 * J, L = 10 * J)
      
      # Generate parameters based on the settings. To be used across all replications.
      parameters <- LVHML:::genpar.func(Kstar,N,J,Tp, withZ, ext, gamma_fix )
      px <- parameters$px
      pz <- parameters$pz
      
      # Initialize storage arrays and variables
      Thetastar = parameters$Thetastar
      Ustar = parameters$Ustar
      #u_len_star: length of u_j. 
      u_len_star <- ncol(Ustar)
      Z = parameters$Z
      X = parameters$X
      Thetahatstore = array(0, dim = c(N, Kstar, nsim))
      Uhatstore = array(0, dim = c(J, u_len_star, nsim))
      Rstore = array(0, dim = c(N, Tp, nsim))
      niterstore = rep(0, nsim)
      Improvementsobjstore = rep(0, nsim)
      objstore = rep(0, nsim)
      ncrit = 10  # Number of IC criteria
      ICstore = array(0, dim = c(nsim, length(Kset), ncrit))
      Khatstore = matrix(0, nrow = nsim, ncol = ncrit)
      Kin = which(Kset == Kstar)  # Index of true K in Kset
      
      # Storing loss information
      Totalloss = rep(0, nsim)
      Totallossfull = matrix(0, nrow = nsim, ncol = Tp)
      Thetaloss = rep(0, nsim)
      Uloss = rep(0, nsim)
      Ulossstore = matrix(0, nrow = nsim, ncol = Tp)
      
      # Store success status of convergence (0: failure, 1: success)
      successstore = rep(0, nsim)
      
      #Store data for fair comparison
      datalist <- list()
      for(z in 1:nsim){
        datalist[[z]] = LVHML:::gendata.func(parameters$Lambda, Kstar, N, J,Tp)
      }
      
   
      # Timing the simulations for proposed
      start_time <- Sys.time()
      for(z in 1:nsim){
        dataset <- datalist[[z]]
        # Data generation
        Y <- dataset$Y
        R <- dataset$R
        # Estimation process
        object = lvhml_est(Y, R,  X , Kset , par, n.cores, Asymp = FALSE, Silent = FALSE, Z ,full = TRUE, ext,gamma_fix)
      }
      end_time <- Sys.time()
      sim_time$proposed[l] <-  as.numeric( end_time - start_time, units = "mins")
      
      # Timing the simulations for proposed_Kstar
      start_time <- Sys.time()
      for(z in 1:nsim){
        dataset <- datalist[[z]]
        # Data generation
        Y <- dataset$Y
        R <- dataset$R
        # Estimation process
        object = lvhml_est(Y, R,  X , Kset = Kstar , par, n.cores, Asymp = FALSE, Silent = FALSE, Z ,full = TRUE, ext,gamma_fix)
      }
      end_time <- Sys.time()
      sim_time$proposed_Kstar[l] <-  as.numeric( end_time - start_time, units = "mins")
      
    
      #Time for LR method
      # Estimation process
      start_time <- Sys.time()
      for(z in 1:nsim){
        dataset <- datalist[[z]]
        # Data generation
        Y <- dataset$Y
        R <- dataset$R
        #Results under LR method
        coef_LR = t(estglmbeta.func(Y, R, X, Tp, par, n.cores,ext, gamma_fix))
      }
      end_time <- Sys.time()
      sim_time$LR[l] <-  as.numeric( end_time - start_time, units = "mins")
      
      #Time for LRRI method
      # Estimation process
      start_time <- Sys.time()
      for(z in 1:nsim){
        dataset <- datalist[[z]]
        # Data generation
        Y <- dataset$Y
        R <- dataset$R
        #Results under LR method
        coef_LRRI = t(estglmm.func(Y, R, X, Tp, par , n.cores,ext, gamma_fix))
      }
      end_time <- Sys.time()
      sim_time$LRRI[l] <-  as.numeric( end_time - start_time, units = "mins")
    }
    
    
    if(z==nsim){
      filename =sprintf("nset%scomputetimeZ%dfixg%dext%dsimKstar%dnsim%d.Rdata",nset, withZ*1, gamma_fix*1, ext*1,  Kstar, nsim )
      save(sim_time , file = filename)
    }
      
  }
}

