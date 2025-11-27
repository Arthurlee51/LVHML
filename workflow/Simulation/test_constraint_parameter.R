#Script to analyze the impact of changing the constraint parameter from 3 to 7, as in Table S2.
# Main script to run simulations
# --------------------------------------------------------------------------------

# Load required libraries
library(BB)
library(expm)
library(doParallel)
library(Rcpp)
library(truncnorm)
library(glmnet)
library(LVHML)

source("functions_for_experiment.R")
# --------------------------------------------------------------------------------
# Choosing N relative to J (options: "s" for small and "L" for large)
Tp <- 4#Number of time points. Set as 4 across all settings.
withZ = FALSE#False if we do not have Z, as in simulation study.
ext = FALSE #False if we have A and Beta both time-independent, as in simulation study.
gamma_fix = FALSE#False if we allow gamma_t being time-dependent. TRUE when gamma_t = gamma for all t. Default is FALSE. 
proj_const_range = 3:7;#Range of constraint parameter to be tested.
# Run simulations 
for(nset in c("s","L")){
for (J in c(100,200,300,400)) {
  set.seed(123)
  par = FALSE # FALSE for normal, TRUE for parallel computing
  n.cores = 1 # Number of cores for parallel computing
  nsim = 100 # Number of simulations. Usually 100. 
  for(Kstar in c(3,8)){
  #Generate N
  if(nset =="s"){
    N = 5*J
  }
  if(nset =="L"){
    N = 10*J
  }
  
  # Generate parameters based on the settings
  parameters <- LVHML:::genpar.func(Kstar,N,J,Tp, withZ, ext, gamma_fix)
  px <- parameters$px
  pz <- parameters$pz
  

  # Initialize storage arrays and variables
  Thetastar = parameters$Thetastar
  Ustar = parameters$Ustar
  #u_len_star: length of u_j. 
  u_len_star <- ncol(Ustar)
  Z = parameters$Z
  X = parameters$X
  out <- list()
  mean_totalloss <- rep(0, length(proj_const_range))
  avg_projstore <- rep(0, length(proj_const_range))
  
  
  for(proj_count in 1:length(proj_const_range) ){#Test the impact using different projection constants.
    proj_const <- proj_const_range[proj_count]
  Thetahatstore = array(0, dim = c(N, Kstar, nsim))
  Uhatstore = array(0, dim = c(J, u_len_star, nsim))
  Rstore = array(0, dim = c(N, Tp, nsim))
  niterstore = rep(0, nsim)
  Improvementsobjstore = rep(0, nsim)
  objstore = rep(0, nsim)
  ncrit = 10  # Number of IC criteria
  
  # Storing loss information
  Totalloss = rep(0, nsim)
  Totallossfull = matrix(0, nrow = nsim, ncol = Tp)
  Thetaloss = rep(0, nsim)
  Uloss = rep(0, nsim)
  Ulossstore = matrix(0, nrow = nsim, ncol = Tp)
  
  #Matrix storing the total number of projections performed.
  total_proj_store = rep(0, nsim)
  # Store success status of convergence (0: failure, 1: success)
  successstore = rep(0, nsim)
  
  # Timing the simulations
  start_time <- Sys.time()
  
  for (z in 1:nsim) {
    print(sprintf("%dth simulation", z))
    
    # Data generation
    dataset = LVHML:::gendata.func(parameters$Lambda, Kstar, N, J,Tp)
    Y <- dataset$Y
    R <- dataset$R
    Rstore[,, z] <- R
    
    # Estimation process
    object = lvhml_est(Y, R,  X , Kset  = Kstar, par, n.cores, Asymp = FALSE, Silent = FALSE, Z ,full = TRUE, ext,gamma_fix, proj_const)
    
    #Extract info from the estimate.
    toeval = object$all_output[[1]]
    total_proj_store[z] <- toeval$total_proj
    # Store results for evaluation
    successstore[z] = toeval$success
    Thetahatstore[,, z] = toeval$Thetahat
    Uhatstore[,, z] = toeval$Uhat
    
    niterstore[z] = toeval$niter
    Improvementsobjstore[z] = toeval$Improvementsobj[niterstore[z]]
    objstore[z] = toeval$objstore[niterstore[z]]
    
    # Evaluation. 
    evalout = eval.func(toeval$Thetahat, toeval$Uhat, Thetastar, Ustar, X, Z, N, J, Tp, px,pz, Kstar, ext, gamma_fix)
    
    Totalloss[z] = evalout$Totalloss
    Thetaloss[z] = evalout$Thetaloss
    Totallossfull[z, ] = evalout$Totallossstore
    Uloss[z] = evalout$Uloss
    Ulossstore[z, ] = evalout$Ulossstore
  }
  
  # Capture last simulation data and average number of iterations
  lasty = Y
  avgniter = mean(niterstore)
  avg_proj <- mean(total_proj_store)
  # Calculate total simulation time
  end_time <- Sys.time()
  testtime = end_time - start_time
  
  mean_totalloss[proj_count] <- mean(Totalloss)
  avg_projstore[proj_count] <- avg_proj 
  
  #Compute Bloss and MMSE
  g_len <- if(!gamma_fix) Tp else 1
  Beta_range <- (g_len+1):(g_len+px)
  Betastar<- Ustar[,Beta_range]
  Betahatstore <- Uhatstore[,Beta_range,] 
  Bloss = rep(0, nsim)
  Bloss <- sapply(1:nsim, function(z) Fnormsum.func(Betahatstore[,,z] , Betastar) )/sqrt(J)
  
  Difference =Betahatstore - array(Betastar, dim = c(J,px,nsim))
  MSE = apply(Difference^2,c(1,2), mean  )
  MMSE = max(MSE)
  # Prepare output list with all relevant information
  out[[paste("proj_const",proj_const, sep="")]] = list(
    "N" = N, 
    "Tp" = Tp, 
    "Rstore" = Rstore,  
    "lasty" = lasty,
    "X" = X, 
    "Z" = Z,
    "successstore" = successstore,
    "objstore" = objstore, 
    "avgniter" = avgniter, 
    "Thetahatstore" = Thetahatstore, 
    "Thetastar" = Thetastar,
    "Ustar" = Ustar, 
    "Uhatstore" = Uhatstore,  
    "niterstore" = niterstore,
    "Improvementsobjstore" = Improvementsobjstore,
    "n.cores" = n.cores, 
    "Totalloss" = Totalloss, 
    "Totallossfull" = Totallossfull,
    "Thetaloss" = Thetaloss,
    "Uloss" = Uloss, 
    "Bloss"= Bloss,
    "MMSE" = MMSE,
    "n.cores" = n.cores, 
    "nsim" = nsim,
    "Kstar" = Kstar,
    "total_proj_store"= total_proj_store,
    "avg_proj"=avg_proj
  )
  }
  # Save the output if it runs successfully until the last simulation
  if (z == nsim) {
    filename = if (par == 0) {
      sprintf("testprojZ%dfixg%dext%dsimJ%dN%dKstar%dnsim%d.Rdata", withZ*1, gamma_fix*1, ext*1, J, N, Kstar, nsim )
    } else {
      sprintf("testprojZ%dfixg%dext%dpancore%dsimJ%dN%dKstar%dnsim%d.Rdata",withZ*1,gamma_fix*1, ext*1, n.cores,  J, N, Kstar, nsim )
    }
    save(testtime, out,   mean_totalloss,avg_projstore,file = filename)
  }
  
  # Print the total simulation time
  print(testtime)
}
}
}
