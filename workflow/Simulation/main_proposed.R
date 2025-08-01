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

# # Source functions from external R scripts
 source("functions_for_experiment.R")
# --------------------------------------------------------------------------------
# Choosing N relative to J (options: "s" for small and "L" for large)
Tp <- 4#Number of time points. Set as 4 across all settings.
withZ = FALSE#False if we do not have Z, as in simulation study.
ext = FALSE #False if we have A and Beta both time-independent, as in simulation study.
gamma_fix = FALSE#False if we allow gamma_t being time-dependent. TRUE when gamma_t = gamma for all t. Default is FALSE. 
gap = FALSE #FALSE if the setting in the main script is used. TRUE if the setting (Section E.1) with sufficient eigengap for identification of latent factors is applied.
# Run simulations for different values of J
for (nset in c("s")){ 
  for (J in c(100)) {
    set.seed(123)
    par = FALSE # FALSE for normal, TRUE for parallel computing
    n.cores = 1 # Number of cores for parallel computing
    nsim = 2 # Number of simulations. Usually 100
    Kset = 1:10 # Range of K values to choose from. 
    Kstar = 3  # True value of K
    N <- switch(nset, s = 5 * J, L = 10 * J)
    
    # Generate parameters based on the settings. To be used across all replications. (Internal function for simulation purpose only)
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
    
    
    for (z in 1:nsim) {
      print(sprintf("%dth simulation", z))
      
      # Data generation
      dataset = LVHML:::gendata.func(parameters$Lambda, Kstar, N, J,Tp)
      Y <- dataset$Y
      R <- dataset$R
      Rstore[,, z] <- R
      
      # Estimation process
      object = lvhml_est(Y, R,  X , Kset , par, n.cores, Asymp = FALSE, Silent = FALSE, Z ,full = TRUE, ext,gamma_fix)
      
      # Store results for evaluation
      ICstore[z, , ] = object$ICs
      Khatstore[z, ] = object$Khat
      all_output = object$all_output
      toeval = all_output[[Kin]]  # Get the object with Kstar for evaluation
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
    

    # Prepare output list with all relevant information
    out = list(
      "N" = N, 
      "Tp" = Tp, 
      "Rstore" = Rstore,  
      "Kset" = Kset, 
      "ICstore" = ICstore,  
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
      "Khatstore" = Khatstore,
      "Totalloss" = Totalloss, 
      "Totallossfull" = Totallossfull,
      "Thetaloss" = Thetaloss,
      "Uloss" = Uloss, 
      "n.cores" = n.cores, 
      "nsim" = nsim,
      "Kstar" = Kstar
    )
    
    # Save the output if it runs successfully until the last simulation
    if (z == nsim) {
      filename = if (par == 0) {
        sprintf("Z%dfixg%dext%dsimJ%dN%dKstar%dnsim%d.Rdata", withZ*1, gamma_fix*1, ext*1, J, N, Kstar, nsim )
        if(gap==TRUE){
          sprintf("gapZ%dfixg%dext%dsimJ%dN%dKstar%dnsim%d.Rdata", withZ*1, gamma_fix*1, ext*1, J, N, Kstar, nsim )
        }
      } else {
        # sprintf("pancore%dsetting%dJ%dN%dKstar%dnsim%d.Rdata", n.cores, setting, J, N, Kstar, nsim)
        sprintf("Z%dfixg%dext%dpancore%dsimJ%dN%dKstar%dnsim%d.Rdata",withZ*1,gamma_fix*1, ext*1, n.cores,  J, N, Kstar, nsim )
        if(gap==TRUE){
          sprintf("gapZ%dfixg%dext%dpancore%dsimJ%dN%dKstar%dnsim%d.Rdata",withZ*1,gamma_fix*1, ext*1, n.cores,  J, N, Kstar, nsim )
        }
      }
      save( out, file = filename)
    }
  }
}
