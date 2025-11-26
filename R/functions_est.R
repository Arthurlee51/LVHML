#' @useDynLib LVHML, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' @title Estimation Function for lvhml Model for binary response data.
#' @description Main function to implement the proposed estimator in a binary setting, as in simulation and data analysis. This includes choosing the number of factors (K) and calculating asymptotic variances if required. The function is capable of handling time-dependent covariates.
#'
#' @param Y A N x J x Tp array storing binary responses(1 or 0). Entries for unobserved data should be marked as NA.
#' @param R A N x Tp matrix representing the missing indicator (1 for observed, 0 for missing).
#' @param X A N x px matrix storing static covariates.
#' @param Kset A vector containing the candidate numbers of factors. The default value is 1:10.
#' @param par Logical. Indicates whether to use parallel computing. Default is FALSE.
#' @param n.cores The number of cores to use for parallel computing. Used in conjunction with par = TRUE. Default is 1.
#' @param Asymp Logical. Indicates whether to estimate the asymptotic variance. Default is FALSE.
#' @param Silent Logical. Option to suppress output showing progress. Default is FALSE.
#' @param Z A N x J x Tp array for time-dependent covariates. Defaults as empty array if not provided.
#' @param full Logical. Indicates whether full details are reported. Default is FALSE.
#' @param ext Logical. Indicates whether the extension allowing time-dependent loadings and coefficients on X is used. Default is FALSE.
#' @param gamma_fix  Logical. Indicates whether the restriction gamma_{jt} = t gamma_{j} is imposed. Default is FALSE.
#' @param proj_const A scalar specifying the projection constant, also referred to as the constraint parameter in the manuscript. Default is 5.
#' @return A list containing the results of the estimation, including:
#'         \item{Thetahat}{Estimated Theta}
#'         \item{Gammahat}{Estimated Gamma}
#'         \item{Betahat}{Estimated Beta}
#'         \item{Ahat}{Estimated A}
#'         \item{Khat}{Estimated K}
#'         \item{ICs}{Information criteria for values of K in Kset}
#'         \item{Sigma_E}{Estimated asymptotic variance matrices for sqrt(N)*betahat_j (if Asymp = TRUE)}
#'         \item{SE}{The standard errors of Betahat (if Asymp = TRUE)}
#'         \item{AsymVarhatstore}{The estimated variance matrices for sqrt(N)*u_j (if Asymp = TRUE)}
#'         \itme{all_output}{List containing all_output (if full = TRUE)}
#' @examples
#' # Generate example dataset
#' set.seed(123)
#' dataset <- lvhml_ex_data(J = 100, nset = "s", Kstar = 3,ext = FALSE, gamma_fix = FALSE, gap = FALSE)
#' attach(dataset)
#' out <- lvhml_est(Y, R, X, Kset = 1:10, par = FALSE, n.cores = 1, Asymp = FALSE, Silent = FALSE, Z = array(0, dim = c(0, 0, 0)),ext = FALSE,gamma_fix = FALSE)
#'
#' @export
lvhml_est <- function(Y, R, X , Kset = 1:10, par = FALSE, n.cores = 1, Asymp = FALSE, Silent = FALSE, Z = array(0, dim = c(0, 0, 0)), full= FALSE, ext = FALSE,gamma_fix = FALSE, proj_const = 5) {
  # Initialize variables and lists for the estimation process
  all_output <- list()  # List storing output across all values of K
  lk <- length(Kset)  # Number of K values to be tested
  ICs <- rep(0, lk)
  lastobjs <- rep(0, lk)

  # Get dimensions of input matrices
  px <- ncol(X)
  pz <- dim(Z)[2]
  N <- dim(Y)[1]
  J <- dim(Y)[2]
  Tp <- dim(Y)[3]

  #Create Zmat: matrix version of Z.
  Zmat <- GetZmat.func(Z, pz, N, Tp)

  # Convert Y to logical matrix to reduce memory usage
  Ymatlog <- GetYmat.func(Y, N, J, Tp) == 1

  # Indices of time points where r_{it} = 1
  indices_R <- unlist(apply(R, 1, function(row) which(row == 1)))
  lengths_R <- rowSums(R)  # Lengths of the time points

  # Loop through each K value in Kset for estimation
  for (i in seq_along(Kset)) {
    K <- Kset[i]
    if (!Silent)  print(sprintf("Estimating the model at K = %d", K))

    # Compute initial values using SVD-based approach
    init <- SVDinit.func(Y, R, X, Zmat, indices_R, lengths_R, Tp, K, px,pz, N, J, par, n.cores,ext, gamma_fix)

    # Estimation of parameters
    out <- lvhml_est_fit_K(Y, Ymatlog, R, X, Zmat, indices_R, lengths_R, Tp, K,  px,pz, ncol(init$Uinit), N, J, init$Tinit, init$Uinit, par, n.cores, Silent, ext, gamma_fix,proj_const)
    all_output <- append(all_output, list(out))  # Store results for each K

    # Calculate and store IC
    if(ext){
      ICs[i] <- -2 * out$lastobj * N * J + K * max(N, Tp*J) * log(J*sum(R) / max(N, Tp*J))
    }else{
      ICs[i] <- -2 * out$lastobj * N * J + K * max(N, J) * log(J*sum(R) / max(N, J))
    }
    lastobjs[i] <- out$lastobj
  }

  # Select Khat that minimises IC
  Khat <- if (lk == 1) Kset[1] else Kset[which.min(ICs)]

  # Get the data for the chosen Khat and return relevant parameters
  chosen <- all_output[[which.min(ICs)]]
  Thetahat <- chosen$Thetahat
  Uhat <- chosen$Uhat
  sep_U <- Sep_U.func(Uhat,Tp,px,pz,Khat, ext, gamma_fix)

  #Prepare results to be returned
  results <- list(
    Thetahat = Thetahat,
    Gammahat = sep_U$Gamma,
    Betahat  = sep_U$Beta,
    Ahat     = sep_U$A,
    Vhat     = sep_U$V,
    Khat     = Khat,
    ICs      = ICs
  )

  #Add asymptotic variances related quantities if Asymp == TRUE
  if(Asymp){
    AsymVarhatstore <- calasympvar.func(X, Z, Thetahat, Uhat, R, J, N, Tp,ext, gamma_fix)
    #Get variance related to beta. Denote as Sigma_E here.
    g_len  <- if (gamma_fix) 1 else Tp
    Sigma_E <- AsymVarhatstore[(g_len + 1):(g_len + px), (g_len + 1):(g_len + px), ]
    SE <- t(sapply(1:J, function(j) sqrt(diag(Sigma_E[,,j])))) / sqrt(N)
    results <- c(results, list(Sigma_E = Sigma_E, SE = SE, AsymVarhatstore=AsymVarhatstore ))
  }

  # Include full output if requested
  if (full) results$all_output <- all_output

  return(results)
}

#-----------------------------------------------------------------------------------------------------------
#' @title Normalize Model Parameters
#' @description  This function normalizes the model parameters (A, Theta, Gamma, and Beta) to satisfy identifiability conditions.
#'
#' @param J Number of events or outcome variables.
#' @param N Number of subjects.
#' @param px Number of static covariates.
#' @param tA Transposed loading matrix (K times J if ext = FALSE and Tp*K times J if ext = TRUE)
#' @param Theta (N x K) Matrix of latent factors.
#' @param Gamma (J x g_len) Matrix of intercepts. g_len = 1 if gamma_fix = TRUE and g_len = Tp otherwise.
#' @param tBeta Transposed regression coefficient matrix of X (px times J if ext = FALSE and Tp*px times J if ext = TRUE).
#' @param X (N x px) Matrix of static covariates.
#' @param ext Logical. Indicates whether the extension allowing time-dependent loadings and coefficients on X is used. Default is FALSE.
#' @param gamma_fix  Logical. Indicates whether the restriction gamma_{jt} = t gamma_{j} is imposed. Default is FALSE.
#'
#' @return
#'  A list containing normalized versions of tA (tAstar), Theta (Thetastar), Gamma (Gammastar), and tBeta (tBetastar).
#'
#'#' @return A list with normalized components:
#'   \item{tAstar}{Normalized transposed loadings.}
#'   \item{Thetastar}{Normalized latent factors}
#'   \item{Gammastar}{Normalized time-dependent intercept}
#'   \item{tBetastar}{Normalized transposed regression coefficients}
#'@export
lvhml_norm = function(J,N,px,tA,Theta,Gamma,tBeta,X,Tp, ext = FALSE, gamma_fix = FALSE){
  # Initial checks for NA values in critical parameters
  if (any(is.na(tA)) || any(is.na(Theta))) {
    message("NA detected in tA or Theta, normalization cannot proceed.")
    return(list(tAstar = NaN, Thetastar = NaN, Gammastar = NaN, tBetastar = NaN))
  }

    Kstar <- ncol(Theta)
    IN = matrix(1, nrow =N, ncol =1)
    Gamma_ncol<- ncol(Gamma)
    #Compute projection matrix.
    if(!gamma_fix){
      INX <- if (px > 0) cbind(IN, X) else IN
      Proj_mat <- solve(t(INX)%*%INX) %*% t(INX) %*% Theta
      TILTH = Theta - INX %*% Proj_mat
    }else if (px > 0) {
        Proj_mat <- solve(t(X)%*%X) %*% t(X) %*% Theta
        TILTH = Theta - X %*% Proj_mat
      }else{
        TILTH = Theta
      }


    if(!ext){
      # Projection of Theta and calculation of adjusted Gamma and Beta
      #Project only if we have general Gamma or X.
      if(px>0 | !gamma_fix){
        Q = as.matrix(Proj_mat%*%tA)
        if(!gamma_fix){
          Gammastar = Gamma + matrix(Q[1,], nrow=J, ncol=Tp)
          b_start <- 2
        }else{
          Gammastar = Gamma
          b_start <- 1
        }
        if(px>0){
          tBetastar = tBeta
          tBetastar[1:px,] = tBetastar[1:px,] +Q[b_start:nrow(Q),]
        }else{
          tBetastar=matrix(0, nrow=0, ncol=J)
        }
      }

    }else{
      #Initiate Gammastar and tBetastar
      Gammastar <- matrix(0, nrow = J, ncol = ncol(Gamma) )
      tBetastar <- matrix(0, nrow = Tp*px, ncol = J)
      #Project for each t to cater the extension.
      #Only perform this projection if px>0 or gamma_fix is FALSE.
      if(px>0 | !gamma_fix){
      for( t in 1:Tp){
        Q_t = as.matrix(Proj_mat%*%tA[((t-1)*Kstar+1):((t-1)*Kstar+Kstar),])
        if(!gamma_fix){
          Gammastar[,t] <- Gamma[,t]+ Q_t[1,]
          b_start <- 2
        }else{
          Gammastar <- Gamma
          b_start <- 1
        }

        if(px >0){
          tBetastar[((t-1)*px+1):((t-1)*px+px),] <- tBeta[((t-1)*px+1):((t-1)*px+px),]+ Q_t[b_start:nrow(Q_t),]
        }
      }
      }
    }


    # Normalizing A and Theta
    if(!ext){
      SigmaJA = tA%*%t(tA)/J
    }else{
        tA1 <- matrix(tA[1:Kstar,],nrow=Kstar)
      SigmaJA = tA1%*%t(tA1)/J
    }

    SigmaNT = t(TILTH)%*%TILTH/N
    object = eigen(expm::sqrtm(SigmaJA)%*%SigmaNT%*%expm::sqrtm(SigmaJA))
    HNJ = t(solve(expm::sqrtm(SigmaJA))%*%object$vectors)
    if(!ext){
      tAstar = HNJ%*%tA #tAstar%*%t(tAstar)/J Gives Identity matrix
    }else{
      #initialise tAstar
      tAstar <- matrix(0, nrow= Kstar*Tp, ncol = J)
      for(t in 1:Tp){
        tAstar[((t-1)*Kstar+1):((t-1)*Kstar+Kstar),]<- HNJ%*%tA[((t-1)*Kstar+1):((t-1)*Kstar+Kstar),]
      }
    }
    Thetastar= TILTH%*%solve(HNJ)#t(Thetastar)%*%Thetastar/N gives diagonal matrix
  return(list( "tAstar" = tAstar, "Thetastar"= Thetastar, "Gammastar"=Gammastar, "tBetastar"= tBetastar))
}
#-------------------------------------------------------------------------------------------------------------
