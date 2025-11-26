#R script containing internal functions essential for estimation of parameters
#------------------------------------------------------------------------------------------------------------
#' Internal Estimation Function for LVHML Model
#'
#' This function estimates A, Beta, and Theta for a given value of K. It's designed
#' to be used internally within the package and supports parallel computing, handling
#' of time-dependent covariates, and convergence checking.
#'
#' @param Ymatlog Logical matrix (N x (J*Tp)) derived from Y.
#' @param Zmat Numeric matrix (pz x (N*Tp)) derived from time-dependent covariates Z.
#' @param indices_R Vector of indices of observed entries in R.
#' @param lengths_R Number of observed time points for each subject.
#' @param Tp Number of time points.
#' @param K Number of factors to estimate.
#' @param px Number of static covariates.
#' @param pz Number of time-dependent covariates.
#' @param u_len Length. of the parameter u_j.
#' @param N Number of subjects.
#' @param J Number of events.
#' @param Tinit Matrix (N times K) of initial values for Theta.
#' @param Uinit Matrix (J times u_len) of initial values for U =(Gamma,Beta, V, A).
#' @seealso \code{\link{lvhml_est}} for the description of \code{Y}, \code{R}, \code{X}, \code{par}, \code{n.cores}, \code{Silent}, \code{ext}, \code{gamma_fix} and \code{proj_const}.
#' @keywords internal
#'
#' @return A list including:
#'   \item{success}{Convergence indicator (1 = success, 0 = failure).}
#'   \item{Thetahat}{Estimated Theta after final iteration.}
#'   \item{Uhat}{Estimated U after final iteration.}
#'   \item{Uit}{Array of U iterates (J x u_len x niter).}
#'   \item{Thetait}{Array of Theta iterates (N x K x niter).}
#'   \item{Improvementsobj}{Vector of objective improvements per iteration.}
#'   \item{niter}{Number of iterations executed.}
#'   \item{objstore}{Objective values per iteration.}
#'   \item{lastobj}{Objective at final iteration.}
#'   \item{total_proj}{Total number of projection operations performed.}
lvhml_est_fit_K <- function(Y, Ymatlog, R,  X, Zmat, indices_R, lengths_R, Tp, K, px,pz, u_len,N, J, Tinit, Uinit, par, n.cores, Silent, ext, gamma_fix, proj_const) {
  # Initialization
  success <- 1  # To store whether convergence is successful
  maxiter <- 300  # Maximum number of iterations
  Uit <- array(0, dim = c(J, u_len, maxiter))  # Estimates of Gamma, A, and Beta
  Thetait <- array(0, dim = c(N, K, maxiter))  # Theta estimates
  niter <- 1  # Iteration counter
  tol <- 1e-6  # Tolerance for convergence
  reit <- 0  # Re-iteration counter

  # Initial estimates
  Thetait[, , 1] <- Tinit
  Uit[, , 1] <- Uinit

  # Objective function initialization
  obj.func <- function(tU, Theta) rcpp_objprep(Ymatlog, R, tU, Theta, X,Zmat, N, J, px,pz,u_len, K, Tp, ext, gamma_fix)
  Improvementsobj <- rep(0, maxiter - 1)
  objstore <- rep(0, maxiter)
  objstore[1] <- obj.func(t(Uit[, , 1]), as.matrix(Thetait[, , 1]))
  total_proj <- 0
  #Count total number of projections during estimation.
  # Iteration loop
  while (niter == 1 || Improvementsobj[niter - 1] > tol) {
    # if (!Silent) {
    #   print(niter)
    # }
    niter <- niter + 1

    # Max iteration check
    if (niter > maxiter) {
      if (reit == 0) {
        reit <- 1
        print("Max iteration reached, retrying")
        niter <- 2  # Reset to continue the loop properly
        Thetait[, , 1] <- matrix(truncnorm::rtruncnorm(N * K, a = -1, b = 1, mean = 0, sd = 1), N, K)
        Uit[, , 1] <- matrix(truncnorm::rtruncnorm(J * u_len, a = -1, b = 1, mean = 0, sd = 1), J, u_len)
        objstore[1] <- obj.func(t(Uit[, , 1]), as.matrix(Thetait[, , 1]))
      } else {
        print("Max iteration reached, does not converge")
        success <- 0
        break
      }
    }

    # parameter update
    U_niter <- Uest.func(Y, R, Uit[, , niter - 1], as.matrix(Thetait[, , niter - 1]), X,Zmat, indices_R, lengths_R, Tp, K, px,pz, N, J, par, n.cores, ext, gamma_fix,proj_const)
    Uit[, , niter] <- t(U_niter$Upar)
    total_proj <- total_proj + U_niter$proj_count
    Theta_niter <- Thetaest.func(Y, R, as.matrix(Thetait[, , niter - 1]), Uit[, , niter], X,Z, Tp, K, px,pz, N,J,  par, n.cores, ext, gamma_fix,proj_const)
    Thetait[, , niter] <- t(Theta_niter$Thetapar )
    total_proj <- total_proj + Theta_niter$proj_count

    # Finite check for Theta
    if (any(!is.finite(Thetait[, , niter]))) {
      print(sprintf("K = %d, Theta not finite, stopping", K))
      success <- 0
      break
    }

    # Objective value update
    objstore[niter] <- obj.func(t(Uit[, , niter]), as.matrix(Thetait[, , niter]))
    Improvementsobj[niter - 1] <- objstore[niter] - objstore[niter - 1]

    # Handle NA in objective improvement
    if (is.na(Improvementsobj[niter - 1])) {
      print(sprintf("K = %d, improvement NA, does not converge", K))
      success <- 0
      break
    }
  }

  # Finalize estimates
  final_iter <- ifelse(niter > maxiter, maxiter, niter)
  sep_out <- Sep_U.func(Uit[, , final_iter], Tp, px, pz,K, ext, gamma_fix)
  hats <- lvhml_norm(J, N, px, t(sep_out$A), as.matrix(Thetait[, , final_iter]), sep_out$Gamma, t(sep_out$Beta),X,Tp,ext, gamma_fix)
  Thetahat <- hats$Thetastar
  Uhat <- GetU.func(hats$Gammastar,t(hats$tBetastar), sep_out$V, t(hats$tAstar))

  # Return results
  return(list(
    success = success,
    Thetahat = Thetahat,
    Uhat = Uhat,
    Uit = Uit,
    Thetait = Thetait,
    Improvementsobj = Improvementsobj,
    niter = final_iter,
    objstore = objstore,
    lastobj = objstore[final_iter],
    total_proj = total_proj
  ))
}


#' SVD-based low-rank reconstruction of the transformed responses
#'
#' @description Internal helper used in the SVD initialisation step
#'   to obtain a low–rank approximation of the transformed response
#'   matrix.
#'
#' @param mat Numeric matrix (N times M_nrow) containing the transformed responses (stacked over time
#'   if \code{ext = TRUE}).
#' @param N Number of subjects.
#' @param K Latent dimension.
#' @param M_nrow Number of columns of \code{mat} (Typically \code{J} if ext is FALSE and \code{J * Tp} otherwise).
#' @param threshold Thresold to be compared with singular values.
#'
#' @return Numeric matrix of the same dimension as \code{mat} containing the low–rank reconstruction.
#'
#' @keywords internal
run_svd <- function(mat, N, K,M_nrow, threshold) {
  svdresult <- svd(mat, nu = M_nrow, nv = M_nrow)
  tildeK <- max(K + 1, which(svdresult$d > threshold))
  indices <- 1:tildeK
  out <- as.matrix(svdresult$u[, indices]) %*% diag(svdresult$d[indices], nrow = length(indices)) %*% t(svdresult$v[, indices])
  return(out)
}


#' Compute intermediate M-matrix from truncated transformed values
#'
#' @description Internal helper that applies the inverse link
#'   function to truncated transformed responses to build the
#'   intermediate M–matrix.
#'
#' @param tildelijts Numeric matrix with the same dimension as the
#'   transformed response matrix from \code{\link{run_svd}}
#' @param epsilon Truncation constant.
#'
#' @return Numeric matrix \code{M} with the same dimensions as \code{tildelijts}.
#'
#' @keywords internal
compute_M <- function(tildelijts, epsilon) {
  # Compute intermediate M matrix based on tildelijts
  M <- matrix(0, nrow = nrow(tildelijts), ncol = ncol(tildelijts))
  indices <- which((tildelijts >= -1 + epsilon) & (tildelijts <= 1 - epsilon))
  M[indices] <- invxi.func(0.5 * (tildelijts[indices] + 1))
  M[(tildelijts < -1 + epsilon)] <- invxi.func(epsilon)
  M[(tildelijts > 1 - epsilon)] <- invxi.func(1 - epsilon)
  return(M)
}
#----------------------------------------------------------------------------------------------------------------
#' Construct SVD-based initial values for U and Theta
#'
#' @description Internal initialisation routine called by
#'   \code{\link{lvhml_est_fit_K}} prior to the main optimisation.
#'
#' @seealso \code{lvhml_est_fit_K} for the description of parameters.
#'
#' @return A list with two components:
#'   \itemize{
#'     \item \code{Uinit}: \code{J x u_len} matrix of initial values
#'           for the event–specific parameter blocks.
#'     \item \code{Tinit}: \code{N x K} matrix of initial values for
#'           the subject–specific scores.
#'   }
#'
#' @keywords internal
#' @return A list containing initial values for U (Gamma,Beta, V, A) and T (Theta).
SVDinit.func <- function(Y, R, X, Zmat, indices_R, lengths_R, Tp, K, px,pz, N, J, par, n.cores, ext,gamma_fix) {
  epsilon <- 0.01  # Constant used for bounding values

  #Initialization
  M_nrow <- if (!ext) J else J * Tp
  tildeM <- matrix(0, nrow = N, ncol = M_nrow)
  Gammainit <- matrix(0, nrow = J, ncol = Tp)

  # Compute proportion of observed responses across time periods
  phat <- colSums(R) / N

  # Transform Y into -1, 0, 1 based on value and missingness
  lijts <- 2 * Y - 1
  lijts[is.na(lijts)] <- 0

  if(!ext){
    # Apply SVD to transformed Y for each time point
    if(!gamma_fix){
    for (t in 1:Tp) {
      tildelijts <- run_svd(lijts[, , t], N, K,M_nrow, 2 * sqrt(N * phat[t]))
      # Compute intermediate M matrix based on tildelijts
      M <- compute_M(tildelijts,epsilon)
        Gammainit[, t] <- colSums(M) / N
        tildeM <- tildeM + M - matrix(Gammainit[, t], byrow = TRUE, nrow = N, ncol = J)
      }
    }else{
      M_all<- array(0, dim = c(N,J,Tp) )
      #M divided by t to get initial estimate for gamma.
      M_mat<- matrix(0, nrow= N*Tp, ncol = J)
      time_vec <- rep(0, N*Tp)
      for( t in 1:Tp){
        tildelijts <- run_svd(lijts[, , t], N, K,M_nrow, 2 * sqrt(N * phat[t]))
        # Compute intermediate M matrix based on tildelijts
        M <- compute_M(tildelijts,epsilon)
        M_all[,,t] <- M
        M_mat[((t-1)*N+1):(t*N),] <- M
        time_vec[((t-1)*N+1):(t*N)] <- t
      }
      #Get Gammainit derived from removing trend, using the ols estimator.
      Gammainit <-colSums( matrix(time_vec,nrow = N*Tp, ncol =J)*M_mat)/sum(time_vec^2)
      Gammainit <- matrix(Gammainit, nrow = J, ncol = 1)
      for(t in 1:Tp){
        tildeM <- tildeM + M_all[,,t] - matrix(Gammainit*t, byrow = TRUE, nrow = N, ncol = J)
      }
    }
    tildeM <- tildeM / Tp  # Normalize tildeM by the number of timepoints
  }else{
    #Get lijts_mat, where it is (lijts[,,1], lijts[,,2],...lijts[,,Tp])
    lijts_mat <- matrix(lijts, nrow = N, ncol = M_nrow)
    tildelijts <- run_svd(lijts_mat, N, K,M_nrow, 2 * sqrt(N * mean( phat)))
    M <- compute_M(tildelijts,epsilon)
    if(!gamma_fix){
      for(t in 1:Tp){
        #indices of the part related to that time-point
        t_part <- (J*(t-1)+1):(t*J)
        Gammainit[, t] <- colSums(M[,t_part]) / N
        tildeM[,t_part] <- tildeM[,t_part]+ M[,t_part] - matrix(Gammainit[, t], byrow = TRUE, nrow = N, ncol = J)
      }
    }else{
      #Initialise Gammainit, and get initial estimate later.
      M_mat <- matrix(0, nrow = N*Tp, ncol =J)
      time_vec <- rep(0, N*Tp)
      for(t in 1:Tp){
        M_mat[((t-1)*N+1):(t*N),] <- M[,((t-1)*J+1):(t*J)]
        time_vec[((t-1)*N+1):(t*N)] <- t
      }
      #Comoute Gammainit as the least square estimate
      Gammainit <- colSums( matrix(time_vec,nrow = N*Tp, ncol =J)*M_mat)/sum(time_vec^2)
      Gammainit <- matrix(Gammainit, nrow = J, ncol = 1)
      for(t in 1:Tp){
        #indices of the part related to that time-point
        t_part <- (J*(t-1)+1):(t*J)
        tildeM[,t_part] <- tildeM[,t_part]+ M[,t_part] - matrix(t*Gammainit, byrow = TRUE, nrow = N, ncol = J)
      }
    }
  }


  # SVD of tildeM to obtain initial Theta and A
  svdresult <- svd(tildeM, nu = M_nrow, nv = M_nrow)
  Tinit <- matrix(sqrt(N) * svdresult$u[, 1:K], nrow = N, ncol = K)

  #Get a version of Ainit first. Need restrcturing to cater for ext.
  Ainit_temp <- matrix(0, nrow = M_nrow, ncol = K)
  for (k in 1:K) {
    Ainit_temp[, k] <- svdresult$d[k] * svdresult$v[, k]
  }
  Ainit_temp <- Ainit_temp / sqrt(N)  # Normalize Ainit

  if(ext){
    Ainit <- matrix(nrow = J, ncol = K*Tp)
    for(t in 1:Tp){
      Ainit[,((t-1)*K+1):(t*K) ] <- Ainit_temp[((t-1)*J+1):(t*J) , ]
    }
  }else{
    Ainit <- Ainit_temp
  }

  # Obtain initial Beta values through GLM, if applicable
  if ( (px+pz) > 0) {
    BetaandVinit <- t(Betaestglm.func(Y, R, Gammainit, Ainit, Tinit, X, Zmat, indices_R, lengths_R, Tp, K, px,pz, N, J, par, n.cores, ext, gamma_fix))
    g_len <- 0
    if(px>0){
      #Get lengths corresponding to Betainit.
      b_len <- if (!ext) px else px * Tp
      Betainit <- BetaandVinit[,(1+g_len):(b_len+g_len)]
      if (b_len == 1) {
        Betainit <- t(Betainit)
      }
    }else{
        b_len <- 0
      }

    if(pz>0){
      Vinit <- BetaandVinit[,(g_len+b_len+1):(g_len+b_len+pz)]
      if (pz == 1) {
        Vinit <- t(Vinit)
      }
    }
  }
  #Now handle the case of empty X or empty Z.
  if(px==0){
    Betainit <- matrix(0, nrow = J, ncol = 0)
  }
  if(pz==0){
    Vinit<- matrix(0, nrow = J, ncol = 0)
  }

  Uinit <- GetU.func(Gammainit, Betainit,Vinit, Ainit)

  return(list(Uinit = Uinit, Tinit = Tinit))
}
#----------------------------------------------------------------------------------------------------------------
#' Estimate Beta (and V) via GLM for all events
#'
#' @description
#' Internal helper that calls \code{Betaparglm.func} for each event to obtain
#' GLM estimates of the regression parameters, given current values of
#' \code{Gamma}, \code{A} and \code{Theta}.
#'

#' @param Gamma (J x g_len) intercept matrix, where g_len depends on gamma_fix
#' @param A Matrix of Loading parameters (J times K if ext = FALSE and J times Tp*K if ext = TRUE)
#' @param Theta Matrix (N times K) of latent factors.
#' #' @param Y,R,X,Zmat,indices_R,lengths_R,Tp,K,px,pz,N,J,par,n.cores,ext,gamma_fix
#'   See \code{\link{lvhml_est_fit_K}} for the meaning and dimensions of these
#'   arguments.
#'
#' @return
#' A numeric matrix of dimension \code{b_len x J} containing the estimated
#' regression coefficients for each event, where \code{b_len} depends on
#' \code{px}, \code{pz} and \code{ext}.
#'
#' @seealso \code{\link{Betaparglm.func}} for the per-event GLM update.
#' @keywords internal
Betaestglm.func <- function(Y, R, Gamma, A, Theta, X, Zmat, indices_R, lengths_R, Tp, K, px,pz, N, J, par, n.cores, ext, gamma_fix) {
  if (par) {
    out <- parallel::mcmapply(
      function(j) as.vector(Betaparglm.func(t(Y[, j, ]), indices_R, lengths_R, Gamma[j, ], A[j, ], Theta, X, Zmat, N, px,pz, K, Tp, ext, gamma_fix)),
      1:J, mc.cores = n.cores
    )
  } else {
    out <- sapply(
      1:J, function(j) as.vector(Betaparglm.func(t(Y[, j, ]), indices_R, lengths_R, Gamma[j, ], A[j, ], Theta, X, Zmat, N, px,pz, K, Tp, ext, gamma_fix))
    )
  }
  return(matrix(out, nrow = length(out) / J, ncol = J))
}

#' Beta (and V) Parameter Estimation with GLM
#'
#' This internal function performs the estimation of Beta parameters for a single event using GLM.
#' It's utilized by \code{Betaestglm.func} for either sequential or parallel processing.
#'
#' @param yris Transposed Y (Tp times N) matrix for the jth event.
#' @param gammaj Vector of gamma coefficients for the jth event.
#' @param aj Vector of loading coefficients for the jth event.
#' @seealso \code{lvhml_est_fit_K} for the description of other parameters.
#'
#' @keywords internal
#'
#' @return A vector of coefficients estimated by the GLM.
Betaparglm.func <- function(yris, indices_R, lengths_R, gammaj, aj, Theta, X, Zmat, N, px,pz, K, Tp, ext, gamma_fix) {
  yj <- yris[!is.na(yris)]
  id <- rep(1:N, times = lengths_R)

  if(!ext){
    Thetaaj <- Theta %*% aj
    Thetaaj_off <- Thetaaj[id]
  }else{
    Thetaaj_mat <- matrix(0, nrow = N, ncol =Tp)
    for(t in 1:Tp){
      Thetaaj_mat[,t] <- Theta%*%aj[ ((t-1)*K + 1):(t*K) ]
    }
    #part of offset to to used.
    Thetaaj_off <- Thetaaj_mat[t(!is.na(yris))]
  }

  Eforj <- rcpp_BetagetEforj(indices_R - 1, lengths_R, X, Zmat, N, Tp, px,pz, ext,gamma_fix)$Eforj
  #If gamma_fix =TRUE, the first column of Eforj correspond to tIN.
  if(!gamma_fix){
    glmfit <- glm.fit(Eforj, yj, family = binomial(), offset = Thetaaj_off + gammaj[indices_R])
  }else{
    glmfit <- glm.fit(Eforj[,-1], yj, family = binomial(), offset = Thetaaj_off + gammaj*indices_R )
  }


  coef <- glmfit$coefficients

  return(coef)
}

#--------------------------------------------------------------------------------------------------------------
#' Update the Values of U (Gamma,Beta, V, A)
#'
#' This function updates the values of U based on the current estimates.
#' It supports parallel processing for efficiency.
#'
#' @param U Current (J x u_len) matrix of event-specific parameters,
#'   where each row contains the stacked vector \code{(Gamma_j, Beta_j, V_j, A_j)}.
#' @param Theta Current (N x K) matrix of latent factors.
#' @seealso \code{lvhml_est_fit_K} for the description of other parameters.
#'
#' @keywords internal
#'
#' @return A matrix containing the updated value of U.
Uest.func <- function(Y, R, U, Theta, X, Zmat, indices_R, lengths_R, Tp, K, px,pz, N, J, par, n.cores, ext, gamma_fix,proj_const) {
  if (par) {
    out <- parallel::mcmapply(
      function(j) as.vector(Upar.func(t(Y[, j, ]), indices_R, lengths_R, U[j, ], Theta, X, Zmat, N, px,pz, K, Tp,ext, gamma_fix,proj_const)),
      1:J, mc.cores = n.cores, SIMPLIFY = FALSE
    )
  } else {
    out <- lapply(
      1:J, function(j) as.vector(Upar.func(t(Y[, j, ]), indices_R, lengths_R, U[j, ], Theta, X, Zmat, N, px,pz, K, Tp,ext, gamma_fix,proj_const))
    )
  }
  #Get back parameter update and projection count updates, respectively.
  Upar <- sapply(1:J,function(j) rbind(out[[j]]$finalPar))#Upar is deliberate set as t(Upar) for simplicit
  proj_count <- sum(sapply(1:J, function(j) out[[j]]$proj_count))
  out <- list("Upar" = Upar, "proj_count"=proj_count)
  return(out)
}

#' Parallel Function for U Parameter Update
#'
#' Performs the update of U parameters for a single event. It calculates gradients and
#' Hessians to find the optimal update direction and magnitude for U parameters.
#'
#' @param yris (Tp x N) Transposed Y matrix for the jth event.
#' @param uj Current estimates of the jth row of U for the event.
#' @param Theta (N times K) Current estimates of Theta.
#'
#'@seealso \code{lvhml_est_fit_K} for the description of other parameters.
#'
#' @keywords internal
#'
#' @return A vector of updated U parameters for the event.
Upar.func <- function(yris, indices_R, lengths_R, uj, Theta, X, Zmat, N, px,pz, K, Tp, ext, gamma_fix,proj_const) {
  yj <- yris[!is.na(yris)]
  uj_len <- length(uj)

  # Adjust indices for Rcpp and calculate necessary matrices for optimization
  touse <- rcpp_getEforjEujXiM0(uj, indices_R - 1, lengths_R, Theta, X, Zmat, N, Tp, K, px,pz, uj_len, ext,gamma_fix)
  Eforj <- touse$Eforj
  Euj <- touse$Euj
  XiM0 <- touse$XiM0
  nonzero_ind <- touse$nonzero_ind#indices of non-zero entries in terms of rcpp (starts from 0)

  # Calculate gradient and Hessian
  grad <- rcpp_cal_grad(yj, XiM0, Eforj, N, uj_len)
  Hes <- rcpp_cal_Hes_foru(XiM0, Eforj, uj, N, length(XiM0),nonzero_ind )

  # Determine descent direction
  drt <- get_drt.func(Hes, grad, uj_len)

  # Update jth row of U.
  out <- rcpp_getnewpar_func_foru(XiM0, yj, N, grad, drt, uj, Eforj, rep(0, length(XiM0)),nonzero_ind ,proj_const)
  return(out)
}

#-----------------------------------------------------------------------------------------------------------------
#' Update Theta Parameters
#'
#' This function updates the Theta parameters based on the current estimates of U, Theta, and other model parameters.
#' It supports parallel processing for efficient computation.
#'
#' @param U Current (J x u_len) matrix of event-specific parameters,
#'   where each row contains the stacked vector \code{(Gamma_j, Beta_j, V_j, A_j)}.
#' @param Theta Current (N x K) matrix of latent factors.
#'
#'@seealso \code{lvhml_est_fit_K} for the description of other parameters.
#'
#' @keywords internal
#'
#' @return A matrix of updated Theta parameters for each subject.
Thetaest.func = function(Y,R,Theta,U,X, Z,Tp,K,px,pz,N,J,par,n.cores, ext, gamma_fix,proj_const){
  out <- Sep_U.func(U, Tp, px, pz,K, ext, gamma_fix)
  Gamma <- out$Gamma
  Beta <- out$Beta
  V <- out$V
  A <- out$A
  if(par){
    out =parallel::mcmapply(function(i) as.vector(Thetapar.func(t(Y[i,,]), R[i,],Theta[i,],A,Gamma,Beta,V,X,Z,K,Tp,px,pz,i, J, ext, gamma_fix,proj_const)) ,1:N,mc.cores = n.cores, SIMPLIFY = FALSE )
  }else{
    out= lapply(1:N, function(i) as.vector(Thetapar.func(t(Y[i,,]), R[i,],Theta[i,],A,Gamma,Beta,V,X,Z,K,Tp,px,pz,i, J, ext, gamma_fix,proj_const)) )
  }
  Thetapar <- sapply(1:N,function(i) rbind(out[[i]]$finalPar))
  proj_count <- sum(sapply(1:N, function(i) out[[i]]$proj_count))
  out <- list("Thetapar" = Thetapar, "proj_count"=proj_count)
  return(out)
}


#' Parallel Function for Theta Parameter Update
#'
#' Computes the update for Theta parameters for a single subject using current model estimates.
#'
#' @param yrjs Transposed Y matrix (Tp x J) for the ith subject.
#' @param Ri Missing indicator vector for the ith subject.
#' @param thetai Current Theta estimate for the ith subject.
#' @param Gamma (J x Tp) intercept matrix
#' @param A Matrix of Loading parameters (J times K if ext = FALSE and J times Tp*K if ext = TRUE)
#' @param Beta Current estimate of Beta (J times px if ext = FALSE and J times Tp*px if ext = TRUE).
#' @param i Subject index.
#'
#'#'@seealso \code{lvhml_est_fit_K} for the description of other parameters.
#'
#' @keywords internal
#'
#' @return A vector of updated Theta parameters for the subject.
Thetapar.func = function(yrjs, Ri,thetai,A,Gamma,Beta,V,X,Z,K,Tp,px,pz,i, J, ext , gamma_fix,proj_const){
  # Extract observed responses
  yi = yrjs[!is.na(yrjs)]

  if(!ext){
    id = rep(1:J, each=sum(Ri))
    # Prepare matrices and offsets for optimization
    tAfori= t(A[id,])
  }else{
    sum_Ri <- sum(Ri)
    #Initiate tAfori
    tAfori <- matrix(0, nrow = K, ncol =sum_Ri*J)
    #get indices of positive ri
    pos_ri <- which(Ri==1)
    for(z in 1:sum_Ri){
      #get t: the corresponding time-point
      t <- pos_ri[z]
      tAfori[,sum_Ri*(0:(J-1))+z ] <- t(A[,((t-1)*K+1):(t*K) ])
    }
  }

  ind = yrjs
  ind[!is.na(ind)]=1
  if(!gamma_fix){
    tGammaprep = t(Gamma)*ind
  }else{
    tGammaprep = t(matrix(Gamma, nrow = J, ncol=Tp) )*matrix(1:Tp, nrow = Tp, ncol = J)*ind
  }

     if(px>0){
      if(!ext){
        tBetaXiprep<-  matrix(Beta%*%X[i,], nrow = Tp, ncol = J, byrow = TRUE)*ind
      }else{
        tBetaXiprep <- ind
        for(z in 1:sum_Ri){
          t <- pos_ri[z]
          tBetaXiprep[t,] <- Beta[,((t-1)*px+1):(t*px)]%*%X[i,]
        }
      }
    }else{
      tBetaXiprep <- ind*0
    }

    if(pz>0){
      tVZiprep =  t(V%*%Z[i,,])*ind
    }else{
      tVZiprep =  ind*0
    }

    tooff = tGammaprep[!is.na(ind)] +tBetaXiprep[!is.na(ind)] + tVZiprep[!is.na(ind)]

  #Compute gradient and Hessian
  eUi  = as.vector(thetai%*%tAfori)+tooff
  XiM0=xi.func(eUi)
  grad<- rcpp_cal_grad(yi,XiM0, t(tAfori), J,K)
  Hes<-rcpp_cal_Hes(as.vector(XiM0), t(tAfori), thetai, J, length(XiM0))
  drt <- get_drt.func(Hes,grad,K)

  #Update parameter by line search
  out = rcpp_getnewpar_func(XiM0,yi,J,grad,drt, thetai,t(tAfori ),tooff,proj_const)

  return(out)
}

#-----------------------------------------------------------------------------------------------------------------
#'Function to get Zmat
#'@keywords internal
#'@param Z  (N x pz x Tp) Array of time-dependent covariates.
#'@seealso \code{lvhml_est_fit_K} for the description of other parameters.
GetZmat.func = function(Z,pz,N,Tp){
  if(pz>0){
    # Initialize an empty matrix with dimensions p * (N*Tp), such that XZmat[, ((i-1)*Tp+1):(i*Tp)] refers to the p \times Tp matrix for the ith individual
    Zmat <- matrix(nrow = pz, ncol = N * Tp)
    # Fill the matrix using a loop over N
    for (i in 1:N){
      Zmat[, ((i-1)*Tp+1):(i*Tp)] <- Z[i,,]
    }
  }else{
    Zmat = matrix(0, nrow=0, ncol=0)
  }
  return(Zmat)
}

#Function to create Ymat
#'@seealso \code{lvhml_est_fit_K} for the description of parameters.
#'@keywords internal
GetYmat.func = function(Y,N,J,Tp){
  Ymat <- matrix(nrow = N, ncol = J * Tp)
  # Fill the matrix using a loop over J
  for (j in 1:J){
    Ymat[, ((j-1)*Tp+1):(j*Tp)] <- Y[,j,]
  }
  return(Ymat)
}


#Function to get descend direction from Hessian matrix and gradient
#' Compute Newton search direction with regularised Hessian
#'
#' @description Internal helper that computes the search direction.
#'
#' @param Hes (par_len x par_len) Hessian matrix.
#' @param grad Numeric vector of length par_len containing the
#'   gradient.
#' @param par_len Integer; length of the parameter vector.
#'
#' @return Numeric vector of length \code{par_len} giving the search
#'   direction.
#'
#' @keywords internal
get_drt.func = function(Hes,grad,par_len){
  drt = NULL

  # Try solving the Hessian directly
  invHes = tryCatch({
    solve(Hes)
  }, error = function(e) {
    NULL
  })

  # If direct inversion failed, regularize the Hessian
  if (is.null(invHes)) {
    message("Matrix inversion failed; regularizing Hessian.")
    Hes_reg = Hes - 0.01 * diag(par_len)
    invHes = solve(Hes_reg)
  }

  drt = as.vector(-invHes %*% grad)

  return(drt)
}

#inverse function of the logistic function. p can be scalar, vector or matrix.
#'@keywords internal
invxi.func = function(p){
  out = log(p/(1-p))
  return(out)
}


#This script contains the function to generate data for simulation
#Function to combine Gamma,Beta, A to form U. Gamma, Beta and A should all have J rows.
#'@keywords internal
GetU.func = function(Gamma,Beta,V,A){
    return(cbind(Gamma,Beta,V,A))
}


#Function to separate U and get back Gamma, Beta, V and A. U should have J rows.
#'@keywords internal
Sep_U.func = function(U, Tp, px, pz,K, ext, gamma_fix){
  J <- nrow(U)
  g_len <- if(gamma_fix) 1 else Tp
  Gamma <- as.matrix(U[,1:g_len])
  if(px>0){
    b_last <- if(ext) g_len+Tp*px else g_len+px
    Beta <- as.matrix(U[,(g_len+1):b_last])
  }else{
    Beta <- matrix(0, nrow=J, ncol=0)
  }

  if(pz >0){
    V <- as.matrix(U[,(b_last+1):(b_last+pz)])
  }else{
    V <- matrix(0, nrow=J, ncol=0)
  }

  if(!ext){
    A <- as.matrix(U[,(b_last+pz+1):(b_last+pz+K)])
  }else{
    A <- as.matrix(U[,(b_last+pz+1):(b_last+pz+Tp*K)])
  }

return(list("Gamma" = Gamma, "Beta" = Beta, "V" = V, "A" =A))
}


#logit function. x can be scalar, vector or matrix.
#'@keywords internal
xi.func = function(x){
  expx = exp(x)
  out=expx/(1 + expx)
  #Make sure it does not explode
  out[x>700]=1
  out[x< -800]=0
  return(out)
}

