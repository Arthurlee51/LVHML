#include <Rcpp.h>
using namespace Rcpp;

//Funrtion to calculate hessian matrix.
//Arguments:
//XiM0: Vector of fitted probabilites with length XiM0_len
//Covariate: (XiM0_len x par_len) Matrix of cavariate
//parvec: paremeter vector of length par_len
//Use cases:
// Covariate: Eforj or Afori, parvec: uj or thetai, len:N or J
// [[Rcpp::export]]
NumericMatrix rcpp_cal_Hes(NumericVector XiM0, NumericMatrix Covariate, NumericVector parvec, int len, int XiM0_len) {
  int par_len = parvec.size();
  NumericMatrix Hes(par_len, par_len);
  NumericVector commonpart = -XiM0*(1 -XiM0 ) ;
  double common;
  double cov_zk;
  //double product;
  for(int z =0; z< XiM0_len; ++z){
    common = commonpart[z];
    for (int k = 0; k < par_len; ++k) {
      cov_zk = Covariate(z, k);
      if(cov_zk==0){
        //skip if cov_zk =0
      }else{
        for (int l = k; l < par_len; ++l) {
          if(Covariate(z, l) ==0){
            //skip if Covariate(z, l) =0
          }else{
            Hes(k, l) += common *cov_zk *Covariate(z, l);
          }

        }
      }

    }
  }

  //Divide value by len
  for (int k = 0; k < par_len; ++k) {
    for (int l = k; l < par_len; ++l) {
      Hes(k, l) = Hes(k, l) / len;
    }
  }

  // Mirror the upper triangle to the lower triangle
  for (int i = 0; i < par_len; ++i) {
    for (int j = i+1; j < par_len; ++j) {
      Hes(j, i) = Hes(i, j);
    }
  }
  return Hes;
}


//Function to calculate hessian matrix for uj with extra parameters to take care of ext.
// nonzero_ind:NumericMatrix of size XiM0_len x ut_len. Row z contains the column indices of the non-zero entries of Covariate for observation z.
//The rest of the arguments are similar to rcpp_cal_Hes.
// [[Rcpp::export]]
NumericMatrix rcpp_cal_Hes_foru(NumericVector XiM0, NumericMatrix Covariate, NumericVector parvec, int len, int XiM0_len, NumericMatrix nonzero_ind) {
  int par_len = parvec.size();
  NumericMatrix Hes(par_len, par_len);
  NumericVector commonpart = -XiM0*(1 -XiM0 ) ;
  double common;
  double cov_zk;
  int ut_len = nonzero_ind.ncol();
  int k;
  int kl;
  for(int z =0; z< XiM0_len; ++z){
    common = commonpart[z];
    for (int j = 0; j < ut_len ; ++j) {
      k = nonzero_ind(z,j);
      cov_zk = Covariate(z, k);
        for (int l = j; l < ut_len; ++l) {
          kl = nonzero_ind(z,l);
            Hes(k, kl) += common *cov_zk *Covariate(z, kl);
        }

    }
  }

  //Divide value by len
  for (int k = 0; k < par_len; ++k) {
    for (int l = k; l < par_len; ++l) {
      Hes(k, l) = Hes(k, l) / len;
    }
  }

  // Mirror the upper triangle to the lower triangle
  for (int i = 0; i < par_len; ++i) {
    for (int j = i+1; j < par_len; ++j) {
      Hes(j, i) = Hes(i, j);
    }
  }
  return Hes;
}


//calculate grad for uj or thetai
//yvec: yi or yj of length len depending on context, Covariate (len x parlen): Eforj or Afori depending on context, len: N or J, parlen: length of aj of thetai
// [[Rcpp::export]]
NumericVector rcpp_cal_grad(NumericVector yvec, NumericVector XiM0, NumericMatrix Covariate, int len, int parlen) {
  NumericVector grad(parlen);

  for (int j = 0; j < parlen; ++j) {
   grad[j] = sum((yvec - XiM0) * Covariate(_, j)) / len;
  }

  return grad;
}


//Function to compute value of objective function related to the par uj or thetai:
//Xi: Fitted probabilities of length len. yvec: yi or yj. len: N or J.
// [[Rcpp::export]]
double rcpp_objpar_prep(NumericVector XiM, NumericVector yvec, int len) {
  double out = 0.0;
  int m = XiM.size();

  for (int j = 0; j < m; ++j) {
    if (yvec[j]) {
      out += log(XiM[j]);
    } else{
      out += log(1 - XiM[j]);
    }
  }

  return out / len;
}


//Logistic link function. x can be NumericVector of any length.
// [[Rcpp::export]]
NumericVector rcpp_xi_func(NumericVector x) {
  return 1 / (1 + exp(-x));
}


// Compute the value of objective function based on likelihood.
//tU (u_len times J)is the tranpose of U. The dimensions and interpretation of the rest of the arguments are consistent with the R function lvhml_est_fit_K.
// [[Rcpp::export]]
double rcpp_objprep(LogicalMatrix Ymatlog, IntegerMatrix R, NumericMatrix tU, NumericMatrix Theta, NumericMatrix X, NumericMatrix Zmat, int N, int J, int px,int pz,int u_len, int K, int Tp, bool ext, bool gamma_fix) {
  double objectiveValue = 0;

  // Extracting submatrices for Gamma and UwlGamma from tU
  NumericMatrix Gamma;
  int g_len ;
  if(!gamma_fix){
    Gamma = tU(Range(0, Tp - 1), Range(0, tU.ncol() - 1));
    g_len = Tp;
  }else{
    Gamma = tU(Range(0,0) , Range(0, tU.ncol() - 1));
    g_len = 1;
  }

  NumericMatrix UwlGamma = tU(Range(g_len, u_len - 1), Range(0, tU.ncol() - 1));

  // Iterate over all subjects, time points, and variables
  for (int i = 0; i < N; i++) {
    for (int t = 0; t < Tp; t++) {
      if (R(i, t) == 1) { // Only proceed if R(i, t) is 1
        for (int j = 0; j < J; j++) {
          //Initialise variables
          double thetaComponent = 0;
          double gammaComponent=0;
          double xzComponent = 0;

          if(!ext){
            for (int k = 0; k < K; k++) {
              thetaComponent += Theta(i, k) * UwlGamma(px+pz + k, j);
            }
          }else{
            for (int k = 0; k < K; k++) {
              thetaComponent += Theta(i, k) * UwlGamma(Tp*px+pz+ t*K+ k, j);
            }
          }



          if (px+pz > 0) {
            //First add the x component, then the z component
            if(px>0){
              if(!ext){
                for (int l = 0; l < px; l++) {
                  xzComponent += X(i,l) * UwlGamma(l, j);
                }
              }else{
                for (int l = 0; l < px; l++) {
                  xzComponent += X(i,l) * UwlGamma(t*px+ l, j);
                }
              }

            }
            if(pz >0){
              if(!ext){
                for (int l = 0; l < pz; l++) {
                  xzComponent += Zmat(l, (i * Tp + t)) * UwlGamma(l+px, j);
                }
              }else{
                for (int l = 0; l < pz; l++) {
                  xzComponent += Zmat(l, (i * Tp + t)) * UwlGamma(l+Tp*px, j);
                }
              }
            }
          }

          if(!gamma_fix){
            gammaComponent=Gamma(t, j);
          }else{
            gammaComponent=Gamma(0,j)*(t+1);
          }

          double logisticFunction = 1 / (1 + exp(-1 * (gammaComponent + xzComponent + thetaComponent)));
          const double eps = 1e-22;
          if (Ymatlog(i, (j * Tp + t))) {
            if(logisticFunction < eps){
              objectiveValue += log(eps);
            }else{
              objectiveValue += log(logisticFunction);
            }
          } else {
            if( (1-logisticFunction)<eps){
              objectiveValue += log(eps);
            }else{
              objectiveValue += log(1 - logisticFunction);
            }

          }
        }
      }
    }
  }

  // Normalize the output
  objectiveValue /= (N * J);
  return objectiveValue;
}


//Add output of indices where covariate in non-zero
//The dimensions and interpretation of the arguments are consistent with the R function lvhml_est_fit_K.
// [[Rcpp::export]]
List rcpp_getEforjEujXiM0(NumericVector uj, IntegerVector indices_R, IntegerVector lengths_R, NumericMatrix Theta, NumericMatrix X,NumericMatrix Zmat, int N, int Tp, int K, int px,int pz, int ujlen, bool ext, bool gamma_fix) {

  int totalObservations = indices_R.length();  // Total number of observations
  //Get number of rows of Gamma.
  int g_len = 0;
  if(!gamma_fix){
    g_len = Tp;
  }else{
    g_len = 1;
  }
  NumericMatrix nonzero_ind;
  nonzero_ind = NumericMatrix(totalObservations,1 + px +pz+K);
  int ujLengthNoD = ujlen - g_len;  // Length of uj excluding gamma part.
  NumericMatrix Eforj(totalObservations, ujlen);  // Matrix E for j
  NumericVector Euj(totalObservations);  // Vector Euj
  NumericVector XiM0(totalObservations);  // Vector XiM0
  int count = 0;

  // Iterate over all subjects
  for (int i = 0; i < N; ++i) {
    // Iterate over the length of R for each subject
    for (int q = 0; q < lengths_R[i]; ++q) {
      //Initiate the vector for additional terms.
      NumericVector additionalTerms(ujLengthNoD);
      int r = indices_R[count];  // Calculate r

      // Add Theta values to additionalTerms
      if(!ext){
        for (int k = 0; k < K; ++k) {
          additionalTerms[px+pz + k] = Theta(i, k);
          nonzero_ind(count,1+px+pz+k) = g_len+px+pz+k;
        }
      }else{
        for (int k = 0; k < K; ++k) {
          additionalTerms[Tp*px+pz +r*K+ k] = Theta(i, k);
          nonzero_ind(count,1+px+pz+k) = g_len+Tp*px+pz+r*K+k;
        }
      }



      // Add X and Zmat values to additionalTerms if px+pz > 0
      if ( (px+pz) > 0) {
        //add x first
        if(px >0){
          if(!ext){
            for (int k = 0; k < px; ++k) {
              additionalTerms[k] = X(i,k) ;
              nonzero_ind(count,1+k) = g_len+k;
            }
          }else{
            for (int k = 0; k < px; ++k) {
              additionalTerms[r*px+k] = X(i,k) ;
              nonzero_ind(count,1+k) = g_len+r*px+k;
            }
          }
        }


        if(pz >0){
          if(!ext){
            for (int k = 0; k < pz; ++k) {
              additionalTerms[k+ px] = Zmat(k, i * Tp + r) ;
              nonzero_ind(count,1+px+k) = g_len+px+k;
            }
          }else{
            for (int k = 0; k < pz; ++k) {
              additionalTerms[Tp*px+k] = Zmat(k, i * Tp + r) ;
              nonzero_ind(count,1+px+k) = g_len+Tp*px+k;
            }
          }


        }
      }
      //Remark: No need to handle the case when px=pz=0, there will be no terms corresponding to it in this case.

      // Update Eforj matrix
      if(!gamma_fix){
        Eforj(count, r) = 1;
        nonzero_ind(count,0) = r;
      }else{
        Eforj(count, 0) = r+1;
        nonzero_ind(count,0)=0;
      }
      for (int l = g_len; l < ujlen; ++l) {
        Eforj(count, l) = additionalTerms[l - g_len];
      }

      // Compute Euj using inner product
      if(!gamma_fix){
        Euj[count] = sum(additionalTerms * uj[Range(g_len, ujlen)]) + uj[r];
      }else{
        Euj[count] = sum(additionalTerms * uj[Range(g_len, ujlen)]) + uj[0]*(r+1);
      }

      count++;
    }
  }

  // Apply xi function on Euj
  XiM0 = rcpp_xi_func(Euj);

  // Return results as a list
  return List::create(
    Named("Eforj") = Eforj,
    Named("Euj") = Euj,
    Named("XiM0") = XiM0,
    Named("nonzero_ind") =nonzero_ind
  );
}

//Projection function
//out; parameter estimate
//K: latent dimension
//proj_const: projection constant.
// [[Rcpp::export]]
List rcpp_proj_func(NumericVector out, int K, double proj_const){
  double constant = proj_const*sqrt(K);
  double value = sqrt(sum(pow(out, 2)));
  int proj_count=0;
  if (value > constant){
    proj_count = 1;
    out = out * constant / value;
  }
  //Return result as a list
  return  List::create(
    Named("out") = out,
    Named("proj_count") = proj_count
  );
}




// Matrix and vector multiplication
// mat and vec can be any matrix and vector as long as vec.length() is the same as mat.ncol().
NumericVector rcpp_matVecMultiply(NumericMatrix mat, NumericVector vec) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  NumericVector out(nrow);
  for(int i = 0; i < nrow; i++) {
    double total = 0;
    for(int j = 0; j < ncol; j++) {
    if(mat(i,j)==0){
      //Skip if mat(i,j)=0
    }else{
      total += mat(i,j) * vec[j];
    }
    }
    out[i] = total;
  }

  return out;
}


// Line search function to update parameter.
//
// Arguments:
//   XiM0:   NumericVector of length len with fitted probabilities at oldpar.
//   yvec: NumericVector of length len with binary responses.
//   grad:NumericVector of length par_len with the gradient at oldpar.
//   oldpar:NumericVector of length par_len with the estimates.
//   drt: NumericVector of length par_len with the search direction.
//   Covariates:NumericMatrix (len x par_len).
//   tooff: NumericVector of length len with the offset term.
//   proj_const: projection constant
// [[Rcpp::export]]
List rcpp_getnewpar_func(NumericVector XiM0, NumericVector yvec, int len, NumericVector grad, NumericVector drt, NumericVector oldpar, NumericMatrix Covariates, NumericVector tooff, double proj_const) {
  // Parameter initialization for the optimization process
  double alpha = 1;
  double zeta = 0.5;
  double chi = 0.9;
  int par_len = oldpar.length();
  double gradSquaredSum = sum(grad * grad);
  double oldObjective = rcpp_objpar_prep(XiM0, yvec, len);
  int iteration = 1;
  double newObjective = -10000;
  NumericVector newPar(par_len);

  // Iterative process to find new parameters
  while (iteration < 30 && newObjective < oldObjective + chi * alpha * gradSquaredSum / zeta) {
    newPar = oldpar + alpha * drt;
    NumericVector XiMnew = rcpp_xi_func(rcpp_matVecMultiply(Covariates, newPar) + tooff);
    newObjective = rcpp_objpar_prep(XiMnew, yvec, len);

    if (Rcpp::traits::is_infinite<REALSXP>(newObjective) || Rcpp::traits::is_na<REALSXP>(newObjective)) {
      newObjective = -10000;
    }
    alpha *= zeta;
    iteration++;
  }

  // Final check and use old parameters if no improvement
  List proj_out  = rcpp_proj_func(newPar, oldpar.size(),proj_const);
  NumericVector finalPar = proj_out["out"];
  int proj_count = proj_out["proj_count"];
  NumericVector XiMfinal = rcpp_xi_func(rcpp_matVecMultiply(Covariates, finalPar) + tooff);
  newObjective = rcpp_objpar_prep(XiMfinal, yvec, len);

  if (newObjective <= oldObjective) {
    finalPar = oldpar;
  }

  return List::create(
    Named("finalPar") = finalPar,
    Named("proj_count") = proj_count
  );
}


// Matrix and vector multiplication with given nonzeroindices to speed up computations relevant to uj.
//Arguments mat and vec are the same as rcpp_matVecMultiply. nonzero_ind has the same number of rows as mat, but is expected to have fewer columns.
NumericVector rcpp_matVecMultiply_foru(NumericMatrix mat, NumericVector vec,NumericMatrix nonzero_ind ) {
  int nrow = mat.nrow();
  int ncol = nonzero_ind.ncol();
  int index;
  NumericVector out(nrow);
  for(int i = 0; i < nrow; i++) {
    double total = 0;
    for(int j = 0; j < ncol; j++) {
      index= nonzero_ind(i,j);
      total += mat(i,index) * vec[index];
    }
    out[i] = total;
  }
  return out;
}


// Line search function to update parameter.
//
// Arguments:
//   nonzero_ind: NumericMatrix of size len x ut_len.
//The interpretations and dimensions of the rest of the arguments are consistent with rcpp_getnewpar_func.
// [[Rcpp::export]]
List rcpp_getnewpar_func_foru(NumericVector XiM0, NumericVector yvec, int len, NumericVector grad, NumericVector drt, NumericVector oldpar, NumericMatrix Covariates, NumericVector tooff,NumericMatrix nonzero_ind, double proj_const) {
  // Parameter initialization for the optimization process
  double alpha = 1;
  double zeta = 0.5;
  double chi = 0.9;
  int par_len = oldpar.length();
  double gradSquaredSum = sum(grad * grad);
  double oldObjective = rcpp_objpar_prep(XiM0, yvec, len);
  int iteration = 1;
  double newObjective = -10000;
  NumericVector newPar(par_len);

  // Iterative process to find new parameters
  while (iteration < 30 && newObjective < oldObjective + chi * alpha * gradSquaredSum / zeta) {
    newPar = oldpar + alpha * drt;
    NumericVector XiMnew = rcpp_xi_func(rcpp_matVecMultiply_foru(Covariates, newPar,nonzero_ind ) + tooff);
    newObjective = rcpp_objpar_prep(XiMnew, yvec, len);

    if (Rcpp::traits::is_infinite<REALSXP>(newObjective) || Rcpp::traits::is_na<REALSXP>(newObjective)) {
      newObjective = -10000;
    }
    alpha *= zeta;
    iteration++;
  }

  // Final check and use old parameters if no improvement
  List proj_out  = rcpp_proj_func(newPar, oldpar.size(), proj_const);
  NumericVector finalPar = proj_out["out"];
  int proj_count = proj_out["proj_count"];
  NumericVector XiMfinal = rcpp_xi_func(rcpp_matVecMultiply_foru(Covariates, finalPar,nonzero_ind) + tooff);
  newObjective = rcpp_objpar_prep(XiMfinal, yvec, len);

  if (newObjective <= oldObjective) {
    finalPar = oldpar;
  }

  return List::create(
    Named("finalPar") = finalPar,
    Named("proj_count") = proj_count
  );
}


//Function to prepare E_J for initial Betahat under SVD approach
//The dimensions and interpretation of the arguments are consistent with the R function lvhml_est_fit_K.
// [[Rcpp::export]]
List rcpp_BetagetEforj( IntegerVector indices_R,IntegerVector lengths_R , NumericMatrix X, NumericMatrix Zmat, int N, int Tp, int px, int pz, bool ext, bool gamma_fix) {

  int sumW =indices_R.length()  ;
//initialise Eforj, based on whether ext is TRUE.
int b_len;
int g_len;
if(ext){
   b_len =px*Tp;
}else{
  b_len =px;
}
if(gamma_fix){
  g_len = 1;
}else{
  g_len = 0;
}
NumericMatrix Eforj(sumW,g_len+ b_len + pz);

    int count = 0;
  for (int i = 0; i < N; ++i) {
    for (int q = 0; q < lengths_R[i]; ++q) {
      //calculate r
      int r = indices_R[count] ;

      if(gamma_fix){
        Eforj(count,0) = r+1;
      }

      if(ext){
        if(px>0){
          for(int k=0; k<px; ++k){
            Eforj(count,g_len+ r*px+ k) = X(i,k);
          }
        }

        if(pz>0){
          for(int k=0; k<pz; ++k){
            Eforj(count,g_len + Tp*px+ k) = Zmat(k, i*Tp + r);
          }
        }

      }else{
        //Add E corresponding to X
        if(px>0){
          for(int k=0; k<px; ++k){
            Eforj(count,g_len +  k) = X(i,k);
          }
        }

        //Add E corresponding to Z
        if(pz >0){
          for(int k=0; k<pz; ++k){
            Eforj(count,g_len + px+ k) = Zmat(k, i*Tp + r);
          }
        }
      }

      count++;
    }
  }
  return List::create(Named("Eforj") = Eforj);
}

// Construct the design matrix Eforj for a single event j in the
// Beta-only model (Gamma and X terms only).
//The dimensions and interpretation of the arguments are consistent with the R function lvhml_est_fit_K.
// [[Rcpp::export]]
List rcpp_BetaonlygetEforj( IntegerVector indices_R,IntegerVector lengths_R , NumericMatrix X, int N, int Tp, int px,bool ext,bool gamma_fix) {

  int sumW =indices_R.length();
  int g_len;
  int b_len;
  if(!gamma_fix){
    g_len = Tp;
  }else{
    g_len = 1;
  }
  if(!ext){
    b_len = px;
  }else{
    b_len = Tp*px;
  }
  int uj_len = g_len+b_len;
  NumericMatrix Eforj(sumW, uj_len);

  int count = 0;
  for (int i = 0; i < N; ++i) {
    for (int q = 0; q < lengths_R[i]; ++q) {
      //calculate r
      int r = indices_R[count] ;
      if(!gamma_fix){
        //set D
        Eforj(count, r) = 1;
      }else{
        Eforj(count, 0) = r+1;
      }

      for(int k=0; k < px; ++k){
        if(!ext){
          Eforj(count, g_len+k) = X(i,k);
        }else{
          Eforj(count, g_len+r*px+ k) =X(i,k);
        }
      }

      count++;
    }
  }
  return List::create(Named("Eforj") = Eforj);
}

