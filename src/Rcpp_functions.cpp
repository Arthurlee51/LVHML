#include <Rcpp.h>
using namespace Rcpp;

//calculate hessian matrix for uj or thetai:
//Mvec: Euj or eUi, Covariate: Eforj or Afori depending on context, parvec: uj or thetai, len:N or J
// [[Rcpp::export]]
NumericMatrix rcpp_cal_Hes(NumericVector XiM0, NumericMatrix Covariate, NumericVector parvec, int len, int XiM0_len) {
  int par_len = parvec.size();
  NumericMatrix Hes(par_len, par_len);
  //NumericVector commonpart = -XiM0 / (1 + exp(Mvec));
  NumericVector commonpart = -XiM0*(1 -XiM0 ) ;
  
  for (int k = 0; k < par_len; ++k) {
    for (int l = k; l < par_len; ++l) {
      Hes(k, l) = sum(commonpart * Covariate(_, k) * Covariate(_, l)) / len;
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
//yvec: yi or yj depending on context, Covariate: Eforj or Afori depending on context, len: N or J, parlen: length of aj of thetai
// [[Rcpp::export]]
NumericVector rcpp_cal_grad(NumericVector yvec, NumericVector XiM0, NumericMatrix Covariate, int len, int parlen) {
  NumericVector grad(parlen);

  for (int j = 0; j < parlen; ++j) {
    grad[j] = sum((yvec - XiM0) * Covariate(_, j)) / len;
  }
  
  return grad;
}




//Function to compute value of objective function of related to the par uj or thetai:
//yvec: yi or yj. len: N or J
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

// [[Rcpp::export]]
NumericVector rcpp_xi_func(NumericVector x) {
  return 1 / (1 + exp(-x));
}

// [[Rcpp::export]]
double rcpp_objprep(LogicalMatrix Ymatlog, IntegerMatrix R, NumericMatrix tU, NumericMatrix Theta, NumericMatrix XZmat, int N, int J, int P, int K, int Tp) {
  double objectiveValue = 0;
  
  // Extracting submatrices for Gamma and UwlGamma from tU
  NumericMatrix Gamma = tU(Range(0, Tp - 1), Range(0, tU.ncol() - 1));
  NumericMatrix UwlGamma = tU(Range(Tp, P + K + Tp - 1), Range(0, tU.ncol() - 1));
  
  // Iterate over all subjects, time points, and variables
  for (int i = 0; i < N; i++) {
    for (int t = 0; t < Tp; t++) {
      if (R(i, t) == 1) { // Only proceed if R(i, t) is 1
        for (int j = 0; j < J; j++) {
          double thetaComponent = 0;
          for (int k = 0; k < K; k++) {
            thetaComponent += Theta(i, k) * UwlGamma(P + k, j);
          }
          
          double xzComponent = 0;
          if (P > 0) {
            for (int l = 0; l < P; l++) {
              xzComponent += XZmat(l, (i * Tp + t)) * UwlGamma(l, j);
            }
          }
          
          double logisticFunction = 1 / (1 + exp(-1 * (Gamma(t, j) + xzComponent + thetaComponent)));
          if (Ymatlog(i, (j * Tp + t))) {
            objectiveValue += log(logisticFunction);
          } else {
            objectiveValue += log(1 - logisticFunction); 
          }
        }
      }
    }
  }
  
  // Normalize the output
  objectiveValue /= (N * J);
  return objectiveValue;
}


// [[Rcpp::export]]
List rcpp_getEforjEujXiM0(NumericVector uj, IntegerVector indices_R, IntegerVector lengths_R, NumericMatrix Theta, NumericMatrix XZmat, int N, int Tp, int K, int P, int ujlen) {
  
  int totalObservations = indices_R.length();  // Total number of observations
  int ujLengthNoD = ujlen - Tp;  // Length of uj excluding timepoints
  NumericMatrix Eforj(totalObservations, ujlen);  // Matrix E for j
  NumericVector Euj(totalObservations);  // Vector Euj
  NumericVector XiM0(totalObservations);  // Vector XiM0
  int count = 0;
  
  // Iterate over all subjects
  for (int i = 0; i < N; ++i) {
    NumericVector additionalTerms(ujLengthNoD);
    
    // Add Theta values to additionalTerms
    for (int k = 0; k < K; ++k) {
      additionalTerms[P + k] = Theta(i, k);
    }
    
    // Iterate over the length of R for each subject
    for (int q = 0; q < lengths_R[i]; ++q) {
      int r = indices_R[count];  // Calculate r
      
      // Add XZmat values to additionalTerms if P > 0
      if (P > 0) {
        for (int k = 0; k < P; ++k) {
          additionalTerms[k] = XZmat(k, i * Tp + r);
        }
      }
      
      // Update Eforj matrix
      Eforj(count, r) = 1;
      for (int l = Tp; l < ujlen; ++l) {
        Eforj(count, l) = additionalTerms[l - Tp];
      }
      
      // Compute Euj using inner product
      Euj[count] = sum(additionalTerms * uj[Range(Tp, ujlen)]) + uj[r];
      count++;
    }
  }
  
  // Apply xi function on Euj
  XiM0 = rcpp_xi_func(Euj);
  
  // Return results as a list
  return List::create(
    Named("Eforj") = Eforj,
    Named("Euj") = Euj,
    Named("XiM0") = XiM0
  );
}

// [[Rcpp::export]]
NumericVector rcpp_proj_func(NumericVector out, int K){
  double constant = 5*sqrt(K);
  double value = sqrt(sum(pow(out, 2)));
  if (value > constant){
    out = out * constant / value;
  }
  return out;
}

// Matrix and vector multiplication
NumericVector rcpp_matVecMultiply(NumericMatrix mat, NumericVector vec) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  NumericVector out(nrow);
  
  for(int i = 0; i < nrow; i++) {
    double total = 0;
    for(int j = 0; j < ncol; j++) {
      total += mat(i,j) * vec[j];
    }
    out[i] = total;
  }
  
  return out;
}



// [[Rcpp::export]]
NumericVector rcpp_getnewpar_func(NumericVector XiM0, NumericVector yvec, int len, NumericVector grad, NumericVector drt, NumericVector oldpar, NumericMatrix Covariates, NumericVector tooff) {
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
  NumericVector finalPar = rcpp_proj_func(newPar, oldpar.size());
  NumericVector XiMfinal = rcpp_xi_func(rcpp_matVecMultiply(Covariates, finalPar) + tooff);
  newObjective = rcpp_objpar_prep(XiMfinal, yvec, len);
  
  if (newObjective <= oldObjective) {
    finalPar = oldpar;
  }
  
  return finalPar;
}


//Function to prepare E_J for initial Betahat under SVD approach
// [[Rcpp::export]]
List rcpp_BetagetEforj( IntegerVector indices_R,IntegerVector lengths_R , NumericMatrix XZmat, int N, int Tp, int P) {
  
  int sumW =indices_R.length()  ;
  NumericMatrix Eforj(sumW, P);
  
  int count = 0;
  for (int i = 0; i < N; ++i) {
    NumericVector toadd(P);
    for (int q = 0; q < lengths_R[i]; ++q) {
      //calculate r
      int r = indices_R[count] ;
      
      for(int k=0; k<P; ++k){
        Eforj(count, k) = XZmat(k, i*Tp + r);
      }
      
      
      count++;
    }
  }
  return List::create(Named("Eforj") = Eforj);
}




