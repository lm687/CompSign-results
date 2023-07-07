// Simple linear regression.
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Y is the matrix of observations and has dimensions (2n,d) and is a matrix of counts
  // x is the covariate matrix and has dimensions (p, 2n)
  // beta is the coefficient matrix for the groups and has dimensions (p, d-1)
  // Z is the random effects matrix and has dimensions (2n, n)
  // u is the corresponding matrix of random effects coefficients and has shape (n,1)
  // logsdu log of the standard deviation of u

  DATA_MATRIX(Y);
  DATA_MATRIX(x);  // Fixed effect design matrix
  DATA_SPARSE_MATRIX(Z);  // Random effect design matrix
  PARAMETER(beta);
  PARAMETER(u);
  PARAMETER(logSigma);
  PARAMETER(n); // number of samples
  PARAMETER(d); // number of features
  ADREPORT(exp(2*logSigma));



  // Distribution of random effect (u):
  Type ans = 0;
  ans -= dnorm(u, Type(0), exp(logsdu), true).sum();


  vector<Type> thetaprime = x * beta + Z * u;
  matrix<Type> theta;
  for(i = 0, i<(2*n), i++){
    denom = 0
    for(j =0, j < d; j++){
      thetaprime[l,k]
    }

    for(j =0, j < d; j++){
    theta[l,j] = thetaprime[l,j]/denom

    }
    theta[l,d] = 1/denom
  }

  ans -= dnorm(x, y, exp(logsd0), true).sum();

  Type nll = -sum(dnorm(Y, a+b*x, exp(logSigma), true));
  return nll;
}