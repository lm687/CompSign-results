
//https://figshare.com/articles/C_code_needed_for_running_R_code/4836503
// Simple multinomial fitting.
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{

  Type nll;           // Define variable that holds the return value (neg. log. lik)
  nll = Type(1.2);    // Assign value 1.2; a cast is needed.

  DATA_MATRIX(Y);
  DATA_INTEGER(d);
  DATA_INTEGER(n);
  PARAMETER_VECTOR(theta_prime); // probability vector in ALR
  int d_min1 = d - 1;

  Type Q; // To get the denominator
  Q = Type(0);
  // Q=0;
  for(int j=0;j<d_min1;j++){
    Q += exp(theta_prime(j));
  }

  vector<Type> theta(d); // The probabilities of each event
  for(int j=0;j<d_min1;j++){
    theta(j) = exp(theta_prime(j))/(Type(1.0)+Q);
  }
  theta(d_min1) = Type(1.0)/(Type(1.0)+Q);

  for(int i=0; i<n; i++){ // Multinomial draws
    vector<Type> ith_row = Y.row(i);
    nll = -(dmultinom(ith_row, theta, true));
  }
  return nll;

}
