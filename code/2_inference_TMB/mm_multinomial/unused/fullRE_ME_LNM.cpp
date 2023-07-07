
// full RE
// additionally, scaled covariance matrix for MVN in LNM
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{

  Type nll;           // Define variable that holds the return value (neg. log. lik)
  nll = Type(0.0);    // Assign value 1.2; a cast is needed.

  DATA_MATRIX(Y); // observations (count matrix)
  int d = Y.cols();
  int n = Y.rows();
  // UNSTRUCTURED_CORR<Type> Sigma(theta);
  DATA_MATRIX(x); // matrix of covariates for fixed effects
  DATA_MATRIX(z); // matrix for random effects
  int num_individuals = z.cols(); // number of patients
  PARAMETER_MATRIX(beta); // coefficients for the fixed effects
  PARAMETER_MATRIX(u_large); // coefficients for the random effects. Even though it is defined as matrix (for TMB matrix multiplication), it is a vector
  PARAMETER_VECTOR(cov_par_RE); // RE
  PARAMETER_VECTOR(logs_sd_RE); // RE
  PARAMETER_VECTOR(log_sd); // log of the standard devs (scaling factor) for MVN
  // PARAMETER_MATRIX(theta_prime);
  PARAMETER_VECTOR(cov_par); //vector<Type> cov_par((d_min1*d_min1-d_min1)/2);

  int d_min1 = d - 1;
  matrix<Type> theta_prime(n,d_min1); // The probabilities of each event (in ALR)

  // PARAMETER_VECTOR(cov_par); // number of parameters in the covariance matrix to infer
  // vector<Type> cov_par((d*d-d)/2); // how I had it before but I think it was wrong
  using namespace density;
  // Type log_sd; // log of the standard devs (scaling factor) for MVN
  UNSTRUCTURED_CORR_t<Type> nll_mvn(cov_par);
  UNSTRUCTURED_CORR_t<Type> nll_mvn_RE(cov_par_RE);
  // int log_sd = 1; // dummy
  // density::VECSCALE_t<density::UNSTRUCTURED_CORR_t<Type>> scnldens = density::VECSCALE(nll_mvn, log_sd);
  matrix<Type> mu(n,d_min1); // mean parameter for multivariate normal
  vector<Type> residual(d_min1);
  vector<Type> Q(n); // The probabilities of each event (marginal of exp of ALR)
  matrix<Type> theta(n,d); // The probabilities of each event

  for(int i=0;i<num_individuals;i++){ // likelihood for the random effects (multivariate normal)
    nll += VECSCALE_t(nll_mvn_RE, exp(logs_sd_RE))(u_large.row(i));
  }

  mu = x * beta + z * u_large;

  using namespace density;
  for(int l=0;l<n;l++){
    residual = vector<Type>(theta_prime.row(l)) - vector<Type>(mu.row(l));
    // nll += MVNORM(Sigma)(residual); // not used here
    nll += VECSCALE_t(nll_mvn, exp(log_sd))(residual);
  }
  
  for(int l=0;l<n;l++){
    Q(l) = 0;
    for(int j=0;j<d_min1;j++){
      Q(l) += exp(theta_prime(l,j));
    }
  }

  for(int l=0;l<n;l++){
    for(int j=0;j<d_min1;j++){
      theta(l,j) = exp(theta_prime(l,j))/(Type(1.0)+Q(l));
    }
    theta(l,d_min1) = Type(1.0)/(Type(1.0)+Q(l));
  }

  REPORT(mu) // mean of MVN 
  REPORT(cov_par) // covariance matrix of MVN

  for(int l=0; l<n; l++){ // Multinomial draws
    nll -= (dmultinom(vector<Type>(Y.row(l)), vector<Type>(theta.row(l)), true));
  }
  return nll;

}
