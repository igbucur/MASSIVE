#include <RcppArmadillo.h>
#include <boost/math/distributions/normal.hpp>
#include <assert.h>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Preprocessor directives for efficiently computing normal prior
#define LOG_PDF_NORM_VEC(p, sd) (- 0.5 * log(2 * M_PI) - log(sd) - (p % p) / (2 * sd % sd))
#define LOG_PDF_NORM(p, sd) (- 0.5 * log(2 * M_PI) - log(sd) - (p * p) / (2 * sd * sd))

// C++ Structure containing parameters
struct Parameters {
  arma::vec sgamma;
  arma::vec salpha;
  double sbeta;
  double skappa_X;
  double skappa_Y;
  double sigma_X;
  double sigma_Y;
};


// Function for reading Rcpp list of parameters.
struct Parameters read_parameter_list(List param_list) {
  
  struct Parameters params;
  
  params.sgamma = as<arma::vec>(param_list["sgamma"]);
  params.salpha = as<arma::vec>(param_list["salpha"]);
  params.sbeta = as<double>(param_list["sbeta"]);
  params.skappa_X = as<double>(param_list["skappa_X"]);
  params.skappa_Y = as<double>(param_list["skappa_Y"]);
  params.sigma_X = as<double>(param_list["sigma_X"]);
  params.sigma_Y = as<double>(param_list["sigma_Y"]);
  
  return params;
}

// Function for reading Rcpp list of prior standard deviations.
struct Parameters read_prior_list(List param_list) {
  
  struct Parameters params;
  
  params.sgamma = as<arma::vec>(param_list["sgamma"]);
  params.salpha = as<arma::vec>(param_list["salpha"]);
  params.sbeta = as<double>(param_list["sbeta"]);
  params.skappa_X = as<double>(param_list["skappa_X"]);
  params.skappa_Y = as<double>(param_list["skappa_Y"]);
  
  return params;
}

//' Rcpp routine for computing the negative log-prior distribution of the MASSIVE model.
//' 
//' @name neg_log_prior
//'
//' @param param_list List of parameter values in Rcpp format.
//' @param prior_sd List of prior standard deviations in Rcpp format.
//'
//' @return numeric value; negative log-prior value for the MASSIVE model.
//' 
//' @keywords internal
double neg_log_prior(List param_list, List prior_sd) {
  
  struct Parameters par = read_parameter_list(param_list);
  struct Parameters sd = read_prior_list(prior_sd);
  
  return - (sum(LOG_PDF_NORM_VEC(par.sgamma, sd.sgamma) + LOG_PDF_NORM_VEC(par.salpha, sd.salpha)) +
            LOG_PDF_NORM(par.sbeta, 10) + LOG_PDF_NORM(par.skappa_X, 10) + LOG_PDF_NORM(par.skappa_Y, 10) - log(par.sigma_X) - log(par.sigma_Y));
}


//' Rcpp routine for computing the negative log-posterior distribution of the MASSIVE model.
//'
//' @param J Integer number of genetic instrumental variables.
//' @param N Integer number of observations.
//' @param SS Numeric matrix containing first- and second-order statistics.
//' @param sigma_G Numeric vector of genetic IV standard deviations.
//' @param param_list List of IV model parameter values.
//' @param prior_sd List of standard deviations for the parameter Gaussian priors.
//' @param n Integer number of alleles (trials) for the binomial genetic variants.
//'
//' @return numeric value; negative log-posterior value for the MASSIVE model.
// [[Rcpp::export]]
double Rcpp_scaled_neg_log_posterior(
    unsigned J, unsigned N, arma::mat SS, arma::vec sigma_G, 
    List param_list, List prior_sd, unsigned n = 2) {
  
  arma::mat f(J + 3, 2, arma::fill::zeros);
  arma::mat sigma(2, 2, arma::fill::zeros);
  
  struct Parameters par = read_parameter_list(param_list);
  struct Parameters sd = read_prior_list(prior_sd);
  
  for (unsigned j = 0; j < J; ++j) {
    f(j+1, 0) = - par.sigma_X * par.sgamma(j) / sigma_G(j);
    f(j+1, 1) = - par.sigma_Y * (par.salpha(j) + par.sbeta * par.sgamma(j)) / sigma_G(j);
  }
  f(J+1, 0) = f(J+2, 1) = 1;
  
  sigma(0, 0) = par.sigma_X * par.sigma_X * (1 + par.skappa_X * par.skappa_X);
  sigma(0, 1) = sigma(1, 0) = par.sigma_X * par.sigma_Y * (par.sbeta * (1 + par.skappa_X * par.skappa_X) + par.skappa_X * par.skappa_Y);
  sigma(1, 1) = par.sigma_Y * par.sigma_Y * (1 + par.sbeta * par.sbeta + (par.skappa_Y + par.sbeta * par.skappa_X) * (par.skappa_Y + par.sbeta * par.skappa_X));
  
  arma::mat S = f.t() * SS * f;
  
  double val;
  double sign;
  
  log_det(val, sign, sigma);
  
  // Check for positive determinant of sigma in Rcpp_scaled_neg_log_posterior.
  assert(sign > 0);
  
  return 0.5 * (log(4 * M_PI * M_PI) + val + trace(inv(sigma) * S)) + 
    neg_log_prior(param_list, prior_sd) / N;
}


//' Rcpp routine for computing the negative log-posterior gradient of the MASSIVE model.
//'
//' @param J Integer number of genetic instrumental variables.
//' @param N Integer number of observations.
//' @param SS Numeric matrix containing first- and second-order statistics.
//' @param sigma_G Numeric vector of genetic IV standard deviations.
//' @param param_list List of IV model parameter values.
//' @param prior_sd List of standard deviations for the parameter Gaussian priors.
//' @param n Integer number of alleles (trials) for the binomial genetic variants.
//'
//' @return numeric vector; negative log-posterior gradient for the MASSIVE model.
// [[Rcpp::export]]
arma::vec Rcpp_scaled_neg_log_gradient(unsigned J, unsigned N, arma::mat SS, arma::vec sigma_G, List param_list, List prior_sd, unsigned n = 2) {
  
  arma::vec gradient(2 * J + 5, arma::fill::zeros);
  arma::mat sf(J + 3, 2, arma::fill::zeros);
  arma::mat iscov(2, 2, arma::fill::zeros);
  
  struct Parameters par = read_parameter_list(param_list);
  struct Parameters sd = read_prior_list(prior_sd);
  
  
  for (unsigned j = 0; j < J; ++j) {
    sf(j+1, 0) = - par.sgamma(j);
    sf(j+1, 1) = - (par.salpha(j) + par.sbeta * par.sgamma(j));
  }
  sf(J+1, 0) = sf(J+2, 1) = 1;
  
  iscov(0, 0) = 1 + par.sbeta * par.sbeta + (par.skappa_Y + par.sbeta * par.skappa_X) * (par.skappa_Y + par.sbeta * par.skappa_X);
  iscov(0, 1) = iscov(1, 0) = - (par.sbeta * (1 + par.skappa_X * par.skappa_X) + par.skappa_X * par.skappa_Y);
  iscov(1, 1) = 1 + par.skappa_X * par.skappa_X;
  
  iscov = iscov / (1 + par.skappa_X * par.skappa_X + par.skappa_Y * par.skappa_Y);
  
  
  arma::mat ihalf_Sigma_Z(J+3, J+3, arma::fill::zeros);
  ihalf_Sigma_Z(0, 0) = 1;
  for (unsigned j = 0; j < J; ++j) ihalf_Sigma_Z(j+1, j+1) = 1 / sigma_G(j); 
  ihalf_Sigma_Z(J+1, J+1) = 1 / par.sigma_X;
  ihalf_Sigma_Z(J+2, J+2) = 1 / par.sigma_Y;
  
  arma::mat common_term = ihalf_Sigma_Z * SS * ihalf_Sigma_Z * sf * iscov;
  arma::mat I_2(2, 2, arma::fill::eye);
  
  // derivatives w.r.t sgamma and salpha
  for (unsigned j = 0; j < J; ++j) {
    gradient(j) = - (common_term(j+1, 0) + par.sbeta * common_term(j+1, 1)) + par.sgamma(j) / (sd.sgamma(j) * sd.sgamma(j) * N);
    gradient(j+J) = - common_term(j+1, 1) + par.salpha(j) / (sd.salpha(j) * sd.salpha(j) * N);
  }
  
  // derivative w.r.t sbeta
  arma::mat d_sbeta_scov(2, 2, arma::fill::zeros);
  d_sbeta_scov(0, 1) = d_sbeta_scov(1, 0) = 1 + par.skappa_X * par.skappa_X;
  d_sbeta_scov(1, 1) = 2 * par.sbeta * (1 + par.skappa_X * par.skappa_X) + 2 * par.skappa_X * par.skappa_Y;
  gradient(2*J) = - 0.5 * trace(d_sbeta_scov * iscov * (sf.t() * common_term - I_2)) + par.sbeta / (100 * N);
  for (unsigned j = 0; j < J; ++j) gradient(2*J) -= common_term(j+1, 1) * par.sgamma(j);
  
  arma::mat d_skappa_X_scov(2, 2, arma::fill::zeros);
  d_skappa_X_scov(0, 0) = 2 * par.skappa_X;
  d_skappa_X_scov(0, 1) = d_skappa_X_scov(1, 0) = 2 * par.sbeta * par.skappa_X + par.skappa_Y;
  d_skappa_X_scov(1, 1) = 2 * par.sbeta * (par.skappa_Y + par.sbeta * par.skappa_X);
  gradient(2*J+1) = - 0.5 * trace(d_skappa_X_scov * iscov * (sf.t() * common_term - I_2)) + par.skappa_X / (100 * N);
  
  arma::mat d_skappa_Y_scov(2, 2, arma::fill::zeros);
  d_skappa_Y_scov(0, 1) = d_skappa_Y_scov(1, 0) = par.skappa_X;
  d_skappa_Y_scov(1, 1) = 2 * (par.skappa_Y + par.sbeta * par.skappa_X);
  gradient(2*J+2) = - 0.5 * trace(d_skappa_Y_scov * iscov * (sf.t() * common_term - I_2)) + par.skappa_Y / (100 * N);
  
  // Huge simplification can be obtained here
  gradient(2*J+3) = - common_term(J+1, 0) / par.sigma_X + (N + 1.0) / (par.sigma_X * N); // + par.sigma_X / (100 * N);
  gradient(2*J+4) = - common_term(J+2, 1) / par.sigma_Y + (N + 1.0) / (par.sigma_Y * N); // + par.sigma_Y / (100 * N);
  
  return gradient;
}


//' Rcpp routine for computing the negative log-posterior Hessian of the MASSIVE model.
//'
//' @param J Integer number of genetic instrumental variables.
//' @param N Integer number of observations.
//' @param SS Numeric matrix containing first- and second-order statistics.
//' @param sigma_G Numeric vector of genetic IV standard deviations.
//' @param param_list List of IV model parameter values.
//' @param prior_sd List of standard deviations for the parameter Gaussian priors.
//' @param n Integer number of alleles (trials) for the binomial genetic variants.
//'
//' @return numeric matrix; negative log-posterior Hessian for the MASSIVE model.
// [[Rcpp::export]]
arma::mat Rcpp_scaled_neg_log_hessian(unsigned J, unsigned N, arma::mat SS, arma::vec sigma_G, List param_list, List prior_sd, unsigned n = 2) {
  
  struct Parameters par = read_parameter_list(param_list);
  struct Parameters sd = read_prior_list(prior_sd);
  
  arma::vec gradient(2 * J + 5, arma::fill::zeros);
  arma::mat hessian(2 * J + 5, 2 * J + 5, arma::fill::zeros);
  
  arma::mat sf(J + 3, 2, arma::fill::zeros);
  arma::mat iscov(2, 2, arma::fill::zeros);
  
  for (unsigned j = 0; j < J; ++j) {
    sf(j+1, 0) = - par.sgamma(j);
    sf(j+1, 1) = - (par.salpha(j) + par.sbeta * par.sgamma(j));
  }
  sf(J+1, 0) = sf(J+2, 1) = 1;
  
  
  iscov(0, 0) = 1 + par.sbeta * par.sbeta + (par.skappa_Y + par.sbeta * par.skappa_X) * (par.skappa_Y + par.sbeta * par.skappa_X);
  iscov(0, 1) = iscov(1, 0) = - (par.sbeta * (1 + par.skappa_X * par.skappa_X) + par.skappa_X * par.skappa_Y);
  iscov(1, 1) = 1 + par.skappa_X * par.skappa_X;
  
  iscov = iscov / (1 + par.skappa_X * par.skappa_X + par.skappa_Y * par.skappa_Y);
  
  
  arma::mat ihalf_Sigma_Z(J+3, J+3, arma::fill::zeros);
  ihalf_Sigma_Z(0, 0) = 1;
  for (unsigned j = 0; j < J; ++j) ihalf_Sigma_Z(j+1, j+1) = 1 / sigma_G(j); 
  ihalf_Sigma_Z(J+1, J+1) = 1 / par.sigma_X;
  ihalf_Sigma_Z(J+2, J+2) = 1 / par.sigma_Y;
  
  arma::mat I_2(2, 2, arma::fill::eye);
  
  arma::mat ihalf_Sigma_X(J+3, J+3, arma::fill::zeros); ihalf_Sigma_X(J+1, J+1) = ihalf_Sigma_Z(J+1, J+1);
  arma::mat ihalf_Sigma_Y(J+3, J+3, arma::fill::zeros); ihalf_Sigma_Y(J+2, J+2) = ihalf_Sigma_Z(J+2, J+2);
  
  arma::mat sSS = ihalf_Sigma_Z * SS * ihalf_Sigma_Z;
  arma::mat common_term = sSS * sf * iscov;
  
  arma::cube d_sgamma_sf(J+3, 2, J, arma::fill::zeros); // derivative of sf w.r.t. sgamma
  arma::cube d_salpha_sf(J+3, 2, J, arma::fill::zeros); // derivative of sf w.r.t. salpha
  arma::mat d_sbeta_sf(J+3, 2, arma::fill::zeros);
  
  for (unsigned j = 0; j < J; ++j) {
    d_sgamma_sf(j+1, 0, j) = -1;
    d_sgamma_sf(j+1, 1, j) = -par.sbeta;
    d_salpha_sf(j+1, 1, j) = -1;
    d_sbeta_sf(j+1, 1) = -par.sgamma(j);
  }
  
  /* 
   * Derivatives of scov
   */
  arma::mat d_sbeta_scov(2, 2, arma::fill::zeros);
  d_sbeta_scov(0, 1) = d_sbeta_scov(1, 0) = 1 + par.skappa_X * par.skappa_X;
  d_sbeta_scov(1, 1) = 2 * par.sbeta * (1 + par.skappa_X * par.skappa_X) + 2 * par.skappa_X * par.skappa_Y;
  
  arma::mat d_skappa_X_scov(2, 2, arma::fill::zeros);
  d_skappa_X_scov(0, 0) = 2 * par.skappa_X;
  d_skappa_X_scov(0, 1) = d_skappa_X_scov(1, 0) = 2 * par.sbeta * par.skappa_X + par.skappa_Y;
  d_skappa_X_scov(1, 1) = 2 * par.sbeta * (par.skappa_Y + par.sbeta * par.skappa_X);
  
  arma::mat d_skappa_Y_scov(2, 2, arma::fill::zeros);
  d_skappa_Y_scov(0, 1) = d_skappa_Y_scov(1, 0) = par.skappa_X;
  d_skappa_Y_scov(1, 1) = 2 * (par.skappa_Y + par.sbeta * par.skappa_X);
  
  arma::mat d2_sbeta_scov(2, 2, arma::fill::zeros);
  d2_sbeta_scov(1, 1) = 2 * (1 + par.skappa_X * par.skappa_X);
  
  arma::mat d2_skappa_X_scov(2, 2, arma::fill::zeros);
  d2_skappa_X_scov(0, 0) = 2;
  d2_skappa_X_scov(0, 1) = d2_skappa_X_scov(1, 0) = 2 * par.sbeta;
  d2_skappa_X_scov(1, 1) = 2 * par.sbeta * par.sbeta;
  
  arma::mat d2_skappa_Y_scov(2, 2, arma::fill::zeros); d2_skappa_Y_scov(1, 1) = 2;
  
  arma::mat d2_sbeta_skappa_X_scov(2, 2, arma::fill::zeros);
  d2_sbeta_skappa_X_scov(0, 1) = d2_sbeta_skappa_X_scov(1, 0) = 2 * par.skappa_X;
  d2_sbeta_skappa_X_scov(1, 1) = 4 * par.sbeta * par.skappa_X + 2 * par.skappa_Y;
  
  arma::mat d2_sbeta_skappa_Y_scov(2, 2, arma::fill::zeros); d2_sbeta_skappa_Y_scov(1, 1) = 2 * par.skappa_X;
  
  arma::mat d2_skappa_X_skappa_Y_scov(2, 2, arma::fill::zeros); 
  d2_skappa_X_skappa_Y_scov(0, 1) = d2_skappa_X_skappa_Y_scov(1, 0) = 1;
  d2_skappa_X_skappa_Y_scov(1, 1) = 2 * par.sbeta;
  
  /*
   * Derivatives of Sigma_Z
   */
  arma::mat d_sigma_X_Sigma_Z(J + 3, J + 3, arma::fill::zeros); d_sigma_X_Sigma_Z(J + 1, J + 1) = 1;
  arma::mat d_sigma_Y_Sigma_Z(J + 3, J + 3, arma::fill::zeros); d_sigma_Y_Sigma_Z(J + 2, J + 2) = 1;
  
  
  for (unsigned j = 0; j < J; ++j) {
    for (unsigned k = 0; k < J; ++k) {
      hessian(j, k) = trace(iscov * d_sgamma_sf.slice(j).t() * sSS * d_sgamma_sf.slice(k)) + double(j == k) / (sd.sgamma(j) * sd.sgamma(j) * N);
      hessian(j, J + k) = hessian(J + k, j) = trace(iscov * d_sgamma_sf.slice(j).t() * sSS * d_salpha_sf.slice(k));
      hessian(J + j, J + k) = trace(iscov * d_salpha_sf.slice(j).t() * sSS * d_salpha_sf.slice(k)) + double(j == k) / (sd.salpha(j) * sd.salpha(j) * N);
    }
    
    hessian(j, 2 * J) = hessian(2 * J, j) = - trace(d_sbeta_scov * iscov * d_sgamma_sf.slice(j).t() * common_term) + trace(iscov * d_sgamma_sf.slice(j).t() * sSS * d_sbeta_sf) + trace(d_salpha_sf.slice(j).t() * common_term);
    hessian(j, 2 * J + 1) = hessian(2 * J + 1, j) = - trace(d_skappa_X_scov * iscov * d_sgamma_sf.slice(j).t() * common_term);
    hessian(j, 2 * J + 2) = hessian(2 * J + 2, j) = - trace(d_skappa_Y_scov * iscov * d_sgamma_sf.slice(j).t() * common_term);
    hessian(j, 2 * J + 3) = hessian(2 * J + 3, j) = - trace(iscov * sf.t() * ihalf_Sigma_X * sSS * d_sgamma_sf.slice(j));
    hessian(j, 2 * J + 4) = hessian(2 * J + 4, j) = - trace(iscov * sf.t() * ihalf_Sigma_Y * sSS * d_sgamma_sf.slice(j));
    hessian(J + j, 2 * J) = hessian(2 * J, J + j) = - trace(d_sbeta_scov * iscov * d_salpha_sf.slice(j).t() * common_term) + trace(iscov * d_salpha_sf.slice(j).t() * sSS * d_sbeta_sf);
    hessian(J + j, 2 * J + 1) = hessian(2 * J + 1, J + j) = - trace(d_skappa_X_scov * iscov * d_salpha_sf.slice(j).t() * common_term);
    hessian(J + j, 2 * J + 2) = hessian(2 * J + 2, J + j) = - trace(d_skappa_Y_scov * iscov * d_salpha_sf.slice(j).t() * common_term);
    hessian(J + j, 2 * J + 3) = hessian(2 * J + 3, J + j) = - trace(iscov * sf.t() * ihalf_Sigma_X * sSS * d_salpha_sf.slice(j));
    hessian(J + j, 2 * J + 4) = hessian(2 * J + 4, J + j) = - trace(iscov * sf.t() * ihalf_Sigma_Y * sSS * d_salpha_sf.slice(j));
    
  }
  
  hessian(2*J, 2*J) = trace(d_sbeta_scov * iscov * d_sbeta_scov * iscov * sf.t() * common_term) -
    trace(d2_sbeta_scov * iscov * sf.t() * common_term) / 2.0 -
    trace(d_sbeta_scov * iscov * d_sbeta_sf.t() * common_term) * 2.0 +
    trace(iscov * d_sbeta_sf.t() * sSS * d_sbeta_sf) -
    trace(iscov * d_sbeta_scov * iscov * d_sbeta_scov) / 2.0 +
    trace(iscov * d2_sbeta_scov) / 2.0 + 1.0 / (100 * N);
  
  hessian(2 * J + 1, 2 * J + 1) = - trace(d2_skappa_X_scov * iscov * sf.t() * common_term) / 2.0 +
    trace(d_skappa_X_scov * iscov * d_skappa_X_scov * iscov * sf.t() * common_term) -
    trace(d_skappa_X_scov * iscov * d_skappa_X_scov * iscov) / 2.0 +
    trace(d2_skappa_X_scov * iscov) / 2.0+ 1.0 / (100 * N);
  
  hessian(2 * J + 2, 2 * J + 2) = - trace(d2_skappa_Y_scov * iscov * sf.t() * common_term) / 2.0 +
    trace(d_skappa_Y_scov * iscov * d_skappa_Y_scov * iscov * sf.t() * common_term) -
    trace(d_skappa_Y_scov * iscov * d_skappa_Y_scov * iscov) / 2.0 +
    trace(d2_skappa_Y_scov * iscov) / 2.0 + 1.0 / (100 * N);
  
  hessian(2 * J + 3, 2 * J + 3) = trace(sf.t() * ihalf_Sigma_Z * d_sigma_X_Sigma_Z * ihalf_Sigma_Z * d_sigma_X_Sigma_Z * common_term) * 2.0 +
    trace(iscov * sf.t() * ihalf_Sigma_Z * d_sigma_X_Sigma_Z * sSS * d_sigma_X_Sigma_Z * ihalf_Sigma_Z * sf) - (N + 1.0) / (par.sigma_X * par.sigma_X * N);
  
  hessian(2 * J + 4, 2 * J + 4) = trace(sf.t() * ihalf_Sigma_Z * d_sigma_Y_Sigma_Z * ihalf_Sigma_Z * d_sigma_Y_Sigma_Z * common_term) * 2.0 +
    trace(iscov * sf.t() * ihalf_Sigma_Z * d_sigma_Y_Sigma_Z * sSS * d_sigma_Y_Sigma_Z * ihalf_Sigma_Z * sf) - (N + 1.0) / (par.sigma_Y * par.sigma_Y * N);
  
  hessian(2 * J, 2 * J + 1) = hessian(2 * J + 1, 2 * J) = 
    - trace(d2_sbeta_skappa_X_scov * iscov * sf.t() * common_term) / 2.0 +
    trace(d_sbeta_scov * iscov * d_skappa_X_scov * iscov * sf.t() * common_term) -
    trace(iscov * d_skappa_X_scov * iscov * d_sbeta_sf.t() * sSS * sf) -
    trace(iscov * d_skappa_X_scov * iscov * d_sbeta_scov) / 2.0 +
    trace(iscov * d2_sbeta_skappa_X_scov) / 2.0;
  
  hessian(2 * J , 2 * J + 2) = hessian(2 * J + 2, 2 * J) = 
    - trace(d2_sbeta_skappa_Y_scov * iscov * sf.t() * common_term) / 2.0 +
    trace(d_sbeta_scov * iscov * d_skappa_Y_scov * iscov * sf.t() * common_term) -
    trace(iscov * d_skappa_Y_scov * iscov * d_sbeta_sf.t() * sSS * sf) -
    trace(iscov * d_skappa_Y_scov * iscov * d_sbeta_scov) / 2.0 +
    trace(iscov * d2_sbeta_skappa_Y_scov) / 2.0;
  
  hessian(2 * J, 2 * J + 3) = hessian(2 * J + 3, 2 * J) =
    trace(d_sbeta_scov * iscov * sf.t() * ihalf_Sigma_Z * d_sigma_X_Sigma_Z * common_term) +
    trace(d_sbeta_sf.t() * ihalf_Sigma_Z * d_sigma_X_Sigma_Z * common_term) -
    trace(iscov * d_sbeta_sf.t() * sSS * d_sigma_X_Sigma_Z * ihalf_Sigma_Z * sf);
  
  hessian(2 * J, 2 * J + 4) = hessian(2 * J + 4, 2 * J) = 
    trace(d_sbeta_scov * iscov * sf.t() * ihalf_Sigma_Z * d_sigma_Y_Sigma_Z * common_term) +
    trace(d_sbeta_sf.t() * ihalf_Sigma_Z * d_sigma_Y_Sigma_Z * common_term) -
    trace(iscov * d_sbeta_sf.t() * sSS * d_sigma_Y_Sigma_Z * ihalf_Sigma_Z * sf);
  
  hessian(2 * J + 1, 2 * J + 2) = hessian(2 * J + 2, 2 * J + 1) = 
    - trace(d2_skappa_X_skappa_Y_scov * iscov * sf.t() * common_term) / 2.0 +
    trace(d_skappa_X_scov * iscov * d_skappa_Y_scov * iscov * sf.t() * common_term) +
    trace(d2_skappa_X_skappa_Y_scov * iscov) / 2.0 -
    trace(d_skappa_X_scov * iscov * d_skappa_Y_scov * iscov) / 2.0;
  
  hessian(2 * J + 1, 2 * J + 3) = hessian(2 * J + 3, 2 * J + 1) =
    trace(d_skappa_X_scov * iscov * sf.t() * ihalf_Sigma_Z * d_sigma_X_Sigma_Z * common_term);
  
  hessian(2 * J + 1, 2 * J + 4) = hessian(2 * J + 4, 2 * J + 1) = 
    trace(d_skappa_X_scov * iscov * sf.t() * ihalf_Sigma_Z * d_sigma_Y_Sigma_Z * common_term);
  
  hessian(2 * J + 2, 2 * J + 3) = hessian(2 * J + 3, 2 * J + 2) =
    trace(d_skappa_Y_scov * iscov * sf.t() * ihalf_Sigma_Z * d_sigma_X_Sigma_Z * common_term);
  
  hessian(2 * J + 2, 2 * J + 4) = hessian(2 * J + 4, 2 * J + 2) =
    trace(d_skappa_Y_scov * iscov * sf.t() * ihalf_Sigma_Z * d_sigma_Y_Sigma_Z * common_term);
  
  hessian(2 * J + 3, 2 * J + 4) = hessian(2 * J + 4, 2 * J + 3) =
    trace(iscov * sf.t() * ihalf_Sigma_Z * d_sigma_X_Sigma_Z * sSS * d_sigma_Y_Sigma_Z * ihalf_Sigma_Z * sf);
  
  return hessian;
}
