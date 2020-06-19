#include "MASSIVE_model.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

using namespace Rcpp;

boost::math::normal standard_gaussian(0.0, 1.0);

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

struct Parameters read_prior_list(List param_list) {
  
  struct Parameters params;
  
  params.sgamma = as<arma::vec>(param_list["sgamma"]);
  params.salpha = as<arma::vec>(param_list["salpha"]);
  params.sbeta = as<double>(param_list["sbeta"]);
  params.skappa_X = as<double>(param_list["skappa_X"]);
  params.skappa_Y = as<double>(param_list["skappa_Y"]);
  
  return params;
}

// [[Rcpp::export]]
double neg_log_prior_MR(List param_list, List prior_sd) {
  
  struct Parameters par = read_parameter_list(param_list);
  struct Parameters sd = read_prior_list(prior_sd);
  
  return - (sum(LOG_PDF_NORM_VEC(par.sgamma, sd.sgamma) + LOG_PDF_NORM_VEC(par.salpha, sd.salpha)) +
    LOG_PDF_NORM(par.sbeta, 10) + LOG_PDF_NORM(par.skappa_X, 10) + LOG_PDF_NORM(par.skappa_Y, 10) - log(par.sigma_X) - log(par.sigma_Y));
    //LOG_PDF_NORM(log(par.sigma_X), 10) + LOG_PDF_NORM(log(par.sigma_Y), 10) 
              
              // LOG_PDF_NORM(par.sigma_X, 10) + LOG_PDF_NORM(par.sigma_Y, 10) + log(4.0));
}

// [[Rcpp::export]]
double scaled_nl_posterior_MR(unsigned J, unsigned N, arma::mat SS, arma::vec sigma_G, List param_list, List prior_sd, unsigned n = 2) {
  
  arma::mat f(J + 3, 2, arma::fill::zeros);
  arma::mat sigma(2, 2, arma::fill::zeros);
  // arma::vec sigma_G(J, arma::fill::zeros);
  
  struct Parameters par = read_parameter_list(param_list);
  struct Parameters sd = read_prior_list(prior_sd);
  
  for (unsigned j = 0; j < J; ++j) {
    // sigma_G(j) = sqrt((1 - SS(0, j+1) / n) * SS(0, j+1)); 
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
  
  if (sign < 0) std::cerr << "WARNING: negative determinant of sigma in scaled_nl_posterior_MR." << std::endl;
  
  double neg_log_prior = neg_log_prior_MR(param_list, prior_sd);
  
  return 0.5 * (log(4 * M_PI * M_PI) + val + trace(inv(sigma) * S)) + neg_log_prior / N;
}


// [[Rcpp::export]]
arma::vec scaled_nl_gradient_MR(unsigned J, unsigned N, arma::mat SS, arma::vec sigma_G, List param_list, List prior_sd, unsigned n = 2) {
  
  arma::vec gradient(2 * J + 5, arma::fill::zeros);
  arma::mat sf(J + 3, 2, arma::fill::zeros);
  arma::mat iscov(2, 2, arma::fill::zeros);
  // arma::vec sigma_G(J, arma::fill::zeros);
  
  struct Parameters par = read_parameter_list(param_list);
  struct Parameters sd = read_prior_list(prior_sd);
  
  
  for (unsigned j = 0; j < J; ++j) {
    // sigma_G(j) = sqrt((1 - SS(0, j+1) / n) * SS(0, j+1)); 
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
  
  // std::cout << "ihalf_Sigma_Z: " << std::endl;
  // std::cout << ihalf_Sigma_Z << std::endl << std::endl;
  // std::cout << "sf: " << std::endl;
  // std::cout << sf << std::endl << std::endl;
  // std::cout << "iscov: " << std::endl;
  // std::cout << iscov << std::endl << std::endl;
  // std::cout << common_term << std::endl;
  
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


// [[Rcpp::export]]
arma::mat scaled_nl_hessian_MR(unsigned J, unsigned N, arma::mat SS, arma::vec sigma_G, List param_list, List prior_sd, unsigned n = 2) {
  
  struct Parameters par = read_parameter_list(param_list);
  struct Parameters sd = read_prior_list(prior_sd);
  
  arma::vec gradient(2 * J + 5, arma::fill::zeros);
  arma::mat hessian(2 * J + 5, 2 * J + 5, arma::fill::zeros);
  
  arma::mat sf(J + 3, 2, arma::fill::zeros);
  arma::mat iscov(2, 2, arma::fill::zeros);
  // arma::vec sigma_G(J, arma::fill::zeros);

  for (unsigned j = 0; j < J; ++j) {
    // sigma_G(j) = sqrt((1 - SS(0, j+1) / n) * SS(0, j+1)); 
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
      // if (j == 0 && k == 0) {
      //   std::cout << trace(iscov * d_sgamma_sf.slice(j).t() * sSS * d_sgamma_sf.slice(k)) << std::endl;
      //   std::cout << (j == k) / (sd.sgamma(j) * sd.sgamma(j)) << std::endl;
      // }
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
                  
  // std::cout << "ihalf_Sigma_Z: " << std::endl;
  // std::cout << ihalf_Sigma_Z << std::endl << std::endl;
  // std::cout << "sf: " << std::endl;
  // std::cout << sf << std::endl << std::endl;
  // std::cout << "iscov: " << std::endl;
  // std::cout << iscov << std::endl << std::endl;
  // std::cout << common_term << std::endl;
  
  // // derivatives w.r.t sgamma and salpha
  // for (unsigned j = 0; j < J; ++j) {
  //   gradient(j) = - (common_term(j+1, 0) + par.sbeta * common_term(j+1, 1)) + par.sgamma(j) / (sd.sgamma(j) * sd.sgamma(j) * N);
  //   gradient(j+J) = - common_term(j+1, 1) + par.salpha(j) / (sd.salpha(j) * sd.salpha(j) * N);
  // }
  // 
  // // derivative w.r.t sbeta
  // arma::mat d_sbeta_scov(2, 2, arma::fill::zeros);
  // d_sbeta_scov(0, 1) = d_sbeta_scov(1, 0) = 1 + par.skappa_X * par.skappa_X;
  // d_sbeta_scov(1, 1) = 2 * par.sbeta * (1 + par.skappa_X * par.skappa_X) + 2 * par.skappa_X * par.skappa_Y;
  // gradient(2*J) = - 0.5 * trace(d_sbeta_scov * iscov * (sf.t() * common_term - I_2)) + par.sbeta / (sd.sbeta * sd.sbeta * N);
  // for (unsigned j = 0; j < J; ++j) gradient(2*J) -= common_term(j+1, 1) * par.sgamma(j);
  // 
  // arma::mat d_skappa_X_scov(2, 2, arma::fill::zeros);
  // d_skappa_X_scov(0, 0) = 2 * par.skappa_X;
  // d_skappa_X_scov(0, 1) = d_skappa_X_scov(1, 0) = 2 * par.sbeta * par.skappa_X + par.skappa_Y;
  // d_skappa_X_scov(1, 1) = 2 * par.sbeta * (par.skappa_Y + par.sbeta * par.skappa_X);
  // gradient(2*J+1) = - 0.5 * trace(d_skappa_X_scov * iscov * (sf.t() * common_term - I_2)) + par.skappa_X / (sd.skappa_X * sd.skappa_X * N);
  // 
  // arma::mat d_skappa_Y_scov(2, 2, arma::fill::zeros);
  // d_skappa_Y_scov(0, 1) = d_skappa_Y_scov(1, 0) = par.skappa_X;
  // d_skappa_Y_scov(1, 1) = 2 * (par.skappa_Y + par.sbeta * par.skappa_X);
  // gradient(2*J+2) = - 0.5 * trace(d_skappa_Y_scov * iscov * (sf.t() * common_term - I_2)) + par.skappa_Y / (sd.skappa_Y * sd.skappa_Y * N);
  // 
  // // Huge simplification can be obtained here
  // gradient(2*J+3) = - (common_term(J+1, 0) - 1) / par.sigma_X; // + par.sigma_X / (100 * N);
  // gradient(2*J+4) = - (common_term(J+2, 1) - 1) / par.sigma_Y; // + par.sigma_Y / (100 * N);
  
  return hessian;
}
