#' Function to extract first-order and second-order estimated moments from
#' summary data.
#'
#' @param J Integer number of instrumental variables.
#' @param beta_XG Numeric vector of estimated G -> X regression terms.
#' @param sigma_XG Numeric vector of estimated G -> X standard errors.
#' @param nobs_XG Integer sample size in G -> X regression (GWA) study.
#' @param beta_YG Numeric vector of estimated G -> Y regression terms.
#' @param sigma_YG Numeric vector of estimated G -> Y standard errors.
#' @param nobs_YG Integer sample size in G -> Y regression (GWA) study.
#' @param beta_YX Numeric vector of estimated X -> Y regression value.
#' @param EAF Numeric vector of effect allele frequencies
#' @param n Integer number of alleles (trials) for binomial genetic variant.
#'
#' @return Numeric matrix of size (J+3)x(J+3) containing the first-order and 
#' second-order moments of the (J+2) vector (G, X, Y).
#' @export
#'
#' @examples
#' MR_regression_coefficients_to_moments(1, 1, 1e-3, 1000, 1, 1e-3, 1000, 1, 1e-3)
MR_regression_coefficients_to_moments <- function(J, beta_XG, sigma_XG, nobs_XG, beta_YG, sigma_YG, nobs_YG, beta_YX, EAF = rep(0.5, J), n = 2) {
  
  SSS <- matrix(0, J + 3, J + 3) # EV{[1 G X Y] [1 G X Y]^T}
  
  SSS[1, 1] <- 1
  SSS[1, 2:(J+1)] <- SSS[2:(J+1), 1] <- n * EAF
  SSS[2:(J+1), 2:(J+1)] <- n * n * EAF %*% t(EAF) + diag(n * EAF * (1 - EAF), J) # E{GG^T}
  
  SSS[1, J+2] <- SSS[J+2, 1] <- t(beta_XG) %*% SSS[1, 2:(J+1)]
  SSS[1, J+3] <- SSS[J+3, 1] <- t(beta_YG) %*% SSS[1, 2:(J+1)]
  
  SSS[J+2, 2:(J+1)] <- SSS[2:(J+1), J+2] <- SSS[2:(J+1), 2:(J+1)] %*% beta_XG
  SSS[J+3, 2:(J+1)] <- SSS[2:(J+1), J+3] <- SSS[2:(J+1), 2:(J+1)] %*% beta_YG
  
  est_var_X <- stats::median((beta_XG^2 + sigma_XG^2 * nobs_XG) * n * EAF * (1 - EAF))
  est_var_Y <- stats::median((beta_YG^2 + sigma_YG^2 * nobs_YG) * n * EAF * (1 - EAF))
  
  SSS[J+2, J+2] <- est_var_X + (n * t(beta_XG) %*% EAF)^2
  SSS[J+2, J+3] <- SSS[J+3, J+2] <- est_var_X * beta_YX + SSS[1, J+2] * SSS[J+3, 1]
  SSS[J+3, J+3] <- est_var_Y + (n * t(beta_YG) %*% EAF)^2
  
  SSS
}


#' Function to generate data from a MASSIVE generating model.
#'
#' @param N Integer number of observations.
#' @param n Integer number of alleles (trials) for the binomial genetic variables.
#' @param p Numeric expected allele frequencies for the genetic variables.
#' @param par List of MASSIVE model parameters
#' @param seed Integer representing seed for random number generation.
#' @param log_scale Logical flag indicating whether scale parameters sigma_X
#' and sigma_Y are given on the log-scale or normal scale.
#'
#' @return List containing generate data vector as well as the scatter matrix
#' of first-order and second-order statistics.
#' @export
#'
#' @examples
#' J <- 5 # number of instruments
#' par <- random_Gaussian_parameters(J)
#' generate_data_MASSIVE_model(N = 1000, n = 2, p = rep(0.3, J), par)
generate_data_MASSIVE_model <- function(N, n, p, par, seed = NULL, log_scale = FALSE) {
  
  set.seed(seed)
  
  if (log_scale) {
    par$sigma_X <- exp(par$sigma_X)
    par$sigma_Y <- exp(par$sigma_Y)
  }
  
  sigma_G <- sqrt(n * p * (1-p))
  gamma <- par$sgamma * par$sigma_X / sigma_G
  alpha <- par$salpha * par$sigma_Y / sigma_G
  beta <- par$sbeta * par$sigma_Y / par$sigma_X
  kappa_X <- par$skappa_X * par$sigma_X
  kappa_Y <- par$skappa_Y * par$sigma_Y
  
  G <- sapply(p, function(t) stats::rbinom(N, n, t)) # genetic variant
  U <- stats::rnorm(N) # confounder
  X <- G %*% gamma + kappa_X * U + stats::rnorm(N, sd = par$sigma_X)
  Y <- G %*% alpha + kappa_Y * U + beta * X + stats::rnorm(N, sd = par$sigma_Y)
  Z <- cbind(1, G, X, Y)
  SS <- t(Z) %*% Z / N # first and second-order moments
  
  list(Z = Z, SS = SS)
}



#' Functioning for automatically determining MASSIVE hyperparameters from data,
#' where we assume a N(0, 10) prior on skappa_X.
#'
#' @param J Integer number of candidate instruments.
#' @param N Integer number of observations.
#' @param SS Numeric matrix containing first- and second-order statistics.
#' @param sigma_G Numeric vector of instrument standard deviations.
#'
#' @return List containing sd_slab and sd_spike hyperparameters, indicating
#' the standard deviation of the slab and spike component, respectively.
#' @export
#'
#' @examples
#' J <- 5 # number of instruments
#' N <- 1000 # number of samples
#' par <- random_Gaussian_parameters(J)
#' dat <- generate_data_MASSIVE_model(N, n = 2, p = rep(0.3, J), par)
#' determine_hyperparameters(J, N, dat$SS, binomial_sigma_G(dat$SS))
determine_hyperparameters <- function(J, N, SS, sigma_G) {
  
  ML <- get_ML_solution(SS, 0, 0, sigma_G); 
  
  var_kappa_X <- 1 + 100 # 1 + sd_skappa_X^2, derived from N(0, 10) prior on skappa_X

  est_sd_slab <- sqrt(sum(ML$sgamma^2) * var_kappa_X / J)
  
  c_root <- stats::uniroot(function(C) {
    (N + 1 - C) * (min(ML$sgamma^2) * var_kappa_X) / est_sd_slab^2 + log(C)
  }, interval = c(N+1, 1e12))
  
  est_sd_spike <- est_sd_slab / sqrt(c_root$root)
  
  list(sd_slab = est_sd_slab, sd_spike = est_sd_spike)
}


#' Functioning for automatically determining MASSIVE hyperparameters from data.
#' where we assume a N(0, sd_slab) prior on skappa_X.
#'
#' @param J Integer number of candidate instruments.
#' @param N Integer number of observations.
#' @param SS Numeric matrix containing first- and second-order statistics.
#' @param sigma_G Numeric vector of instrument standard deviations.
#'
#' @return List containing sd_slab and sd_spike hyperparameters, indicating
#' the standard deviation of the slab and spike component, respectively.
#' @export
#'
#' @examples
#' J <- 5 # number of instruments
#' N <- 1000 # number of samples
#' par <- random_Gaussian_parameters(J)
#' dat <- generate_data_MASSIVE_model(N, n = 2, p = rep(0.3, J), par)
#' alternative_determine_hyperparameters(J, N, dat$SS, binomial_sigma_G(dat$SS))
alternative_determine_hyperparameters <- function(J, N, SS, sigma_G) {
  
  ML <- get_ML_solution(SS, 0, 0, sigma_G); 
  est_sd_slab <- sqrt(sum(ML$sgamma^2) / J)

  c_root <- stats::uniroot(function(C) {
    (N + 1 - C) * min(ML$sgamma^2) * (1 + est_sd_slab^2) / est_sd_slab^2 + log(C)
  }, interval = c(N+1, 1e12))
  
  est_sd_spike <- est_sd_slab / sqrt(c_root$root)
  
  list(sd_slab = est_sd_slab, sd_spike = est_sd_spike)
}


#' Function for deriving the estimated standard deviation of the genetic variants 
#' from the first-order and second-order statistics, assuming these are 
#' binomially distributed.
#'
#' @param SS Numeric matrix of first-order and second-order statistics.
#' @param n Integer number of alleles (trials) for the binomial genetic variables. 
#'
#' @return Numeric vector containing the standard deviations of the genetic variants.
#' @export
#'
#' @examples
#' J <- 5 # number of instruments
#' par <- random_Gaussian_parameters(J)
#' dat <- generate_data_MASSIVE_model(N = 1000, n = 2, p = rep(0.3, J), par)
#' binomial_sigma_G(dat$SS)
binomial_sigma_G <- function(SS, n = 2) {
  J <- nrow(SS) - 3
  sqrt(SS[1, 2:(J+1)] * (1 - SS[1, 2:(J+1)] / n))
}


