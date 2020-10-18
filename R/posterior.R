#' Routine for computing the negative log-posterior distribution of the MASSIVE model
#'
#' @param J Integer number of genetic instrumental variables.
#' @param N Integer number of observations.
#' @param SS Numeric matrix containing first- and second-order statistics.
#' @param sigma_G Numeric vector of genetic IV standard deviations.
#' @param param_list List of IV model parameter values.
#' @param prior_sd List of standard deviations for the parameter Gaussian priors.
#' @param n Integer number of alleles (trials) for the binomial genetic variants.
#' @param log_scale Logical flag indicating whether scale parameters (sigma_X, sigma_Y)
#' are given on the log-scale or on the normal scale.
#' 
#' @return numeric value; negative log-posterior value for the MASSIVE model.
#' @export
scaled_neg_log_posterior <- function(J, N, SS, sigma_G, param_list, prior_sd, n = 2, log_scale = TRUE) {
  
  if (log_scale == TRUE) {
    param_list$sigma_X <- exp(param_list$sigma_X)
    param_list$sigma_Y <- exp(param_list$sigma_Y)
  }
  
  post <- Rcpp_scaled_neg_log_posterior(J, N, SS, sigma_G, param_list, prior_sd, n)
  post
}

#' Routine for computing the negative log-posterior gradient for the MASSIVE model
#'
#' @param J Integer number of genetic instrumental variables.
#' @param N Integer number of observations.
#' @param SS Numeric matrix containing first- and second-order statistics.
#' @param sigma_G Numeric vector of genetic IV standard deviations.
#' @param param_list List of IV model parameter values.
#' @param prior_sd List of standard deviations for the parameter Gaussian priors.
#' @param n Integer number of alleles (trials) for the binomial genetic variants.
#' @param log_scale Logical flag indicating whether scale parameters (sigma_X, sigma_Y)
#' are given on the log-scale or on the normal scale.
#'
#' @return numeric vector; negative log-posterior gradient for the MASSIVE model.
#' @export
scaled_neg_log_gradient <- function(J, N, SS, sigma_G, param_list, prior_sd, n = 2, log_scale = TRUE) {
  
  if (log_scale == TRUE) {
    param_list$sigma_X <- exp(param_list$sigma_X)
    param_list$sigma_Y <- exp(param_list$sigma_Y)
  }

  grad <- Rcpp_scaled_neg_log_gradient(J, N, SS, sigma_G, param_list, prior_sd, n)
  
  if (log_scale == TRUE) {
    grad[2*J+4] <- grad[2*J+4] * param_list$sigma_X
    grad[2*J+5] <- grad[2*J+5] * param_list$sigma_Y
  }
  
  grad
}

#' Routine for computing the negative log-posterior Hessian for the MASSIVE model
#'
#' @param J Integer number of genetic instrumental variables.
#' @param N Integer number of observations.
#' @param SS Numeric matrix containing first- and second-order statistics.
#' @param sigma_G Numeric vector of genetic IV standard deviations.
#' @param param_list List of IV model parameter values.
#' @param prior_sd List of standard deviations for the parameter Gaussian priors.
#' @param n Integer number of alleles (trials) for the binomial genetic variants.
#' @param log_scale Logical flag indicating whether scale parameters (sigma_X, sigma_Y)
#' are given on the log-scale or on the normal scale.
#'
#' @return numeric matrix; negative log-posterior Hessian for the MASSIVE model.
#' @export
scaled_neg_log_hessian <- function(J, N, SS, sigma_G, param_list, prior_sd, n = 2, log_scale = TRUE) {
  
  if (log_scale == TRUE) {
    param_list$sigma_X <- exp(param_list$sigma_X)
    param_list$sigma_Y <- exp(param_list$sigma_Y)
  }
  
  grad <- Rcpp_scaled_neg_log_gradient(J, N, SS, sigma_G, param_list, prior_sd, n)
  hess <- Rcpp_scaled_neg_log_hessian(J, N, SS, sigma_G, param_list, prior_sd, n)
  
  if (log_scale == TRUE) {
    hess[2*J+4, ] <- hess[2*J+4, ] * param_list$sigma_X
    hess[2*J+5, ] <- hess[2*J+5, ] * param_list$sigma_Y
    hess[, 2*J+4] <- hess[, 2*J+4] * param_list$sigma_X
    hess[, 2*J+5] <- hess[, 2*J+5] * param_list$sigma_Y
    
    hess[2*J+4, 2*J+4] <- hess[2*J+4, 2*J+4] + grad[2*J+4] * param_list$sigma_X
    hess[2*J+5, 2*J+5] <- hess[2*J+5, 2*J+5] + grad[2*J+5] * param_list$sigma_Y
  }
  
  hess
}
