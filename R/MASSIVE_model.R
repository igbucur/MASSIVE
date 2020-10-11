#' Routine for computing the posterior distribution of the MASSIVE model
#'
#' @param J Integer number of genetic instrumental variables.
#' @param N Integer number of observations.
#' @param SS Numeric matrix containing first- and second-order statistics.
#' @param sigma_G Numeric vector of genetic IV standard deviations.
#' @param param_list List of IV model parameter values.
#' @param prior_sd List of standard deviations for the parameter Gaussian priors.
#' @param n Integer number of alleles (trials) for the binomial genetic variants.
#' 
#' @return numeric value; posterior value for the MASSIVE model.
#' @export
scaled_nl_posterior_log <- function(J, N, SS, sigma_G, param_list, prior_sd, n = 2) {
  param_list$sigma_X <- exp(param_list$sigma_X)
  param_list$sigma_Y <- exp(param_list$sigma_Y)
  post <- scaled_nl_posterior_MR(J, N, SS, sigma_G, param_list, prior_sd, n)
  # print(paste("sbeta", param_list$sbeta, "skappa_X", param_list$skappa_X, "skappa_Y", param_list$skappa_Y, "post", post))
  post
}

#' Routine for computing the posterior gradient for the MASSIVE model
#'
#' @param J Integer number of genetic instrumental variables.
#' @param N Integer number of observations.
#' @param SS Numeric matrix containing first- and second-order statistics.
#' @param sigma_G Numeric vector of genetic IV standard deviations.
#' @param param_list List of IV model parameter values.
#' @param prior_sd List of standard deviations for the parameter Gaussian priors.
#' @param n Integer number of alleles (trials) for the binomial genetic variants.
#'
#' @return numeric vector; posterior gradient for the MASSIVE model.
#' @export
scaled_nl_gradient_log <- function(J, N, SS, sigma_G, param_list, prior_sd, n = 2) {
  param_list$sigma_X <- exp(param_list$sigma_X)
  param_list$sigma_Y <- exp(param_list$sigma_Y)
  grad <- scaled_nl_gradient_MR(J, N, SS, sigma_G, param_list, prior_sd, n)
  grad[2*J+4] <- grad[2*J+4] * param_list$sigma_X
  grad[2*J+5] <- grad[2*J+5] * param_list$sigma_Y
  grad
}

#' Routine for computing the posterior Hessian for the MASSIVE model
#'
#' @param J Integer number of genetic instrumental variables.
#' @param N Integer number of observations.
#' @param SS Numeric matrix containing first- and second-order statistics.
#' @param sigma_G Numeric vector of genetic IV standard deviations.
#' @param param_list List of IV model parameter values.
#' @param prior_sd List of standard deviations for the parameter Gaussian priors.
#' @param n Integer number of alleles (trials) for the binomial genetic variants.
#'
#' @return numeric matrix; posterior Hessian for the MASSIVE model.
#' @export
scaled_nl_hessian_log <- function(J, N, SS, sigma_G, param_list, prior_sd, n = 2) {
  param_list$sigma_X <- exp(param_list$sigma_X)
  param_list$sigma_Y <- exp(param_list$sigma_Y)
  
  grad <- scaled_nl_gradient_MR(J, N, SS, sigma_G, param_list, prior_sd, n)
  hess <- scaled_nl_hessian_MR(J, N, SS, sigma_G, param_list, prior_sd, n)
  
  hess[2*J+4, ] <- hess[2*J+4, ] * param_list$sigma_X
  hess[2*J+5, ] <- hess[2*J+5, ] * param_list$sigma_Y
  hess[, 2*J+4] <- hess[, 2*J+4] * param_list$sigma_X
  hess[, 2*J+5] <- hess[, 2*J+5] * param_list$sigma_Y
  
  hess[2*J+4, 2*J+4] <- hess[2*J+4, 2*J+4] + grad[2*J+4] * param_list$sigma_X
  hess[2*J+5, 2*J+5] <- hess[2*J+5, 2*J+5] + grad[2*J+5] * param_list$sigma_Y
  hess
}
