#' Function for converting vector of parameters to list format.
#'
#' @param param_vector Vector of MASSIVE model parameters.
#'
#' @return List containing same parameters.
#' @export
parameter_vector_to_list <- function(param_vector) {
  J <- (length(param_vector) - 5) / 2
  list(
    sgamma = param_vector[1:J],
    salpha = param_vector[(J+1):(2*J)],
    sbeta = param_vector[2*J+1],
    skappa_X = param_vector[2*J+2],
    skappa_Y = param_vector[2*J+3],
    sigma_X = param_vector[2*J+4],
    sigma_Y = param_vector[2*J+5]
  )
}

#' Function for converting list of parameters to vector format, often necessary 
#' when running optimization routines.
#'
#' @param param_list List of MASSIVE model parameters. 
#'
#' @return Numeric vector containing same parameters.
#' @export
#'
#' @examples
#' par <- random_Gaussian_parameters(10)
#' parameter_list_to_vector(par)
parameter_list_to_vector <- function(param_list) {
  c(
    param_list$sgamma,
    param_list$salpha,
    param_list$sbeta,
    param_list$skappa_X,
    param_list$skappa_Y,
    param_list$sigma_X,
    param_list$sigma_Y
  )
}


#' Generate random Gaussian distributed parameters
#'
#' @param J Integer number of instrumental variables.
#' @param log_scale Logical flag indicating whether scale parameters sigma_X
#' and sigma_Y should be computed on the log-scale (normally distributed), or
#' on the normal scale (log-normally distributed).
#'
#' @return List of random parameters for the IV model generated from Gaussian distributions.
#' @export
#'
#' @examples 
#' random_Gaussian_parameters(10)
#' random_Gaussian_parameters(10, log_scale = TRUE)
random_Gaussian_parameters <- function(J, log_scale = FALSE) {
  sgamma <- stats::rnorm(J)
  salpha <- stats::rnorm(J)
  sbeta <- stats::rnorm(1)
  skappa_X <- stats::rnorm(1)
  skappa_Y <- stats::rnorm(1)
  
  if (log_scale) {
    sigma_X <- stats::rnorm(1)
    sigma_Y <- stats::rnorm(1)
  } else {
    sigma_X <- exp(stats::rnorm(1))
    sigma_Y <- exp(stats::rnorm(1))
  }
  
  list(
    sgamma = sgamma, 
    salpha = salpha, 
    sbeta = sbeta,
    skappa_X = skappa_X,
    skappa_Y = skappa_Y,
    sigma_X = sigma_X,
    sigma_Y = sigma_Y
  )
}


#' Function for deriving maximum likelihood solution given values for scale-free
#' confounding coefficients.
#' @param SS Numeric matrix containing first- and second-order statistics.
#' @param skappa_X Numeric scale-free confounding coefficient on exposure.
#' @param skappa_Y Numeric scale-free confounding coefficient on outcome.
#' @param sigma_G Numeric vector of instrument standard deviations.
#'
#' @return List of parameters on the ML manifold with the same values for
#' the scaled confounding coefficients (skappa_X, skappa_Y).
#' @export
#'
#' @examples
#' J <- 5 # number of instruments
#' par <- random_Gaussian_parameters(J)
#' generate_data_MASSIVE_model(N = 1000, n = 2, p = rep(0.3, J), par)
get_ML_solution <- function(SS, skappa_X, skappa_Y, sigma_G) {
  
  J <- nrow(SS) - 3
  
  cov_XY_G <- SS[(J+2):(J+3), (J+2):(J+3)] - SS[(J+2):(J+3), 2:(J+1)] %*% solve(SS[2:(J+1), 2:(J+1)]) %*% SS[2:(J+1), (J+2):(J+3)]
  
  ML_sigma_X <- sqrt(cov_XY_G[1, 1] / (1 + skappa_X^2))
  ML_sigma_Y <- sqrt((cov_XY_G[1, 1] * cov_XY_G[2, 2] - cov_XY_G[1, 2]^2) * (1 + skappa_X^2) / ((1 + skappa_X^2 + skappa_Y^2) * cov_XY_G[1, 1]))
  ML_sbeta <- (cov_XY_G[1, 2] / (ML_sigma_X * ML_sigma_Y) - skappa_X * skappa_Y) / (1 + skappa_X^2)
  
  ML_sgamma <- solve(SS[2:(J+1), 2:(J+1)], SS[J+2, 2:(J+1)]) / ML_sigma_X * sigma_G
  ML_sGamma <- solve(SS[2:(J+1), 2:(J+1)], SS[J+3, 2:(J+1)]) / ML_sigma_Y * sigma_G
  
  list(
    sgamma = ML_sgamma,
    salpha = ML_sGamma - ML_sbeta * ML_sgamma,
    sbeta = ML_sbeta,
    skappa_X = skappa_X,
    skappa_Y = skappa_Y,
    sigma_X = ML_sigma_X,
    sigma_Y = ML_sigma_Y
  )
}