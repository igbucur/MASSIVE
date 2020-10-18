
#' Routine for more robust local optimum search of MASSIVE posterior, where
#' the search is repeated until a better value cannot be found any more.
#'
#' @param J Integer number of genetic instrumental variables.
#' @param N Integer number of observations.
#' @param SS Numeric matrix containing first- and second-order statistics.
#' @param sigma_G Numeric vector of genetic IV standard deviations.
#' @param prior_sd List of standard deviations for the parameter Gaussian priors.
#' @param skappa_X Scale-free confounding coefficient to exposure used for initialization.
#' @param skappa_Y Scale-free confounding coefficient to outcome used for initialization.
#' @param tol Numeric tolerance value used to decide if a better optimum was found.
#' @param post_fun Function for computing the IV model posterior value.
#' @param gr_fun Function for computing the IV model posterior gradient.
#' @param hess_fun Function for computing the IV model posterior Hessian.
#'
#' @return optim object (see \link[stats]{optim} containing optimum parameters, 
#' the value obtained at the optimum, and the Hessian at the optimum.
#' @export
#'
#' @examples
#' J <- 5 # number of instruments
#' N <- 1000 # number of samples
#' parameters <- random_Gaussian_parameters(J) 
#' EAF <- runif(J, 0.1, 0.9) # EAF random values
#' dat <- generate_data_MASSIVE_model(N, 2, EAF, parameters)
#' robust_find_optimum(
#'   J, N, dat$SS, binomial_sigma_G(dat$SS), 
#'   prior_sd = decode_IV_model(get_random_IV_model(J), 1, 0.01), 
#'   skappa_X = 1, skappa_Y = 1
#' )
robust_find_optimum <- function(J, N, SS, sigma_G, prior_sd, skappa_X, skappa_Y, tol = 1e-6,
                         post_fun = scaled_neg_log_posterior, gr_fun = scaled_neg_log_gradient, hess_fun = scaled_neg_log_hessian) {
  

  
  crt_opt <- parameter_vector_to_list(rep(0, 2*J+5))
  new_opt <- TRUE
  
  while (new_opt == TRUE) {
    
    par <- get_ML_solution(SS, skappa_X, skappa_Y, sigma_G)
    par$sigma_X <- log(par$sigma_X)
    par$sigma_Y <- log(par$sigma_Y)
    
    optim_MAP <- stats::optim(parameter_list_to_vector(par), fn = function(x) {
      params <- list(
        sgamma = x[1:J],
        salpha = x[(J+1):(2*J)],
        sbeta = x[2*J+1],
        skappa_X = x[2*J+2],
        skappa_Y = x[2*J+3],
        sigma_X = x[2*J+4],
        sigma_Y = x[2*J+5]
      )
      post_fun(J, N, SS, sigma_G, params, prior_sd)
    }, gr = function(x) {
      params <- list(
        sgamma = x[1:J],
        salpha = x[(J+1):(2*J)],
        sbeta = x[2*J+1],
        skappa_X = x[2*J+2],
        skappa_Y = x[2*J+3],
        sigma_X = x[2*J+4],
        sigma_Y = x[2*J+5]
      )
      gr_fun(J, N, SS, sigma_G, params, prior_sd)
    }, method = "L-BFGS-B", control = list(maxit = 50000, factr = 1))
  
    optim_par <- parameter_vector_to_list(optim_MAP$par)
    skappa_X <- optim_par$skappa_X
    skappa_Y <- optim_par$skappa_Y
    
    # only accept new optimum if it is significantly better
    new_opt <- (post_fun(J, N, SS, sigma_G, optim_par, prior_sd) - 
      post_fun(J, N, SS, sigma_G, crt_opt, prior_sd)) < - tol
    crt_opt <- optim_par
  }
  
  optim_MAP$hessian <- hess_fun(J, N, SS, sigma_G, parameter_vector_to_list(optim_MAP$par), prior_sd)
  
  optim_MAP
}


#' Routine for finding MASSIVE posterior local optimum.
#'
#' @param J Integer number of genetic instrumental variables.
#' @param N Integer number of observations.
#' @param SS Numeric matrix containing first- and second-order statistics.
#' @param sigma_G Numeric vector of genetic IV standard deviations.
#' @param prior_sd List of standard deviations for the parameter Gaussian priors.
#' @param skappa_X Scale-free confounding coefficient to exposure used for initialization.
#' @param skappa_Y Scale-free confounding coefficient to outcome used for initialization.
#' @param post_fun Function for computing the IV model posterior value.
#' @param gr_fun Function for computing the IV model posterior gradient.
#' @param hess_fun Function for computing the IV model posterior Hessian.
#'
#' @return optim object (see \link[stats]{optim} containing optimum parameters, 
#' the value obtained at the optimum, and the Hessian at the optimum.
#' @export
#'
#' @examples
#' J <- 5 # number of instruments
#' N <- 1000 # number of samples
#' parameters <- random_Gaussian_parameters(J) 
#' EAF <- runif(J, 0.1, 0.9) # EAF random values
#' dat <- generate_data_MASSIVE_model(N, 2, EAF, parameters)
#' find_optimum(
#'   J, N, dat$SS, binomial_sigma_G(dat$SS), 
#'   prior_sd = decode_IV_model(get_random_IV_model(J), 1, 0.01), 
#'   skappa_X = 1, skappa_Y = 1
#' )
find_optimum <- function(J, N, SS, sigma_G, prior_sd, skappa_X, skappa_Y, 
                         post_fun = scaled_neg_log_posterior, gr_fun = scaled_neg_log_gradient, hess_fun = scaled_neg_log_hessian) {
  
  par <- get_ML_solution(SS, skappa_X, skappa_Y, sigma_G)
  optim_MAP <- stats::optim(c(
    par$sgamma, par$salpha, par$sbeta, par$skappa_X, par$skappa_Y, log(par$sigma_X), log(par$sigma_Y)
  ), fn = function(x) {
    params <- list(
      sgamma = x[1:J],
      salpha = x[(J+1):(2*J)],
      sbeta = x[2*J+1],
      skappa_X = x[2*J+2],
      skappa_Y = x[2*J+3],
      sigma_X = x[2*J+4],
      sigma_Y = x[2*J+5]
    )
    post_fun(J, N, SS, sigma_G, params, prior_sd)
  }, gr = function(x) {
    params <- list(
      sgamma = x[1:J],
      salpha = x[(J+1):(2*J)],
      sbeta = x[2*J+1],
      skappa_X = x[2*J+2],
      skappa_Y = x[2*J+3],
      sigma_X = x[2*J+4],
      sigma_Y = x[2*J+5]
    )
    gr_fun(J, N, SS, sigma_G, params, prior_sd)
  }, method = "L-BFGS-B", control = list(maxit = 50000, factr = 1))
  
  optim_MAP$hessian <- hess_fun(J, N, SS, sigma_G, parameter_vector_to_list(optim_MAP$par), prior_sd)
  
  optim_MAP
}

