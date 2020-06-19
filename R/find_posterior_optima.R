
#' Routine for (more robustly) finding a MASSIVE posterior local optimum
#'
#' @param J 
#' @param N 
#' @param SS 
#' @param sigma_G 
#' @param prior_sd 
#' @param skappa_X 
#' @param skappa_Y 
#' @param tol 
#' @param post_fun 
#' @param gr_fun 
#' @param hess_fun 
#'
#' @return
#' @export
#'
#' @examples
robust_find_optimum <- function(J, N, SS, sigma_G, prior_sd, skappa_X, skappa_Y, tol = 1e-6,
                         post_fun = scaled_nl_posterior_log, gr_fun = scaled_nl_gradient_log, hess_fun = scaled_nl_hessian_log) {
  

  
  crt_opt <- parameter_vector_to_list(rep(0, 2*J+5))
  new_opt <- TRUE
  
  while (new_opt == TRUE) {
    
    par <- get_ML_solution(SS, skappa_X, skappa_Y, sigma_G)
    par$sigma_X <- log(par$sigma_X)
    par$sigma_Y <- log(par$sigma_Y)
    
    optim_MAP <- optim(parameter_list_to_vector(par), fn = function(x) {
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


find_optimum <- function(J, N, SS, sigma_G, prior_sd, skappa_X, skappa_Y, 
                         post_fun = scaled_nl_posterior_log, gr_fun = scaled_nl_gradient_log, hess_fun = scaled_nl_hessian_log) {
  
  par <- get_ML_solution(SS, skappa_X, skappa_Y, sigma_G)
  optim_MAP <- optim(c(
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
