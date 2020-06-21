#' Routines for computing the posterior distribution of the MASSIVE model

scaled_nl_posterior_log <- function(J, N, SS, sigma_G, param_list, prior_sd, n = 2) {
  param_list$sigma_X <- exp(param_list$sigma_X)
  param_list$sigma_Y <- exp(param_list$sigma_Y)
  post <- scaled_nl_posterior_MR(J, N, SS, sigma_G, param_list, prior_sd, n)
  post
}

scaled_nl_gradient_log <- function(J, N, SS, sigma_G, param_list, prior_sd, n = 2) {
  param_list$sigma_X <- exp(param_list$sigma_X)
  param_list$sigma_Y <- exp(param_list$sigma_Y)
  grad <- scaled_nl_gradient_MR(J, N, SS, sigma_G, param_list, prior_sd, n)
  grad[2*J+4] <- grad[2*J+4] * param_list$sigma_X
  grad[2*J+5] <- grad[2*J+5] * param_list$sigma_Y
  grad
}

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
