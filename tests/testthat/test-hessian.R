test_that("Rcpp_scaled_neg_log_hessian is the Hessian of Rcpp_scaled_neg_log_posterior", {
  
  J <- 20
  N <- 1000
  
  true_par <- random_Gaussian_parameters(J)
  sd_spike <- runif(1, 1e-3, 1e-1)
  sd_slab <- runif(1, 1e-1, 1e1)
  
  data <- generate_data_MASSIVE_model(N, 2, rep(0.3, J), true_par)
  
  SS <- data$SS
  sigma_G <- binomial_sigma_G(SS, 2)
  
  rand_par <- random_Gaussian_parameters(J)
  
  model <- get_random_IV_model(J)
  prior <- decode_IV_model(model, sd_slab, sd_spike)

  hess <- Rcpp_scaled_neg_log_hessian(J, N, SS, sigma_G, rand_par, prior)
  hess_benchmark <- numDeriv::hessian(function(x) {
    Rcpp_scaled_neg_log_posterior(J, N, SS, sigma_G, parameter_vector_to_list(x), prior)
  }, parameter_list_to_vector(rand_par))
  
  expect_equal(hess, hess_benchmark, tolerance = 1e-6)
})

test_that("Rcpp_scaled_neg_log_hessian is the Jacobian of Rcpp_scaled_neg_log_gradient", {
  
  J <- 20
  N <- 1000
  
  true_par <- random_Gaussian_parameters(J)
  sd_spike <- runif(1, 1e-3, 1e-1)
  sd_slab <- runif(1, 1e-1, 1e1)
  
  data <- generate_data_MASSIVE_model(N, 2, rep(0.3, J), true_par)
  
  SS <- data$SS
  sigma_G <- binomial_sigma_G(SS, 2)
  
  rand_par <- random_Gaussian_parameters(J)
  
  model <- get_random_IV_model(J)
  prior <- decode_IV_model(model, sd_slab, sd_spike)
  
  hess <- Rcpp_scaled_neg_log_hessian(J, N, SS, sigma_G, rand_par, prior)
  hess_benchmark <- numDeriv::jacobian(function(x) {
    Rcpp_scaled_neg_log_gradient(J, N, SS, sigma_G, parameter_vector_to_list(x), prior)
  }, parameter_list_to_vector(rand_par))
  
  expect_equal(hess, hess_benchmark, tolerance = 1e-6)
})


test_that("scaled_neg_log_hessian is the Hessian of scaled_neg_log_posterior", {
  
  test <- function(log_scale) {
    J <- 20
    N <- 1000
    
    true_par <- random_Gaussian_parameters(J, log_scale = log_scale)
    sd_spike <- runif(1, 1e-3, 1e-1)
    sd_slab <- runif(1, 1e-1, 1e1)
    
    data <- generate_data_MASSIVE_model(N, 2, rep(0.3, J), true_par, log_scale = log_scale)
    
    SS <- data$SS
    sigma_G <- binomial_sigma_G(SS, 2)
    
    rand_par <- random_Gaussian_parameters(J, log_scale = log_scale)
    
    model <- get_random_IV_model(J)
    prior <- decode_IV_model(model, sd_slab, sd_spike)
    
    hess <- scaled_neg_log_hessian(J, N, SS, sigma_G, rand_par, prior, log_scale = log_scale)
    hess_benchmark <- numDeriv::hessian(function(x) {
      scaled_neg_log_posterior(J, N, SS, sigma_G, parameter_vector_to_list(x), prior, log_scale = log_scale)
    }, parameter_list_to_vector(rand_par))
    
    expect_equal(hess, hess_benchmark, tolerance = 1e-6)
  }
  
  test(log_scale = TRUE)
  test(log_scale = FALSE)

})

test_that("scaled_neg_log_hessian is the Jacobian of scaled_neg_log_gradient", {
  
  test <- function(log_scale) {
    J <- 20
    N <- 1000
    
    true_par <- random_Gaussian_parameters(J, log_scale = log_scale)
    sd_spike <- runif(1, 1e-3, 1e-1)
    sd_slab <- runif(1, 1e-1, 1e1)
    
    data <- generate_data_MASSIVE_model(N, 2, rep(0.3, J), true_par, log_scale = log_scale)
    
    SS <- data$SS
    sigma_G <- binomial_sigma_G(SS, 2)
    
    rand_par <- random_Gaussian_parameters(J, log_scale = log_scale)
    
    model <- get_random_IV_model(J)
    prior <- decode_IV_model(model, sd_slab, sd_spike)
    
    hess <- scaled_neg_log_hessian(J, N, SS, sigma_G, rand_par, prior, log_scale = log_scale)
    hess_benchmark <- numDeriv::jacobian(function(x) {
      scaled_neg_log_gradient(J, N, SS, sigma_G, parameter_vector_to_list(x), prior, log_scale = log_scale)
    }, parameter_list_to_vector(rand_par))
    
    expect_equal(hess, hess_benchmark, tolerance = 1e-6)
  }
  
  test(log_scale = TRUE)
  test(log_scale = FALSE)
})

