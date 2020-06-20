test_that("scaled_nl_gradient_MR is the gradient of scaled_nl_posterior_MR", {
  
  J <- 20
  N <- 1000
  
  true_par <- random_Gaussian_parameters(J)
  sd_spike <- runif(1, 1e-3, 1e-1)
  sd_slab <- runif(1, 1e-1, 1e1)
  
  # NOTE: fails with single parameter
  data <- gen_data_miv_sem(N, 2, rep(0.3, J), true_par)
  
  SS <- data$ESS
  sigma_G <- binomial_sigma_G(SS, 2)
  
  rand_par <- random_Gaussian_parameters(J)
  
  model <- get_random_model(J)
  prior <- decode_model(model, sd_slab, sd_spike)
  
  grad <- c(scaled_nl_gradient_MR(J, N, SS, sigma_G, rand_par, prior))
  grad_benchmark <- numDeriv::grad(function(x) {
    scaled_nl_posterior_MR(J, N, SS, sigma_G, parameter_vector_to_list(x), prior)
  }, parameter_list_to_vector(rand_par))
  
  expect_equal(grad, grad_benchmark, tolerance = 1e-6)
})


test_that("scaled_nl_gradient_log is the gradient of scaled_nl_posterior_log", {
  
  J <- 20
  N <- 1000
  
  true_par <- random_Gaussian_parameters_log(J)
  sd_spike <- runif(1, 1e-3, 1e-1)
  sd_slab <- runif(1, 1e-1, 1e1)
  
  # NOTE: fails with single parameter
  data <- gen_data_miv_sem_log(N, 2, rep(0.3, J), true_par)
  
  SS <- data$ESS
  sigma_G <- binomial_sigma_G(SS, 2)
  
  rand_par <- random_Gaussian_parameters_log(J)
  
  model <- get_random_model(J)
  prior <- decode_model(model, sd_slab, sd_spike)
  
  grad <- c(scaled_nl_gradient_log(J, N, SS, sigma_G, rand_par, prior))
  grad_benchmark <- numDeriv::grad(function(x) {
    scaled_nl_posterior_log(J, N, SS, sigma_G, parameter_vector_to_list(x), prior)
  }, parameter_list_to_vector(rand_par))
  
  expect_equal(grad, grad_benchmark, tolerance = 1e-6)
})
