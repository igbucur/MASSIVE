test_that("scaled_neg_log_posterior is equal to of Rcpp_scaled_neg_log_posterior", {
  
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
  
  post <- scaled_neg_log_posterior(J, N, SS, sigma_G, rand_par, prior, log_scale = FALSE)
  post_benchmark <- Rcpp_scaled_neg_log_posterior(J, N, SS, sigma_G, rand_par, prior)
  
  expect_equal(post, post_benchmark, tolerance = 1e-6)
  
})