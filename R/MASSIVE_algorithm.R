# source('R/utils.R')
# Rcpp::sourceCpp('src/massive_model.cpp', embeddedR = FALSE)
# source('R/MR_MH_sampler.R')
# source('R/MR_Laplace_approximation.R')

#' Title
#'
#' @param J 
#' @param N 
#' @param SS 
#' @param sd_slab 
#' @param sd_spike 
#' @param max_iter Maximum number of stochastic search steps
#' @param greedy_search Greedy search function
#' @param LA_function Laplace approximation function
#'
#' @return
#' @export
#'
#' @examples
MASSIVE <- function(J, N, SS, sigma_G, sd_slab, sd_spike, max_iter = 1000,
                    greedy_search = parallel_greedy_search, Laplace_approximation = safe_smart_LA_log) {
  
  # WARNING: this must be set before determine hyperparameters
  # if (is.null(sigma_G)) {
  #   # Assuming binomial distribution
  #   sigma_G <- binomial_sigma_G(SS)
  # }
  # 
  # if (is.null(sd_slab) || is.null(sd_spike)) {
  #   est_sd <- determine_hyperparameters(J, N, SS, sigma_G)
  # } else {
  #   est_sd <- list(sd_slab = sd_slab, sd_spike = sd_spike)
  # }

  
  greedy_start <- greedy_search(J, N, SS, sigma_G, sd_slab = sd_slab, sd_spike = sd_spike, LA_function = Laplace_approximation)
  model_approximations <- find_causal_models(J, N, SS, sigma_G, greedy_start = greedy_start, LA_function = Laplace_approximation,
                                             keep_greedy_approximations = TRUE, sd_slab = sd_slab, sd_spike = sd_spike, stop_after = max_iter)
  model_approximations
}

analyze_MASSIVE_output <- function(MASSIVE_result, J, N, pruning_threshold = 3e-3) {
  MASSIVE_models <- compare_and_select_models(MASSIVE_result, prob_threshold = pruning_threshold)
  beta_samples <- get_LA_posterior_samples(J, N, MASSIVE_models, 1e5)
  plot(density(beta_samples$betas, bw = 0.01))
  analyze_model_approximations(MASSIVE_models)
}