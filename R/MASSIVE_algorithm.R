# source('R/utils.R')
# Rcpp::sourceCpp('src/massive_model.cpp', embeddedR = FALSE)
# source('R/MR_MH_sampler.R')
# source('R/MR_Laplace_approximation.R')

#' MASSIVE algorithm (Model Assessment and Stochastic Search for Instrumental Variable Estimation)
#'
#' @param J Integer number of genetic instrumental variables.
#' @param N Integer number of observations.
#' @param SS Numeric matrix containing first- and second-order statistics.
#' @param sigma_G Numeric vector of genetic IV standard deviations.
#' @param sd_slab Numeric value representing slab component standard deviation.
#' @param sd_spike Numeric value representing spike component standard deviation.
#' @param max_iter Maximum number of stochastic search steps
#' @param Laplace_approximation Function used to compute Laplace approximation of IV model.
#' @param greedy_search Greedy search function
#'
#' @return List of explored IV model approximations and their evidences
#' @export
#'
#' @examples
#' J <- 5 # number of instruments
#' N <- 1000 # number of samples
#' parameters <- random_Gaussian_parameters(J) 
#' EAF <- runif(J, 0.1, 0.9) # EAF random values
#' dat <- gen_data_miv_sem(N, n, EAF, parameters)
#' MASSIVE(J, N, dat$ESS, binomial_sigma_G(dat$ESS), 1, 0.01)
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
