#' Routine to find the most likely causal models using a simple Metropolis-
#' Hastings IV model stochastic search
#'
#' @param J Integer number of candidate instruments.
#' @param N Integer number of observations.
#' @param SS Numeric matrix containing first- and second-order statistics.
#' @param sigma_G Numeric vector of genetic IV standard deviations.
#' @param sd_slab Numeric scale parameter of slab component.
#' @param sd_spike Numeric scale parameter of spike component.
#' @param max_iter Integer defining maximum number of MCMC iterations.
#' @param propose_model Proposal function indicating new IV model candidates.
#' @param LA_function Function for computing the IV model Laplace approximation.
#' @param greedy_start List returned by greedy search to be used as starting point.
#' @param keep_greedy_approximations Logical flag asking if greedy search history
#' should be stored and used in the new MCMC search.
#'
#' @return A list containing the list of models explored and their approximations,
#' the number of models explored, the total number of MCMC iterations, as well
#' as the acceptance rate.
#' @export
#'
#' @examples
#' J <- 5 # number of instruments
#' N <- 1000 # number of samples
#' parameters <- random_Gaussian_parameters(J) 
#' EAF <- runif(J, 0.1, 0.9) # EAF random values
#' dat <- gen_data_miv_sem(N, n, EAF, parameters)
#' find_causal_models(J, N, dat$ESS, binomial_sigma_G(dat$ESS), 1, 0.01)
find_causal_models <- function(J, N, SS, sigma_G, sd_slab = 1, sd_spike = 0.01, max_iter = 1000, propose_model = smart_proposal,
                               LA_function = safe_smart_LA_log, greedy_start = NULL, keep_greedy_approximations = FALSE) {
  
  if (!is.null(greedy_start)) {
    
    current_model <- greedy_start$greedy_model
    
    if (keep_greedy_approximations) {
      # warning("Incorporating list of Laplace approximations computed during greedy search.")
      approximations <- greedy_start$approximations
      models_visited <- greedy_start$models_visited
    } else {
      current_model_LA <- LA_function(J, N, SS, sigma_G, prior_sd = decode_model(current_model, sd_slab, sd_spike))
      approximations <- list()
      approximations[[current_model]] <- current_model_LA
      models_visited <- 1
    }
  } else {
    warning("Initial model not specified. Starting from random model.")
    current_model <- get_ply_model(J)
    current_model_LA <- LA_function(J, N, SS, sigma_G, prior_sd = decode_model(current_model, sd_slab, sd_spike))
    approximations <- list()
    approximations[[current_model]] <- current_model_LA
    models_visited <- 1
  }
  
  
  iter <- 0
  iter_since_last_new_model <- 0
  accepted_models <- 0
  
  while(iter < max_iter) {
    
    iter <- iter + 1
    iter_since_last_new_model <- iter_since_last_new_model + 1  
  
    new_model <- propose_model(current_model, J)
    
    if (is.null(approximations[[new_model]])) {
      
      iter_since_last_new_model <- 0
      models_visited <- models_visited + 1
      print(paste("New model found at iteration", iter, ":", new_model))
      new_model_LA <- LA_function(J, N, SS, sigma_G, prior_sd = decode_model(new_model, sd_slab, sd_spike))
      approximations[[new_model]] <- new_model_LA
    }
    
    log_bf <- approximations[[new_model]]$evidence - approximations[[current_model]]$evidence
    # print(log_bf)
    acceptance_prob <- min(1, exp(log_bf))
    
    if (stats::runif(1) < acceptance_prob) {
      accepted_models <- accepted_models + 1
      current_model <- new_model
    }
  }
  
  list(models_visited = models_visited, total_iterations = iter, acceptance_rate = accepted_models / iter, approximations = approximations)
}


simple_proposal <- function(model, J) {
  
  stopifnot(nchar(model) == 2 * J + 7)
  
  # correct for separator
  idx <- sample(c(1:J, (J+2):(2*J+1), 2*J+3, 2*J+5, 2*J+7), 1)
  
  substr(model, idx, idx) <- ifelse(substr(model, idx, idx) == "1", "0", "1")
  
  model
}

smart_proposal <- function(model, J) {
  
  stopifnot(nchar(model) == 2 * J + 7)
  
  # correct for separator
  idx <- sample((J+2):(2*J+1), 1)
  substr(model, idx, idx) <- ifelse(substr(model, idx, idx) == "1", "0", "1")
  
  
  model
}


