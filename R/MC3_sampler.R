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
#' @param ... Extra arguments to pass to Laplace_approximation function.
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
#' dat <- generate_data_MASSIVE_model(N, 2, EAF, parameters)
#' find_causal_models(J, N, dat$SS, binomial_sigma_G(dat$SS), 1, 0.01)
find_causal_models <- function(
  J, N, SS, sigma_G, sd_slab = 1, sd_spike = 0.01, 
  max_iter = 1000, propose_model = propose_neighbor_IV_model,
  LA_function = safe_smart_LA_log, 
  greedy_start = NULL, keep_greedy_approximations = FALSE, ...
) {
  
  if (!is.null(greedy_start)) {
    
    current_model <- greedy_start$greedy_model
    
    if (keep_greedy_approximations) {
      # warning("Incorporating list of Laplace approximations computed during greedy search.")
      approximations <- greedy_start$approximations
      models_visited <- greedy_start$models_visited
    } else {
      current_model_LA <- LA_function(J, N, SS, sigma_G, prior_sd = decode_model(current_model, sd_slab, sd_spike), ...)
      approximations <- list()
      approximations[[current_model]] <- current_model_LA
      models_visited <- 1
    }
  } else {
    warning("Initial model not specified. Starting from random model.")
    current_model <- get_random_IV_model(J)
    current_model_LA <- LA_function(J, N, SS, sigma_G, prior_sd = decode_model(current_model, sd_slab, sd_spike), ...)
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
      new_model_LA <- LA_function(J, N, SS, sigma_G, prior_sd = decode_model(new_model, sd_slab, sd_spike), ...)
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
  
  list(models_visited = models_visited, total_iterations = iter, 
       acceptance_rate = accepted_models / iter, approximations = approximations)
}


#' Proposal function for generating a random neighbor of an IV model.
#'
#' @param IV_model Character vector encoding IV model.
#' @param J Integer number of candidate instruments.
#'
#' @return Character vector encoding random neighbor of IV model.
#' 
#' @keywords internal
propose_neighbor_IV_model <- function(IV_model, J) {
  
  stopifnot(nchar(IV_model) == 2 * J + 7)
  
  idx <- sample((J+2):(2*J+1), 1)
  substr(IV_model, idx, idx) <- ifelse(substr(IV_model, idx, idx) == "1", "0", "1")
  
  IV_model
}


