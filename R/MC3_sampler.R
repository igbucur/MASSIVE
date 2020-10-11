#' Routine to find the most likely causal models using a simple Metropolis-
#' Hastings IV model stochastic search
#'
#' @param J Integer number of genetic instrumental variables.
#' @param N Integer number of observations.
#' @param SS Numeric matrix containing first- and second-order statistics.
#' @param sigma_G Numeric vector of genetic IV standard deviations.
#' @param sd_slab Numeric value representing slab component standard deviation.
#' @param sd_spike Numeric value representing spike component standard deviation.
#' @param stop_after Integer defining maximum number of MCMC iterations.
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
find_causal_models <- function(J, N, SS, sigma_G, sd_slab = 1, sd_spike = 0.01, stop_after = 1000, propose_model = smart_proposal,
                               LA_function = safe_smart_LA_log, greedy_start = NULL, keep_greedy_approximations = FALSE) {
  
  if (!is.null(greedy_start)) {
    
    current_model <- greedy_start$greedy_model
    
    if (keep_greedy_approximations) {
      # warning("Incorporating list of Laplace approximations computed during greedy search.")
      approximations <- greedy_start$approximations
      models_visited <- greedy_start$models_visited
    } else {
      current_model_LA <- LA_function(J, N, SS, sigma_G, par = random_Gaussian_parameters(J), prior_sd = decode_model(current_model, sd_slab, sd_spike))
      approximations <- list()
      approximations[[current_model]] <- current_model_LA
      models_visited <- 1
    }
  } else {
    warning("Initial model not specified. Starting from random model.")
    current_model <- get_ply_model(J)
    # current_model <- get_random_model(J)
    current_model_LA <- LA_function(J, N, SS, sigma_G, par = random_Gaussian_parameters(J), prior_sd = decode_model(current_model, sd_slab, sd_spike))
    approximations <- list()
    approximations[[current_model]] <- current_model_LA
    models_visited <- 1
  }
  
  
  iter <- 0
  iter_since_last_new_model <- 0
  accepted_models <- 0
  
  while(iter < stop_after) {
    
    iter <- iter + 1
    iter_since_last_new_model <- iter_since_last_new_model + 1  
  
    new_model <- propose_model(current_model, J)
    
    if (is.null(approximations[[new_model]])) {
      
      iter_since_last_new_model <- 0
      models_visited <- models_visited + 1
      print(paste("New model found at iteration", iter, ":", new_model))
      new_model_LA <- LA_function(J, N, SS, sigma_G, par = current_model_LA$MAP, prior_sd = decode_model(new_model, sd_slab, sd_spike))
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




get_LA_posterior_samples <- function(J, N, list_models, num_samples = 10000, log_sd = TRUE) {
  
  models <- names(list_models)
  model_num_samples <- round(sapply(list_models, '[[', 'posterior_probability') * num_samples)
  # correct for inaccuracies due to rounding
  model_num_samples[1] <- model_num_samples[1] + num_samples - sum(model_num_samples)
  
  spars <- do.call('rbind', lapply(1:length(model_num_samples), function(i) {
    
    LA <- list_models[[models[i]]]
    LA$optima <- LA$optima[order(sapply(LA$optima, '[[', 'mixture_prob'), decreasing = T)] # TODO: move
    
    mixture_num_samples <- round(sapply(LA$optima, '[[', 'mixture_prob') * model_num_samples[i])
    mixture_num_samples[1] <- mixture_num_samples[1] + model_num_samples[i] - sum(mixture_num_samples)
    
    do.call('rbind', lapply(1:length(mixture_num_samples), function(j) {
      if ((n <- mixture_num_samples[j]) == 0) {
        matrix(0, 0, 2*J+5)
      } else {
        MASS::mvrnorm(n = mixture_num_samples[j], mu = LA$optima[[j]]$par, Sigma = solve(LA$optima[[j]]$hessian * N))
      }
    }))
  }))

  if (log_sd) {
    betas <-  spars[, 2*J+1] / exp(spars[, 2*J+4]) * exp(spars[, 2*J+5])
  } else {
    betas <- spars[, 2*J+1] / spars[, 2*J+4] * spars[, 2*J+5]
  }
  
  
  originating_model <- unlist(sapply(names(model_num_samples), function(model) {
    rep(model, model_num_samples[[model]])
  }))
  
  names(originating_model) <- NULL
  originating_model <- as.factor(originating_model)
  
  list(spars = spars, betas = betas, originating_model = originating_model)
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
  # idx <- sample(c(1:J, (J+2):(2*J+1), 2*J+3, 2*J+5, 2*J+7), 1)
  substr(model, idx, idx) <- ifelse(substr(model, idx, idx) == "1", "0", "1")
  
  # if (idx <= J) {
  #   idx <- idx + J + 1 # also change corresponding alpha
  #   substr(model, idx, idx) <- ifelse(substr(model, idx, idx) == "1", "0", "1")
  # } else if (idx <= 2 * J + 1) {
  #   idx <- idx - J - 1 # also change corresponding gamma
  #   substr(model, idx, idx) <- ifelse(substr(model, idx, idx) == "1", "0", "1")
  # } else 
  
  # if (idx == 2 * J + 5) {
  #   idx <- 2 * J + 7 # also change skappa_Y
  #   substr(model, idx, idx) <- ifelse(substr(model, idx, idx) == "1", "0", "1")
  # } else if (idx == 2 * J + 7) {
  #   idx <- 2 * J + 5 # also change skappa_X
  #   substr(model, idx, idx) <- ifelse(substr(model, idx, idx) == "1", "0", "1")
  # }
  
  model
}


