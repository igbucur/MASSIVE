#' Model Assessment and Stochastic Search for Instrumental Variable Estimation
#'
#' @param J Integer number of candidate instruments
#' @param N Integer number of observations
#' @param SS Numeric moments matrix
#' @param sigma_G Numeric vector of instrument variances
#' @param sd_slab Numeric scale parameter of slab component
#' @param sd_spike Numeric scale parameter of spike component
#' @param max_iter Maximum number of stochastic search steps
#' @param Laplace_approximation Function computing the approximated model log-evidence
#' @param greedy_search Greedy search function
#'
#' @return
#' @export
#'
#' @examples
#' set.seed(2020)
#' J <- 10
#' N <- 10000
#' G <- matrix(rbinom(N * J, 2, 0.3), N, J)
#' U <- rnorm(N)
#' X <- G %*% runif(J, 0.3, 0.5) + U + rnorm(N)
#' Y <- G[, 1:5] %*% runif(5, 0.1, 0.3) + X + U + rnorm(N)
#' 
#' Z <- cbind(1, G, X, Y)
#' SS <- t(Z) %*% Z / N
#' sigma_G <- apply(G, 2, sd)
#' 
#' 
#' samples <- MASSIVE::MASSIVE(J, N, SS, sigma_G, sd_slab = 1, sd_spike = 0.01, max_iter = 1000)
#' plot(density(samples$betas))
#' median(samples$betas)
MASSIVE <- function(J, N, SS, sigma_G, sd_slab, sd_spike, max_iter = 1000,
                    greedy_search = parallel_greedy_search, Laplace_approximation = safe_smart_LA_log,
                    pruning_threshold = 3e-3, posterior_samples = 10000) {
  
  greedy_start <- greedy_search(J, N, SS, sigma_G, sd_slab = sd_slab, sd_spike = sd_spike, LA_function = Laplace_approximation)
  model_approximations <- find_causal_models(J, N, SS, sigma_G, greedy_start = greedy_start, LA_function = Laplace_approximation,
                                             keep_greedy_approximations = TRUE, sd_slab = sd_slab, sd_spike = sd_spike, max_iter = max_iter)
  pruned_list <- prune_model_list(model_approximations, pruning_threshold) # at least thirty samples in each component
  BMA_posterior_samples <- sample_from_BMA_posterior(J, N, pruned_list, posterior_samples)
  BMA_posterior_samples
}

prune_model_list <- function(discovered_models, prob_threshold = 3e-3) {
  
  log_evidences <- sapply(discovered_models$approximations, '[[', 'evidence')
  print(paste("Best model found has log evidence:", max(log_evidences)))
  bayes_factors <- exp(log_evidences - max(log_evidences))
  probabilities <- sort(bayes_factors / sum(bayes_factors), decreasing = T)
  
  # Prune models until the one with the lowest probability makes up for at least 
  # 3e-3 probability, equivalent to 30 out of 10000 samples
  index <- length(probabilities)
  while(probabilities[index] / sum(probabilities[1:(index-1)]) < prob_threshold) {
    index <- index - 1
  }
  
  final_probabilities <- probabilities[1:index] / sum(probabilities[1:index])
  
  promising_models <- discovered_models$approximations[names(final_probabilities)]
  for (model in names(final_probabilities)) {
    promising_models[[model]]$posterior_probability <- final_probabilities[[model]]
  }
  promising_models
}

#' Sample from BMA posterior derived using the MASSIVE approach
#'
#' @param J Integer number of candidate instruments
#' @param N Integer number of observations
#' @param list_models List of IV models found by running MASSIVE
#' @param num_samples 
#' @param log_sd 
#'
#' @return List containing samples from BMA posterior
sample_from_BMA_posterior <- function(J, N, list_models, num_samples = 10000, log_sd = TRUE) {
  
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

analyze_MASSIVE_output <- function(MASSIVE_result, J, N, pruning_threshold = 3e-3) {
  MASSIVE_models <- compare_and_select_models(MASSIVE_result, prob_threshold = pruning_threshold)
  beta_samples <- get_LA_posterior_samples(J, N, MASSIVE_models, 1e5)
  plot(density(beta_samples$betas, bw = 0.01))
  analyze_model_approximations(MASSIVE_models)
}