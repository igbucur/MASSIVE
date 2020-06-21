
#' Routine to perform greedy search at the beginning of MASSIVE algorithm.
#' Here all the neighbors of a model are evaluated and the best one is selected
#' for the next step if it improves the log-evidence.
#'
#' @param J Integer number of candidate instruments
#' @param N Integer number of observations
#' @param SS Numeric moments matrix
#' @param sigma_G Numeric vector of instrument variances
#' @param sd_slab Numeric scale parameter of slab component
#' @param sd_spike Numeric scale parameter of spike component
#' @param init_model Character string describing initial model
#' @param tol Numeric tolerance stopping value for optimization
#' @param post_fun Function computing the negative log-posterior
#' @param get_neighbors Function computing the negative log-posterior gradie
#' @param LA_function Function computing the approximated model log-evidence
#'
#' @return Character string indicating model found using greedy search
parallel_greedy_search <- function(J, N, SS, sigma_G, sd_slab = 1, sd_spike = 0.01, init_model = NULL, get_neighbors = IV_neighbor_models, LA_function = safe_smart_LA_log) {
  
  approximations <- list()
  scores <- list()
  
  if (is.null(init_model)) {
    
    best_start <- which.max(c(
      LA_function(J, N, SS, sigma_G, prior_sd = decode_model(get_full_model(J), sd_slab, sd_spike))$evidence, 
      LA_function(J, N, SS, sigma_G, prior_sd = decode_model(get_ivar_model(J), sd_slab, sd_spike))$evidence
    ))
    
    if (best_start == 1) init_model <- get_full_model(J)
    else if (best_start == 2) init_model <- get_ivar_model(J)
    
  } 
  
  current_model <- init_model
  current_model_LA <- LA_function(J, N, SS, sigma_G, prior_sd = decode_model(current_model, sd_slab, sd_spike))
  
  print(paste("Current best model:", current_model))
  print(paste("Evidence of best model:", current_model_LA$evidence))
  
  approximations[[current_model]] <- current_model_LA
  scores[[current_model]] <- current_model_LA$evidence
  
  models_visited <- 1
  better_model_found <- TRUE
  
  while (better_model_found) {
    
    better_model_found <- FALSE
    
    candidates <- get_neighbors(current_model, J)
    
    # min(parallel::detectCores() / 2, length(candidates))
    new_approximations <- parallel::mclapply(
      candidates, function(new_model) {
        if (is.null(scores[[new_model]])) {
          #models_visited <- models_visited + 1
          #print(paste("New model found during greedy search:", new_model))
          new_model_LA <- LA_function(J, N, SS, sigma_G, prior_sd = decode_model(new_model, sd_slab, sd_spike))
          new_score <- new_model_LA$evidence
          
          return (list(name = new_model, LA = new_model_LA))
        }
        
        NULL
      }, mc.cores = 1)
    
    best_score <- scores[[current_model]]
    
    for (new_eval in new_approximations) {
      if (!is.null(new_eval)) {
        models_visited <- models_visited + 1
        approximations[[new_eval$name]] <- new_eval$LA
        scores[[new_eval$name]] <- new_eval$LA$evidence
        
        if (scores[[new_eval$name]] > best_score) {
          better_model_found <- TRUE
          current_model <- new_eval$name
          current_model_LA <- new_eval$LA
          best_score <- scores[[new_eval$name]]
        }
      }
    }
    
    print(paste("Current best model:", current_model))
    print(paste("Evidence of best model:", current_model_LA$evidence))
    
    # if (any(new_scores > scores[[current_model]]))
  }
  
  list(greedy_model = current_model, approximations = approximations, models_visited = models_visited)
}

#' Routine to perform greedy search at the beginning of MASSIVE algorithm.
#' Here neighbors are examined in a random order until the first one that
#' improves the log-evidence is found and selected.
#'
#' @param J Integer number of candidate instruments
#' @param N Integer number of observations
#' @param SS Numeric moments matrix
#' @param sigma_G Numeric vector of instrument variances
#' @param sd_slab Numeric scale parameter of slab component
#' @param sd_spike Numeric scale parameter of spike component
#' @param init_model Character string describing initial model
#' @param tol Numeric tolerance stopping value for optimization
#' @param post_fun Function computing the negative log-posterior
#' @param get_neighbors Function computing the negative log-posterior gradie
#' @param LA_function Function computing the approximated model log-evidence
#'
#' @return Character string indicating model found using greedy search
stochastic_greedy_search <- function(J, N, SS, sigma_G = NULL, sd_slab = 1, sd_spike = 0.01, init_model = NULL, get_neighbors = IV_neighbor_models, LA_function = safe_smart_LA_log) {
  
  approximations <- list()
  
  if (is.null(init_model)) {
    
    best_start <- which.max(c(
      LA_function(J, N, SS, sigma_G, prior_sd = decode_model(get_full_model(J), sd_slab, sd_spike))$evidence, 
      LA_function(J, N, SS, sigma_G, prior_sd = decode_model(get_ivar_model(J), sd_slab, sd_spike))$evidence
    ))
    
    if (best_start == 1) init_model <- get_full_model(J)
    else if (best_start == 2) init_model <- get_ivar_model(J)
    
  } 
  
  current_model <- init_model
  current_model_LA <- LA_function(J, N, SS, sigma_G, prior_sd = decode_model(current_model, sd_slab, sd_spike))
  
  print(paste("Current best model:", current_model))
  print(paste("Evidence of best model:", current_model_LA$evidence))
  
  approximations[[current_model]] <- current_model_LA
  
  models_visited <- 1
  better_model_found <- TRUE
  
  while (better_model_found) {
    
    better_model_found <- FALSE
    
    candidates <- sample(get_neighbors(current_model, J)) # get_neighbors returned a fixed order, so permute first
    
    for (new_model in candidates) {
      if (is.null(approximations[[new_model]])) {
        models_visited <- models_visited + 1
        
        new_model_LA <- LA_function(J, N, SS, sigma_G, prior_sd = decode_model(new_model, sd_slab, sd_spike))
        approximations[[new_model]] <- new_model_LA
        # print(paste("New model found during greedy search:", new_model))
        # print(paste("Evidence of new model:", new_model_LA$evidence))
      }
      
      if (approximations[[new_model]]$evidence > approximations[[current_model]]$evidence) {
        better_model_found <- TRUE
        current_model <- new_model
        current_model_LA <- new_model_LA
        break
      }
    }
    
    print(paste("Current best model:", current_model))
    print(paste("Evidence of best model:", current_model_LA$evidence))
  }
  
  list(greedy_model = current_model, approximations = approximations, models_visited = models_visited)
}

#' Helper function that returns neighbors of a particular IV model.
#'
#' @param model Character string description of model.
#' @param J Integer number of candidate instruments
#'
#' @return List of model neighbors.
#' 
#' @keywords internal
IV_neighbor_models <- function(model, J) {
  # ids <- c(1:J, (J+2):(2*J+1), 2*J+3, 2*J+5)
  
  lapply((J+2):(2*J+1), function(idx) {
    new_model <- model
    substr(new_model, idx, idx) <- ifelse(substr(new_model, idx, idx) == "1", "0", "1")
    # if (idx == 2*J+5) {
    #   substr(new_model, idx+2, idx+2) <- substr(new_model, idx, idx) # couple confounding coefficients
    # }
    new_model
  })
}