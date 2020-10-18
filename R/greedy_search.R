
#' Routine to perform greedy search at the beginning of MASSIVE algorithm
#' 
#' @description At each greedy step, we examine all of the IV models in the
#' neighbor set in parallel. The neighbor with the highest model evidence is
#' chosen and the algorithm moves to it if the evidence is higher compared to
#' the current IV model.
#'
#' @param J Integer number of candidate instruments..
#' @param N Integer number of observations.
#' @param SS Numeric matrix containing first- and second-order statistics.
#' @param sigma_G Numeric vector of genetic IV standard deviations.
#' @param sd_slab Numeric scale parameter of slab component.
#' @param sd_spike Numeric scale parameter of spike component.
#' @param init_model Character vector describing starting IV model in search.
#' @param get_neighbors Function to get neighbor IV models at every step.
#' @param LA_function Function for computing the IV model Laplace approximation
#' @param ... Extra arguments to pass to Laplace approximation function.
#'
#' @return A list containing the optimum found with greedy search, the list of 
#' IV models visited and their approximations, and the number of visited models.
#' 
#' @export
#'
#' @examples
#' J <- 5 # number of instruments
#' N <- 1000 # number of samples
#' parameters <- random_Gaussian_parameters(J) 
#' EAF <- runif(J, 0.1, 0.9) # EAF random values
#' dat <- generate_data_MASSIVE_model(N, 2, EAF, parameters)
#' parallel_greedy_search(J, N, dat$SS, binomial_sigma_G(dat$SS), 1, 0.01)
parallel_greedy_search <- function(
  J, N, SS, sigma_G, sd_slab = 1, sd_spike = 0.01, 
  init_model = NULL, get_neighbors = neighbor_IV_models, 
  LA_function = safe_Laplace_approximation, ...
) {
  
  approximations <- list()
  scores <- list()
  
  # If no initial IV model given, start from either the empty or full IV model
  if (is.null(init_model)) {
    
    best_start <- which.max(c(
      LA_function(J, N, SS, sigma_G, prior_sd = decode_IV_model(get_full_IV_model(J), sd_slab, sd_spike), ...)$evidence, 
      LA_function(J, N, SS, sigma_G, prior_sd = decode_IV_model(get_empty_IV_model(J), sd_slab, sd_spike), ...)$evidence
    ))
    
    if (best_start == 1) init_model <- get_full_IV_model(J)
    else if (best_start == 2) init_model <- get_empty_IV_model(J)
    
  } 
  
  current_model <- init_model
  current_model_LA <- LA_function(J, N, SS, sigma_G, prior_sd = decode_IV_model(current_model, sd_slab, sd_spike), ...)
  
  print(paste("Current best model:", current_model))
  print(paste("Evidence of best model:", current_model_LA$evidence))
  
  approximations[[current_model]] <- current_model_LA
  scores[[current_model]] <- current_model_LA$evidence
  
  models_visited <- 1
  better_model_found <- TRUE
  
  while (better_model_found) {
    
    better_model_found <- FALSE
    
    candidates <- get_neighbors(current_model, J)
    
    new_approximations <- parallel::mclapply(
      candidates, function(new_model) {
        if (is.null(scores[[new_model]])) {
          new_model_LA <- LA_function(J, N, SS, sigma_G, prior_sd = decode_IV_model(new_model, sd_slab, sd_spike), ...)
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
    
  }
  
  list(greedy_model = current_model, approximations = approximations, models_visited = models_visited)
}

#' Greedy search function in the IV model space.
#' 
#' @description At each greedy step, we explore the IV models in the neighbor
#' set one by one. As soon as one has higher evidence than the current IV model,
#' we greedily move to the better neighbor.
#'
#' @param J Integer number of candidate instruments..
#' @param N Integer number of observations.
#' @param SS Numeric matrix containing first- and second-order statistics.
#' @param sigma_G Numeric vector of genetic IV standard deviations.
#' @param sd_slab Numeric scale parameter of slab component.
#' @param sd_spike Numeric scale parameter of spike component.
#' @param init_model Character vector describing starting IV model in search.
#' @param get_neighbors Function to get neighbor IV models at every step.
#' @param LA_function Function for computing the IV model Laplace approximation
#' @param ... Extra arguments to pass to Laplace approximation function.
#'
#' @return A list containing the optimum found with greedy search, the list of 
#' IV models visited and their approximations, and the number of visited models.
#' 
#' @export
#'
#' @examples
#' J <- 5 # number of instruments
#' N <- 1000 # number of samples
#' parameters <- random_Gaussian_parameters(J) 
#' EAF <- runif(J, 0.1, 0.9) # EAF random values
#' dat <- generate_data_MASSIVE_model(N, 2, EAF, parameters)
#' stochastic_greedy_search(J, N, dat$SS, binomial_sigma_G(dat$SS), 1, 0.01)
stochastic_greedy_search <- function(
  J, N, SS, sigma_G = NULL, sd_slab = 1, sd_spike = 0.01, 
  init_model = NULL, get_neighbors = neighbor_IV_models, 
  LA_function = safe_Laplace_approximation, ...
  ) {
  
  approximations <- list()
  
  # If no initial IV model given, start from either the empty or full IV model
  if (is.null(init_model)) {
    
    best_start <- which.max(c(
      LA_function(J, N, SS, sigma_G, prior_sd = decode_IV_model(get_full_IV_model(J), sd_slab, sd_spike), ...)$evidence, 
      LA_function(J, N, SS, sigma_G, prior_sd = decode_IV_model(get_empty_IV_model(J), sd_slab, sd_spike), ...)$evidence
    ))
    
    if (best_start == 1) init_model <- get_full_IV_model(J)
    else if (best_start == 2) init_model <- get_empty_IV_model(J)
  } 
  
  current_model <- init_model
  current_model_LA <- LA_function(J, N, SS, sigma_G, prior_sd = decode_IV_model(current_model, sd_slab, sd_spike), ...)
  
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
        
        new_model_LA <- LA_function(J, N, SS, sigma_G, prior_sd = decode_IV_model(new_model, sd_slab, sd_spike), ...)
        approximations[[new_model]] <- new_model_LA
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

# 
# neighbor_models <- function(model, J) {
#   ids <- c(1:J, (J+2):(2*J+1), 2*J+3, 2*J+5, 2*J+7)
#   
#   lapply(ids, function(idx) {
#     new_model <- model
#     substr(new_model, idx, idx) <- ifelse(substr(new_model, idx, idx) == "1", "0", "1")
#     new_model
#   })
# }

#' Function for returning the neighboring IV model, i.e., those for which a
#' signel component changes from 0 (slab) to 1 (spike) or vice versa.
#'
#' @param model Character vector describing IV model.
#' @param J Integer number of candidate instruments.
#'
#' @return List of IV models neighboring (one component changes) given IV model.
#' 
#' @keywords internal
neighbor_IV_models <- function(model, J) {
  
  lapply((J+2):(2*J+1), function(idx) {
    new_model <- model
    substr(new_model, idx, idx) <- ifelse(substr(new_model, idx, idx) == "1", "0", "1")
    new_model
  })
}
