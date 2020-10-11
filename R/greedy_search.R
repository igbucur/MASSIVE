
#' Routine to perform greedy search at the beginning of MASSIVE algorithm
#' 
#' @description At each greedy step, we examine all of the IV models in the
#' neighbor set in parallel. The neighbor with the highest model evidence is
#' chosen and the algorithm moves to it if the evidence is higher compared to
#' the current IV model.
#'
#' @param J Integer number of genetic instrumental variables.
#' @param N Integer number of observations.
#' @param SS Numeric matrix containing first- and second-order statistics.
#' @param sigma_G Numeric vector of genetic IV standard deviations.
#' @param sd_slab Numeric value representing slab component standard deviation.
#' @param sd_spike Numeric value representing spike component standard deviation.
#' @param init_model Character vector describing starting IV model in search.
#' @param get_neighbors Function to get neighbor IV models at every step.
#' @param LA_function Function for computing the IV model Laplace approximation
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
#' dat <- gen_data_miv_sem(N, n, EAF, parameters)
#' parallel_greedy_search(J, N, dat$ESS, binomial_sigma_G(dat$ESS), 1, 0.01)
parallel_greedy_search <- function(J, N, SS, sigma_G, sd_slab = 1, sd_spike = 0.01, init_model = NULL, get_neighbors = smart_neighbor_models, LA_function = safe_smart_LA_log) {
  
  approximations <- list()
  scores <- list()
  
  if (is.null(init_model)) {
    
    best_start <- which.max(c(
      LA_function(J, N, SS, sigma_G, par = random_Gaussian_parameters(J), prior_sd = decode_model(get_full_model(J), sd_slab, sd_spike))$evidence, 
      # LA_function(J, N, SS, sigma_G, par = random_Gaussian_parameters(J), prior_sd = decode_model(get_null_model(J), sd_slab, sd_spike))$evidence, 
      LA_function(J, N, SS, sigma_G, par = random_Gaussian_parameters(J), prior_sd = decode_model(get_ivar_model(J), sd_slab, sd_spike))$evidence
    ))
    
    if (best_start == 1) init_model <- get_full_model(J)
    # else if (best_start == 2) init_model <- get_null_model(J)
    else if (best_start == 2) init_model <- get_ivar_model(J)
    
  } 
  
  current_model <- init_model
  current_model_LA <- LA_function(J, N, SS, sigma_G, par = random_Gaussian_parameters(J), prior_sd = decode_model(current_model, sd_slab, sd_spike))
  
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
          new_model_LA <- LA_function(J, N, SS, sigma_G, par = current_model_LA$MAP, prior_sd = decode_model(new_model, sd_slab, sd_spike))
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

#' Greedy search function in the IV model space.
#' 
#' @description At each greedy step, we explore the IV models in the neighbor
#' set one by one. As soon as one has higher evidence than the current IV model,
#' we greedily move to the better neighbor.
#'
#' @param J Integer number of genetic instrumental variables.
#' @param N Integer number of observations.
#' @param SS Numeric matrix containing first- and second-order statistics.
#' @param sigma_G Numeric vector of genetic IV standard deviations.
#' @param sd_slab Numeric value representing slab component standard deviation.
#' @param sd_spike Numeric value representing spike component standard deviation.
#' @param init_model Character vector describing starting IV model in search.
#' @param get_neighbors Function to get neighbor IV models at every step.
#' @param LA_function Function for computing the IV model Laplace approximation
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
#' dat <- gen_data_miv_sem(N, n, EAF, parameters)
#' stochastic_greedy_search(J, N, dat$ESS, binomial_sigma_G(dat$ESS), 1, 0.01)
stochastic_greedy_search <- function(J, N, SS, sigma_G = NULL, sd_slab = 1, sd_spike = 0.01, init_model = NULL, get_neighbors = smart_neighbor_models, LA_function = safe_smart_LA_log) {
  
  approximations <- list()
  
  if (is.null(init_model)) {
    
    best_start <- which.max(c(
      LA_function(J, N, SS, sigma_G, par = random_Gaussian_parameters(J), prior_sd = decode_model(get_full_model(J), sd_slab, sd_spike))$evidence, 
      # LA_function(J, N, SS, sigma_G, par = random_Gaussian_parameters(J), prior_sd = decode_model(get_null_model(J), sd_slab, sd_spike))$evidence, 
      LA_function(J, N, SS, sigma_G, par = random_Gaussian_parameters(J), prior_sd = decode_model(get_ivar_model(J), sd_slab, sd_spike))$evidence
    ))
    
    if (best_start == 1) init_model <- get_full_model(J)
    # else if (best_start == 2) init_model <- get_null_model(J)
    else if (best_start == 2) init_model <- get_ivar_model(J)
    
  } 
  
  current_model <- init_model
  current_model_LA <- LA_function(J, N, SS, sigma_G, par = random_Gaussian_parameters(J), prior_sd = decode_model(current_model, sd_slab, sd_spike))
  
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
        
        new_model_LA <- LA_function(J, N, SS, sigma_G, par = current_model_LA$MAP, prior_sd = decode_model(new_model, sd_slab, sd_spike))
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


neighbor_models <- function(model, J) {
  ids <- c(1:J, (J+2):(2*J+1), 2*J+3, 2*J+5, 2*J+7)
  
  lapply(ids, function(idx) {
    new_model <- model
    substr(new_model, idx, idx) <- ifelse(substr(new_model, idx, idx) == "1", "0", "1")
    new_model
  })
}

smart_neighbor_models <- function(model, J) {
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