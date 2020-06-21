safe_smart_LA_log <- function(J, N, SS, sigma_G, prior_sd, 
                              post_fun = scaled_nl_posterior_log, 
                              gr_fun = scaled_nl_gradient_log, 
                              hess_fun = scaled_nl_hessian_log, 
                              opt_fun = find_optimum, starting_points = "smart") {
  
  tryCatch(
    smart_LA_log(J, N, SS, sigma_G, prior_sd, post_fun, gr_fun, hess_fun, opt_fun, starting_points),
    error = function(e) {
      smart_LA_log(J, N, SS, sigma_G, prior_sd, post_fun, gr_fun, hess_fun, opt_fun, "guess")
    }
  )
}

robust_safe_smart_LA_log <- function(J, N, SS, sigma_G, prior_sd, 
                              post_fun = scaled_nl_posterior_log, 
                              gr_fun = scaled_nl_gradient_log, 
                              hess_fun = scaled_nl_hessian_log, 
                              opt_fun = robust_find_optimum, starting_points = "smart") {

  tryCatch(
    smart_LA_log(J, N, SS, sigma_G, prior_sd, post_fun, gr_fun, hess_fun, opt_fun, starting_points),
    error = function(e) {
      smart_LA_log(J, N, SS, sigma_G, prior_sd, post_fun, gr_fun, hess_fun, opt_fun,  "guess")
    }
  )
 #  out
}

#' Routine for computing the Laplace approximation of the MASSIVE posterior
#'
#' @param J Integer number of candidate instruments
#' @param N Integer number of observations
#' @param SS Numeric moments matrix
#' @param sigma_G Numeric vector of instrument variances
#' @param prior_sd List containing prior hyperparameters
#' @param post_fun Function computing the negative log-posterior
#' @param gr_fun Function computing the negative log-posterior gradient
#' @param hess_fun Function computing the negative log-posterior Hessian
#' @param opt_fun Function for optimizing over the negative log-posterior
#' @param starting_points Character string describing the preset optimization starting points list
#'
#' @return List containing mixture of Laplace approximations at the discovered optima.
smart_LA_log <- function(J, N, SS, sigma_G, prior_sd, 
                         post_fun = scaled_nl_posterior_log, 
                         gr_fun = scaled_nl_gradient_log, 
                         hess_fun = scaled_nl_hessian_log, 
                         opt_fun = robust_find_optimum, starting_points = "smart") {
  
  if (starting_points == "guess") {
    starting_points <- list(
      c(prior_sd$sd_spike, prior_sd$sd_spike),
      c(prior_sd$sd_spike, -prior_sd$sd_spike),
      c(2 * prior_sd$sd_slab, 2 * prior_sd$sd_slab),
      c(2 * prior_sd$sd_slab, - 2 * prior_sd$sd_slab),
      c(prior_sd$sd_slab, 2 * prior_sd$sd_slab),
      c(prior_sd$sd_slab, -2 * prior_sd$sd_slab),
      c(2 * prior_sd$sd_slab, -prior_sd$sd_slab),
      c(2 * prior_sd$sd_slab, -prior_sd$sd_slab)
    )
  } else if (starting_points == "smart") {
    starting_points <- smarting_points(SS, sigma_G)
  } else {
    stop("Unknown value for parameter starting_points.")
  }
  
  results <- lapply(starting_points, function(cc) {
    opt_fun(J, N, SS, sigma_G, prior_sd, cc[1], cc[2],
                 post_fun = post_fun, gr_fun = gr_fun, hess_fun = hess_fun)
  })
  
  for (i in 1:length(results)) {
    # force all solutions to have positive skappa_X
    if (results[[i]]$par[2*J+2] < 0) {
      results[[i]]$par[c(2*J+2, 2*J+3)] <- - results[[i]]$par[c(2*J+2, 2*J+3)]
    }
  }
  
  remove_duplicates <- function(results) {
    
    if (length(results) == 1) {
      return (results)
    }
    
    unique_results <- results[1]
    
    for (i in 2:length(results)) {
      r <- results[[i]]
      
      unique <- TRUE
      
      for (j in 1:length(unique_results)) {
        
        u <- unique_results[[j]]
        
        # Various situations could occur
        diff_par <- abs(u$par - r$par)
        
        if (mean(diff_par) < 1e-3 || mean(diff_par[1:J] < 1e-6)) {
          #stopifnot(round(u$value, 3) == round(r$value, 3))
          unique <- FALSE
          break
        }
      }
      
      if (unique) {
        unique_results[[length(unique_results) + 1]] <- r
      } else {
        
        # choose the one with the smaller nl_post value
        if (u$value > r$value) {
          unique_results[[j]] <- r
        }
      }
    }
    
    unique_results
  }  


  # compute the mixture of Laplace approximations
  for (i in 1:length(results)) {
    

    if(is_positive_definite(results[[i]]$hessian)) {
      det_optim <- determinant(results[[i]]$hessian)
      attr(det_optim$modulus, "logarithm") <- NULL
      results[[i]]$LA <- - N * results[[i]]$value - 0.5 * (det_optim$modulus + (2 * J + 5) * log(N)) + (2 * J + 5) / 2 * log(2 * pi)
      
      if (abs(results[[i]]$par[2*J+2]) < prior_sd$sd_spike && abs(results[[i]]$par[2*J+3]) < prior_sd$sd_spike) {
        results[[i]]$origin <- TRUE # the optimum is at the origin
      } else {
        results[[i]]$origin <- FALSE
        results[[i]]$LA <- results[[i]]$LA + log(2) # the optima are symmetric around the origin
      }
    } else {
      # warning("Found saddle point!")
      results[[i]] <- NA
    }
  }
  
  
  results <- results[!is.na(results)] # remove saddle points
  # stopifnot(length(results) > 0)
  if (length(results) == 0) {
    stop("No local optima found!")
    # return(optima = list(), num_optima = 0, evidence = -Inf)
  }
  
  results <- remove_duplicates(results)
  results <- results[order(sapply(results, '[[', 'value'))] # sort results by decreasing posterior value
  
  total_evidence <- matrixStats::logSumExp(sapply(results, '[[', 'LA'))
  
  for (i in 1:length(results)) {
    results[[i]]$mixture_prob <- exp(results[[i]]$LA - total_evidence)
  }
  
  list(optima = results, num_optima = length(results), evidence = total_evidence)
}

is_positive_definite <- function(hessian, tol = 1e-08) {
  eigenvalues <- eigen(hessian, only.values = TRUE)$values
  
  eigenvalues <- ifelse(abs(eigenvalues) < tol, 0, eigenvalues)
  
  all(eigenvalues > 0)
}
