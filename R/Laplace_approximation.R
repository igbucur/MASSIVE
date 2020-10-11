

#' Routine for computing the Laplace approximation of the MASSIVE posterior
#'
#' @param J Integer number of genetic instrumental variables.
#' @param N Integer number of observations.
#' @param SS Numeric matrix containing first- and second-order statistics.
#' @param sigma_G Numeric vector of genetic IV standard deviations.
#' @param prior_sd List of standard deviations for the parameter Gaussian priors.
#' @param post_fun Function for computing the IV model posterior value.
#' @param gr_fun Function for computing the IV model posterior gradient.
#' @param hess_fun Function for computing the IV model posterior Hessian.
#' @param opt_fun Function for finding the IV model posterior optima.
#' @param starting_points Character vector indicating how to pick the starting
#' points: "smart" or "guess" strategy?
#'
#' @return A list containing the Laplace approximation.
#' \itemize{
#'   \item optima - List of optima found.
#'   \item num_optima - Number of optima found on the posterior surface.
#'   \item evidence - Numeric value of total approximation model evidence.
#' }
#' 
#' @export
#'
#' @examples
#' J <- 5 # number of instruments
#' N <- 1000 # number of samples
#' parameters <- random_Gaussian_parameters(J) 
#' EAF <- runif(J, 0.1, 0.9) # EAF random values
#' dat <- gen_data_miv_sem(N, n, EAF, parameters)
#' safe_smart_LA_log(J, N, dat$ESS, binomial_sigma_G(dat$ESS), decode_model(get_ply_model(J), 1, 0.01))
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
  #  out
}

robust_safe_smart_LA_log <- function(J, N, SS, sigma_G, prior_sd, 
                              post_fun = scaled_nl_posterior_log, 
                              gr_fun = scaled_nl_gradient_log, 
                              hess_fun = scaled_nl_hessian_log, 
                              opt_fun = robust_find_optimum, starting_points = "smart") {

  tryCatch(
    smart_LA_log(J, N, SS, sigma_G, prior_sd, post_fun, gr_fun, hess_fun, opt_fun, starting_points),
    error = function(e) {
      smart_LA_log(J, N, SS, sigma_G, prior_sd, post_fun, gr_fun, hess_fun, opt_fun, "guess")
    }
  )
 #  out
}


smart_LA_log <- function(J, N, SS, sigma_G, prior_sd, 
                         post_fun = scaled_nl_posterior_log, 
                         gr_fun = scaled_nl_gradient_log, 
                         hess_fun = scaled_nl_hessian_log, 
                         opt_fun = robust_find_optimum, starting_points = "smart") {
  
  # TODO: maybe add optimization using 'par' as starting point
  
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
        # print(str(u))
        
        # found a scenario where the first two statements did not hold, but it seemed to be the same minimum
        
        # if (round(u$par[2*J+2], 3) == round(r$par[2*J+2], 3) &&
        #     round(u$par[2*J+3], 3) == round(r$par[2*J+3], 3) &&
        #     round(u$value, 3) == round(r$value, 3)) {
        #   unique <- FALSE
        #   break
        # }
        
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
        
        
        # sometimes two optima are found very close to each other with significantly different hessians
        # in this case we pick the one with less curvature (higher evidence), since it is probably less
        # prone to numerical errors and encompasses the other optimum as well
        # if (u$LA < r$LA) {
        #   unique_results[[j]] <- r
        # }
        # unique_results[[j]]$LA <- max(u$LA, r$LA)
      }
    }
    
    unique_results
  }  

  is_posdef <- function(h, tol = 1e-08) {
    eigenvalues <- eigen(h, only.values = TRUE)$values
    n <- nrow(h)
    for (i in 1:n) {
      if (abs(eigenvalues[i]) < tol) {
        eigenvalues[i] <- 0
      }
    }
    if (any(eigenvalues <= 0)) {
      return(FALSE)
    }
    return(TRUE)
  }
  
  # print(length(results))
  
  # compute the mixture of Laplace approximations
  for (i in 1:length(results)) {
    
    # symmetrize again because is.positive.definite is stupid
    # results[[i]]$hessian <- (results[[i]]$hessian + t(results[[i]]$hessian)) / 2
    # print(is_posdef(results[[i]]$hessian, 1e-8))
    
    if(is_posdef(results[[i]]$hessian, 1e-8)) {
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
  stopifnot(length(results) > 0)
  # if (length(results) == 0) {
  #   error("No optima found")
  # return(optima = list(), num_optima = 0, evidence = -Inf)
  # }
  
  
  results <- remove_duplicates(results)
  results <- results[order(sapply(results, '[[', 'value'))] # sort results by decreasing posterior value
  
  total_evidence <- matrixStats::logSumExp(sapply(results, '[[', 'LA'))
  
  for (i in 1:length(results)) {
    results[[i]]$mixture_prob <- exp(results[[i]]$LA - total_evidence)
  }
  
  # posterior_Gaussian_MR(J, N, SS, MAP, prior_sd, n)
  list(optima = results, num_optima = length(results), evidence = total_evidence)
}
