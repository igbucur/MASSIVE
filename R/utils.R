
c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Parameter container conversion ------------------------------------------

#' Function for converting vector of parameters to list format.
#'
#' @param param_vector Vector of MASSIVE model parameters.
#'
#' @return List containing same parameters.
#' @export
parameter_vector_to_list <- function(param_vector) {
  # TODO: maybe mention J explicitly?
  J <- (length(param_vector) - 5) / 2
  list(
    sgamma = param_vector[1:J],
    salpha = param_vector[(J+1):(2*J)],
    sbeta = param_vector[2*J+1],
    skappa_X = param_vector[2*J+2],
    skappa_Y = param_vector[2*J+3],
    sigma_X = param_vector[2*J+4],
    sigma_Y = param_vector[2*J+5]
  )
}

#' Function for converting list of parameters to vector format, often necessary 
#' when running optimization routines.
#'
#' @param param_list List of MASSIVE model parameters. 
#'
#' @return Numeric vector containing same parameters.
#' @export
#'
#' @examples
#' par <- random_Gaussian_parameters(10)
#' parameter_list_to_vector(par)
parameter_list_to_vector <- function(param_list) {
  c(
    param_list$sgamma,
    param_list$salpha,
    param_list$sbeta,
    param_list$skappa_X,
    param_list$skappa_Y,
    param_list$sigma_X,
    param_list$sigma_Y
  )
}

random_Gaussian_parameters_log <- function(J) {
  sgamma <- stats::rnorm(J)
  salpha <- stats::rnorm(J)
  sbeta <- stats::rnorm(1)
  skappa_X <- stats::rnorm(1)
  skappa_Y <- stats::rnorm(1)
  sigma_X <- stats::rnorm(1)
  sigma_Y <- stats::rnorm(1)
  
  list(
    sgamma = sgamma, 
    salpha = salpha, 
    sbeta = sbeta,
    skappa_X = skappa_X,
    skappa_Y = skappa_Y,
    sigma_X = sigma_X,
    sigma_Y = sigma_Y
  )
}

#' Generate random Gaussian distributed parameters
#'
#' @param J Integer number of instrumental variables.
#'
#' @return List of random parameters for the IV model generated from Gaussian distributions.
#' @export
#'
#' @examples random_Gaussian_parameters(10)
random_Gaussian_parameters <- function(J) {
  sgamma <- stats::rnorm(J)
  salpha <- stats::rnorm(J)
  sbeta <- stats::rnorm(1)
  skappa_X <- stats::rnorm(1)
  skappa_Y <- stats::rnorm(1)
  sigma_X <- exp(stats::rnorm(1))
  sigma_Y <- exp(stats::rnorm(1))
  
  list(
    sgamma = sgamma, 
    salpha = salpha, 
    sbeta = sbeta,
    skappa_X = skappa_X,
    skappa_Y = skappa_Y,
    sigma_X = sigma_X,
    sigma_Y = sigma_Y
  )
}

MR_regression_coefficients_to_moments <- function(J, beta_XG, sigma_XG, nobs_XG, beta_YG, sigma_YG, nobs_YG, beta_YX, EAF = rep(0.5, J), n = 2) {
  
  SSS <- matrix(0, J + 3, J + 3) # EV{[1 G X Y] [1 G X Y]^T}
  
  SSS[1, 1] <- 1
  SSS[1, 2:(J+1)] <- SSS[2:(J+1), 1] <- n * EAF
  SSS[2:(J+1), 2:(J+1)] <- n * n * EAF %*% t(EAF) + diag(n * EAF * (1 - EAF), J) # E{GG^T}
  
  SSS[1, J+2] <- SSS[J+2, 1] <- t(beta_XG) %*% SSS[1, 2:(J+1)]
  SSS[1, J+3] <- SSS[J+3, 1] <- t(beta_YG) %*% SSS[1, 2:(J+1)]
  
  SSS[J+2, 2:(J+1)] <- SSS[2:(J+1), J+2] <- SSS[2:(J+1), 2:(J+1)] %*% beta_XG
  SSS[J+3, 2:(J+1)] <- SSS[2:(J+1), J+3] <- SSS[2:(J+1), 2:(J+1)] %*% beta_YG
  
  est_var_X <- stats::median((beta_XG^2 + sigma_XG^2 * nobs_XG) * n * EAF * (1 - EAF))
  est_var_Y <- stats::median((beta_YG^2 + sigma_YG^2 * nobs_YG) * n * EAF * (1 - EAF))
  
  SSS[J+2, J+2] <- est_var_X + (n * t(beta_XG) %*% EAF)^2
  SSS[J+2, J+3] <- SSS[J+3, J+2] <- est_var_X * beta_YX + SSS[1, J+2] * SSS[J+3, 1]
  SSS[J+3, J+3] <- est_var_Y + (n * t(beta_YG) %*% EAF)^2
  
  SSS
}

generate_data_MASSIVE_model_log <- function(N, n, p, par, seed = NULL) {
  par$sigma_X <- exp(par$sigma_X)
  par$sigma_Y <- exp(par$sigma_Y)
  generate_data_MASSIVE_model(N, n, p, par, seed)
}

#' Function to generate data from a MASSIVE generating model.
#'
#' @param N Integer number of observations.
#' @param n Integer number of alleles (trials) for the binomial genetic variables.
#' @param p Numeric expected allele frequencies for the genetic variables.
#' @param par List of MASSIVE model parameters
#' @param seed Integer representing seed for random number generation.
#'
#' @return List containing generate data vector as well as the scatter matrix
#' of first-order and second-order statistics.
#' @export
#'
#' @examples
#' J <- 5 # number of instruments
#' par <- random_Gaussian_parameters(J)
#' generate_data_MASSIVE_model(N = 1000, n = 2, p = rep(0.3, J), par)
generate_data_MASSIVE_model <- function(N, n, p, par, seed = NULL) {
  
  set.seed(seed)
  
  sigma_G <- sqrt(n * p * (1-p))
  gamma <- par$sgamma * par$sigma_X / sigma_G
  alpha <- par$salpha * par$sigma_Y / sigma_G
  beta <- par$sbeta * par$sigma_Y / par$sigma_X
  kappa_X <- par$skappa_X * par$sigma_X
  kappa_Y <- par$skappa_Y * par$sigma_Y
  
  G <- sapply(p, function(t) stats::rbinom(N, n, t)) # genetic variant
  U <- stats::rnorm(N) # confounder
  X <- G %*% gamma + kappa_X * U + stats::rnorm(N, sd = par$sigma_X)
  Y <- G %*% alpha + kappa_Y * U + beta * X + stats::rnorm(N, sd = par$sigma_Y)
  Z <- cbind(1, G, X, Y)
  SS <- t(Z) %*% Z / N # first and second-order moments
  
  list(Z = Z, SS = SS)
}




get_ML_solution <- function(SS, skappa_X, skappa_Y, sigma_G) {
  
  J <- nrow(SS) - 3
  
  cov_XY_G <- SS[(J+2):(J+3), (J+2):(J+3)] - SS[(J+2):(J+3), 2:(J+1)] %*% solve(SS[2:(J+1), 2:(J+1)]) %*% SS[2:(J+1), (J+2):(J+3)]
  
  ML_sigma_X <- sqrt(cov_XY_G[1, 1] / (1 + skappa_X^2))
  ML_sigma_Y <- sqrt((cov_XY_G[1, 1] * cov_XY_G[2, 2] - cov_XY_G[1, 2]^2) * (1 + skappa_X^2) / ((1 + skappa_X^2 + skappa_Y^2) * cov_XY_G[1, 1]))
  ML_sbeta <- (cov_XY_G[1, 2] / (ML_sigma_X * ML_sigma_Y) - skappa_X * skappa_Y) / (1 + skappa_X^2)
  
  ML_sgamma <- solve(SS[2:(J+1), 2:(J+1)], SS[J+2, 2:(J+1)]) / ML_sigma_X * sigma_G
  ML_sGamma <- solve(SS[2:(J+1), 2:(J+1)], SS[J+3, 2:(J+1)]) / ML_sigma_Y * sigma_G

  list(
    sgamma = ML_sgamma,
    salpha = ML_sGamma - ML_sbeta * ML_sgamma,
    sbeta = ML_sbeta,
    skappa_X = skappa_X,
    skappa_Y = skappa_Y,
    sigma_X = ML_sigma_X,
    sigma_Y = ML_sigma_Y
  )
}





# Model codecs and random generators --------------------------------------



get_full_IV_model <- function(J) {
  paste(c(
    paste(rep("1", J), collapse = ""),
    paste(rep("1", J), collapse = ""),
    "1", "1", "1"
  ), collapse = "|")
}

get_empty_IV_model <- function(J) {
  paste(c(
    paste(rep("1", J), collapse = ""),
    paste(rep("0", J), collapse = ""),
    "1", "1", "1"
  ), collapse = "|")
}


#' Function for generating a random IV model (combination of slab and spike priors).
#'
#' @param J Integer number of candidate instruments.
#'
#' @return Character vector representing random IV model.
#' @export
#'
#' @examples
#' get_random_IV_model(5)
get_random_IV_model <- function(J) {
  sp_or_sl <- as.character(sample(c(0, 1), J, replace = TRUE))
  paste(c(
    paste(rep("1", J), collapse = ""),
    paste(sp_or_sl, collapse = ""),
    "1", "1", "1"
  ), collapse = "|")
}

#' Function for decoding a character vector representing IV model.
#'
#' @param code Character vector encoding IV model.
#' @param sd_slab Standard deviation of slab component.
#' @param sd_spike Standard deviation of spike component.
#'
#' @return List describing prior for IV model.
#' @export
#'
#' @examples
#' decode_model(code = get_random_IV_model(5), sd_slab = 1, sd_spike = 0.01)
decode_model <- function(code, sd_slab, sd_spike) {
  tokens <- strsplit(code, "\\|")[[1]]
  IV_model <- stringr::str_extract_all(tokens, "[01]")
  names(IV_model) <- c('sgamma', 'salpha', 'sbeta', 'skappa_X', 'skappa_Y')
  IV_model <- lapply(IV_model, function(x) ifelse(x == "1", sd_slab, sd_spike))
  IV_model$sd_slab <- sd_slab
  IV_model$sd_spike <- sd_spike
  
  IV_model
}

# encode_model <- function(prior) {
#   paste0(c(
#     paste0(sapply(prior$sgamma, function(sd) ifelse(sd == prior$sd_slab, 1, 0)), collapse = ""),
#     paste0(sapply(prior$salpha, function(sd) ifelse(sd == prior$sd_slab, 1, 0)), collapse = ""),
#     ifelse(prior$sbeta == prior$sd_slab, 1, 0),
#     ifelse(prior$skappa_X == prior$sd_slab, 1, 0),
#     ifelse(prior$skappa_Y == prior$sd_slab, 1, 0)
#   ), collapse = "|")
# }

determine_hyperparameters <- function(J, N, SS, sigma_G, fact = 101) {
  
  ML <- get_ML_solution(SS, 0, 0, sigma_G); 

  est_sd_slab <- sqrt(sum(ML$sgamma^2) * fact / J)
  
  # min_sgamma <- sqrt(min(ML$sgamma^2) * fact)
  # P <- N * (dnorm(0, mean = min_sgamma, sd = est_sd_slab, log = TRUE) - dnorm(min_sgamma, mean = min_sgamma, sd = est_sd_slab, log = TRUE))
  # G <- dnorm(min_sgamma, sd = est_sd_spike, log = TRUE) - dnorm(min_sgamma, sd = est_sd_slab, log = TRUE)

  c_root <- stats::uniroot(function(C) {
    (N + 1 - C) * (min(ML$sgamma^2) * fact) / est_sd_slab^2 + log(C) # + 5 to avoid numerical issues? probably fine since just an approximation
  }, interval = c(N+1, 1e12))
  
  est_sd_spike <- est_sd_slab / sqrt(c_root$root)
  
  list(sd_slab = est_sd_slab, sd_spike = est_sd_spike)
}

old_determine_hyperparameters <- function(J, N, SS, sigma_G) {
  
  ML <- get_ML_solution(SS, 0, 0, sigma_G); 
  # ss <- sum(ML$sgamma^2)
  
  est_sd_slab <- sqrt(sum(ML$sgamma^2) / J)

  c_root <- stats::uniroot(function(C) {
    (N + 1 - C) * min(ML$sgamma^2) * (1 + est_sd_slab^2) / est_sd_slab^2 + log(C) # + 5 to avoid numerical issues? probably fine since just an approximation
  }, interval = c(N+1, 1e12))
  
  est_sd_spike <- est_sd_slab / sqrt(c_root$root)
  
  list(sd_slab = est_sd_slab, sd_spike = est_sd_spike)
}


#' Function for deriving the estimated standard deviation of the genetic variants 
#' from the first-order and second-order statistics, assuming these are 
#' binomially distributed.
#'
#' @param SS Numeric matrix of first-order and second-order statistics.
#' @param n Integer number of alleles (trials) for the binomial genetic variables. 
#'
#' @return Numeric vector containing the standard deviations of the genetic variants.
#' @export
#'
#' @examples
#' J <- 5 # number of instruments
#' par <- random_Gaussian_parameters(J)
#' dat <- generate_data_MASSIVE_model(N = 1000, n = 2, p = rep(0.3, J), par)
#' binomial_sigma_G(dat$SS)
binomial_sigma_G <- function(SS, n = 2) {
  J <- nrow(SS) - 3
  sqrt(SS[1, 2:(J+1)] * (1 - SS[1, 2:(J+1)] / n))
}




table_view <- function(list_models, J) {
  
  par_names <- c(
    paste0("sgamma[", 1:J, "]"),
    paste0("salpha[", 1:J, "]"),
    "sbeta",
    "skappa_X",
    "skappa_Y",
    "sigma_X",
    "sigma_Y"
  )
  
  model_names <- names(list_models)
  
  do.call('rbind.data.frame', lapply(1:length(list_models), function(i) {
    
    LA <- list_models[[i]]
    
    y <- do.call('rbind.data.frame', lapply(1:length(LA$optima), function(j) {
      optimum <- LA$optima[[j]]
      
      c(list(model_names[i], LA$evidence, LA$num_optima, j, optimum$origin, optimum$LA,
             optimum$mixture_prob * LA$posterior_probability), as.list(optimum$par))
    }))
    
    names(y) <- c(
      "model",
      "log_evidence",
      "num_optima",
      "id_optimum",
      "at_origin",
      "lev_optimum",
      "probability",
      par_names
    )
    y
  }))
}

