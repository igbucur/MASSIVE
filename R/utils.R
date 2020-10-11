
c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Parameter container conversion ------------------------------------------

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

gen_data_miv_sem_log <- function(N, n, p, par, seed = NULL) {
  par$sigma_X <- exp(par$sigma_X)
  par$sigma_Y <- exp(par$sigma_Y)
  gen_data_miv_sem(N, n, p, par, seed)
}

gen_data_miv_sem <- function(N, n, p, par, seed = NULL) {
  
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
  ESS <- t(Z) %*% Z / N # first and second-order moments
  
  list(Z = Z, ESS = ESS, U = U)
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

get_null_model <- function(J) {
  paste(c(
    paste(rep("0", J), collapse = ""),
    paste(rep("0", J), collapse = ""),
    "0", "0", "0"
  ), collapse = "|")
}

get_full_model <- function(J) {
  paste(c(
    paste(rep("1", J), collapse = ""),
    paste(rep("1", J), collapse = ""),
    "1", "1", "1"
  ), collapse = "|")
}

get_ivar_model <- function(J) {
  paste(c(
    paste(rep("1", J), collapse = ""),
    paste(rep("0", J), collapse = ""),
    "1", "1", "1"
  ), collapse = "|")
}

get_random_model <- function(J) {
  # note, we bind prior on confounding coefficients (either both slab or spike)
  sp_or_sl <- as.character(sample(c(0, 1), 2 * J + 2, replace = TRUE))
  paste(c(
    paste(sp_or_sl[1:J], collapse = ""),
    paste(sp_or_sl[(J+1):(2*J)], collapse = ""),
    sp_or_sl[2*J+1], sp_or_sl[2*J+2], sp_or_sl[2*J+2]
  ), collapse = "|")
}

get_ply_model <- function(J) {
  sp_or_sl <- as.character(sample(c(0, 1), J, replace = TRUE))
  paste(c(
    paste(rep("1", J), collapse = ""),
    paste(sp_or_sl, collapse = ""),
    "1", "1", "1"
  ), collapse = "|")
}

decode_model <- function(code, sd_slab, sd_spike) {
  tokens <- strsplit(code, "\\|")[[1]]
  digits <- stringr::str_extract_all(tokens, "[01]")
  names(digits) <- c('sgamma', 'salpha', 'sbeta', 'skappa_X', 'skappa_Y')
  digits <- lapply(digits, function(x) ifelse(x == "1", sd_slab, sd_spike))
  digits$sd_slab <- sd_slab
  digits$sd_spike <- sd_spike
  
  digits
}

encode_model <- function(prior) {
  paste0(c(
    paste0(sapply(prior$sgamma, function(sd) ifelse(sd == prior$sd_slab, 1, 0)), collapse = ""),
    paste0(sapply(prior$salpha, function(sd) ifelse(sd == prior$sd_slab, 1, 0)), collapse = ""),
    ifelse(prior$sbeta == prior$sd_slab, 1, 0),
    ifelse(prior$skappa_X == prior$sd_slab, 1, 0),
    ifelse(prior$skappa_Y == prior$sd_slab, 1, 0)
  ), collapse = "|")
}

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

spike_slab_boundary <- function(sd_slab, sd_spike) {
  C <- sd_slab / sd_spike
  sqrt(2*log(C) * sd_slab^2 / (C^2 - 1))
}

binomial_sigma_G <- function(SS, n = 2) {
  J <- nrow(SS) - 3
  sqrt(SS[1, 2:(J+1)] * (1 - SS[1, 2:(J+1)] / n))
}

# Routine for computing smart starting points
smarting_points <- function(SS, sigma_G = binomial_sigma_G(SS)) {
  
  J <- ncol(SS) - 3
  
  # Zero kappa solution
  zero_kappa <- c(0, 0)
  # zero_kappa <- get_ML_solution(SS, 0, 0, sigma_G)
  
  # Zero beta solution
  # also works here, because the means are zero
  cov_XY_G <- SS[(J+2):(J+3), (J+2):(J+3)] - SS[(J+2):(J+3), 2:(J+1)] %*% solve(SS[2:(J+1), 2:(J+1)]) %*% SS[2:(J+1), (J+2):(J+3)]
  #cov_XY <- SS[-1, -1] - SS[-1, 1, drop = FALSE] %*% SS[1, -1, drop = FALSE]
  #cov_XY_G <- cov_XY[(J+1):(J+2), (J+1):(J+2)] - cov_XY[(J+1):(J+2), 1:J] %*% solve(cov_XY[1:J, 1:J]) %*% cov_XY[1:J, (J+1):(J+2)]
  
  corr <- abs(cov_XY_G[1, 2] / sqrt(cov_XY_G[1, 1] * cov_XY_G[2, 2]))
  skappa_X <- sqrt(corr / (1 - corr))
  skappa_Y <- skappa_X * sign(cov_XY_G[1, 2])
  
  zero_beta <- c(skappa_X, skappa_Y)
  # zero_beta <- get_ML_solution(SS, skappa_X, skappa_Y, sigma_G)

  # (Close to) Zero alpha solution
  r_y <- solve(SS[2:(J+1), 2:(J+1)], SS[2:(J+1), J+3])
  r_x <- solve(SS[2:(J+1), 2:(J+1)], SS[2:(J+1), J+2])
  
  # instead of (X, Y), we basically work with (X, Y - beta_ML * X)
  beta_ML <- sum(r_y * r_x) / sum(r_x^2)
  adj_corr <- abs(
    (cov_XY_G[1, 2] - beta_ML * cov_XY_G[1, 1]) /
      sqrt(cov_XY_G[1, 1] * (cov_XY_G[2, 2] + beta_ML^2 * cov_XY_G[1, 1] - 2 * beta_ML * cov_XY_G[1, 2]))
  )
  
  skappa_X <- sqrt(adj_corr / (1 - adj_corr))
  skappa_Y <- skappa_X * sign(cov_XY_G[1, 2] - beta_ML * cov_XY_G[1, 1])
  
  min_alpha <- c(skappa_X, skappa_Y)
  
  # grid_skappa_X <- seq(100, 105, 0.1)
  # grid_skappa_Y <- seq(-1, 1, 0.1)
  # 
  # grid_cc <- expand.grid(grid_skappa_X, grid_skappa_Y)
  # 
  # alpha_min <- function(cc) {
  #   
  #   ML <- get_ML_solution(SS, cc[1], cc[2], sigma_G)
  #   # mean((ML$salpha * ML$sigma_Y)^2)
  #   mean((ML$salpha * ML$sigma_Y)^2)
  # }
  # 
  # alpha_z <- apply(grid_cc, 1, function(cc) {
  #   
  #   ML <- get_ML_solution(SS, cc[1], cc[2], sigma_G)
  #   # mean((ML$salpha * ML$sigma_Y)^2)
  #   mean((ML$salpha * ML$sigma_Y)^2)
  # })
  # 
  # contour(grid_skappa_X, grid_skappa_Y, matrix(alpha_z, length(grid_skappa_X), length(grid_skappa_Y)), nlevels = 50)
  # 
  # f <- function(beta) {
  #   mean((solve(SS[2:(J+1), 2:(J+1)], SS[2:(J+1), J+3]) - beta * solve(SS[2:(J+1), 2:(J+1)], SS[2:(J+1), J+2]))^2)
  # }
  # 
  
  list(zero_kappa, zero_beta, min_alpha)
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


analyze_MASSIVE_output <- function(MASSIVE_result, J, N, pruning_threshold = 3e-3) {
  MASSIVE_models <- compare_and_select_models(MASSIVE_result, prob_threshold = pruning_threshold)
  beta_samples <- get_LA_posterior_samples(J, N, MASSIVE_models, 1e5)
  plot(density(beta_samples$betas, bw = 0.01))
  analyze_model_approximations(MASSIVE_models)
}