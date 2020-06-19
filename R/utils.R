
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
  sgamma <- rnorm(J)
  salpha <- rnorm(J)
  sbeta <- rnorm(1)
  skappa_X <- rnorm(1)
  skappa_Y <- rnorm(1)
  sigma_X <- rnorm(1)
  sigma_Y <- rnorm(1)
  
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
  sgamma <- rnorm(J)
  salpha <- rnorm(J)
  sbeta <- rnorm(1)
  skappa_X <- rnorm(1)
  skappa_Y <- rnorm(1)
  sigma_X <- exp(rnorm(1))
  sigma_Y <- exp(rnorm(1))
  
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
  
  est_var_X <- median((beta_XG^2 + sigma_XG^2 * nobs_XG) * n * EAF * (1 - EAF))
  est_var_Y <- median((beta_YG^2 + sigma_YG^2 * nobs_YG) * n * EAF * (1 - EAF))
  
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
  
  G <- sapply(p, function(t) rbinom(N, n, t)) # genetic variant
  U <- rnorm(N) # confounder
  X <- G %*% gamma + kappa_X * U + rnorm(N, sd = par$sigma_X)
  Y <- G %*% alpha + kappa_Y * U + beta * X + rnorm(N, sd = par$sigma_Y)
  Z <- cbind(1, G, X, Y)
  ESS <- t(Z) %*% Z / N # first and second-order moments
  
  list(Z = Z, ESS = ESS, U = U)
}


derive_sufficient_statistics <- function(theta, b21, b31, b32, v1, v2, v3, sc24, sc34, estimate_ss = FALSE, no_samples = NULL, seed = NULL) {
  
  source('utils/covariance_decomposition.R')
  
  l0 <- - sqrt(v1 * theta / (1 - theta))
  l1 <- sqrt(v1 * (1 - theta) / theta)
  
  sb21 <- b21 * sqrt(v1) / sqrt(v2)
  sb31 <- b31 * sqrt(v1) / sqrt(v3)
  sb32 <- b32 * sqrt(v2) / sqrt(v3)
  
  # B <- matrix(c(0, b21, b31, 0, 0, b32, 0, 0, 0), 3, 3)
  # V <- diag(c(v1, v2, v3))
  C <- getC(sc24, sc34, 1)
  rC <- C[2:3, 3, drop = FALSE]
  
  I <- diag(3)
  # Sigma <- true_cov <- getSigma(b21, b31, b32, sc24, c34, v1, v2, v3, 1)
  
  V_TL <- diag(c(v2, v3))
  
  smu_TL <- matrix(c(sb21, sb31 + sb21 * sb32), 2, 1)
  sSigma_TL <- matrix(c(1 + sc24 * sc24, sc24 * sc34 + sb32 * (1 + sc24 * sc24),
                        sc24 * sc34 + sb32 * (1 + sc24 * sc24), 1 + sb32 * sb32 + (sc34 + sb32 * sc24)^2),
                      nrow = 2, ncol = 2)
  
  print(smu_TL)
  print(sSigma_TL)
  
  mu_TL <- sqrt(V_TL) %*% smu_TL %*% solve(sqrt(v1))
  Sigma_TL <- sqrt(V_TL) %*% sSigma_TL %*% sqrt(V_TL)
  
  ss <- list()
  ss$l0 <- l0
  ss$l1 <- l1
  
  
  # Compare estimands with true parameters ----------------------------------
  
  if (estimate_ss) {
    library(mvtnorm)
    set.seed(seed)
    
    if (is.null(no_samples)) {
      print("Number of samples not specified. Setting this number to 1000")
      no_samples <- 1000
    }
    
    
    L <- sample(c(l0, l1), prob = c(1 - theta, theta), size = no_samples, replace = T)
    T_l_0 <- rmvnorm(no_samples, l0 * mu_TL, Sigma_TL)
    T_l_1 <- rmvnorm(no_samples, l1 * mu_TL, Sigma_TL)
    
    data <- data.frame(L = L, T_i = ifelse(L == l0, T_l_0[, 1], T_l_1[, 1]), T_j = ifelse(L == l0, T_l_0[, 2], T_l_1[, 2]))
    data_T <- as.matrix(data[, -1])
    
    ss$L <- mean(data$L)
    ss$LL <- mean(data$L * data$L)
    ss$LT <- matrix(c(mean(data$L * data$T_i), mean(data$L * data$T_j)), 2, 1)
    ss$TT <- t(data_T) %*% data_T / no_samples
    ss$samples <- no_samples
    
  } else {
    
    ss$L <- theta * (l1 - l0) + l0 # \EV{L}
    ss$LL <- ss$L^2 + v1 # \EV{L^2}
    ss$LT <- mu_TL * ss$LL
    ss$TT <- Sigma_TL + (mu_TL %*% t(mu_TL)) * ss$LL
    ss$samples <- Inf
  }
  
  ss
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




quick_approximation <- function(MAP, prior_sd) {
  
  MAP_value <- scaled_nl_posterior_MR(J, N, SS, est_sigma_G, MAP, prior_sd) * N
  log_det <- determinant(scaled_nl_hessian_MR(J, N, SS, est_sigma_G, MAP, prior_sd))$modulus
  
  MAP_value - 0.5 * log_det + (2 * J + 5) / 2 * log(2 * pi)
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

  c_root <- uniroot(function(C) {
    (N + 1 - C) * (min(ML$sgamma^2) * fact) / est_sd_slab^2 + log(C) # + 5 to avoid numerical issues? probably fine since just an approximation
  }, interval = c(N+1, 1e12))
  
  est_sd_spike <- est_sd_slab / sqrt(c_root$root)
  
  list(sd_slab = est_sd_slab, sd_spike = est_sd_spike)
}

old_determine_hyperparameters <- function(J, N, SS, sigma_G) {
  
  ML <- get_ML_solution(SS, 0, 0, sigma_G); 
  # ss <- sum(ML$sgamma^2)
  
  est_sd_slab <- sqrt(sum(ML$sgamma^2) / J)

  c_root <- uniroot(function(C) {
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





compare_and_select_models <- function(discovered_models, prob_threshold = 1e-4) {
  
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


