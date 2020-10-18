

# 0. Setup ----------------------------------------------------------------

Rcpp::sourceCpp('src/posterior.cpp', embeddedR = FALSE)
file_list <- grep("RcppExports", dir("R", full.names = TRUE), value = TRUE, invert = TRUE)
for (file in file_list) source(file)

library(rstan)
library(ivbma)
library(tidyverse)
library(R2BGLiMS)


# 4.2 Computing the Approximation - Figure 3 ------------------------------

# Wrapper for computing a local optimum over the MASSIVE posterior surface.
optimize_posterior_log <- function(J, N, SS, sigma_G, par, prior_sd) {
  
  optim_MAP <- optim(c(
    par$sgamma, par$salpha, par$sbeta, par$skappa_X, par$skappa_Y, log(par$sigma_X), log(par$sigma_Y)
  ), function(x) {
    params <- list(
      sgamma = x[1:J],
      salpha = x[(J+1):(2*J)],
      sbeta = x[2*J+1],
      skappa_X = x[2*J+2],
      skappa_Y = x[2*J+3],
      sigma_X = x[2*J+4],
      sigma_Y = x[2*J+5]
    )
    scaled_neg_log_posterior(J, N, SS, sigma_G, params, prior_sd)
  }, gr = function(x) {
    params <- list(
      sgamma = x[1:J],
      salpha = x[(J+1):(2*J)],
      sbeta = x[2*J+1],
      skappa_X = x[2*J+2],
      skappa_Y = x[2*J+3],
      sigma_X = x[2*J+4],
      sigma_Y = x[2*J+5]
    )
    scaled_neg_log_gradient(J, N, SS, sigma_G, params, prior_sd)
  }, method = "L-BFGS-B", control = list(maxit = 50000, factr = 1))
  
  list(
    sgamma = optim_MAP$par[1:J],
    salpha = optim_MAP$par[(J+1):(2*J)],
    sbeta = optim_MAP$par[2*J+1],
    skappa_X = optim_MAP$par[2*J+2],
    skappa_Y = optim_MAP$par[2*J+3],
    sigma_X = optim_MAP$par[2*J+4],
    sigma_Y = optim_MAP$par[2*J+5]
  )
}

# Wrapper for computing a local optimum over the MASSIVE posterior surface when the confounding coefficients are fixed.
optimize_posterior_log_cc <- function(J, N, SS, sigma_G, par, prior_sd) {
  
  optim_MAP <- optim(c(
    par$sgamma, par$salpha, par$sbeta, log(par$sigma_X), log(par$sigma_Y)
  ), function(x) {
    params <- list(
      sgamma = x[1:J],
      salpha = x[(J+1):(2*J)],
      sbeta = x[2*J+1],
      skappa_X = par$skappa_X,
      skappa_Y = par$skappa_Y,
      sigma_X = x[2*J+2],
      sigma_Y = x[2*J+3]
    )
    scaled_neg_log_posterior(J, N, SS, sigma_G, params, prior_sd)
  }, gr = function(x) {
    params <- list(
      sgamma = x[1:J],
      salpha = x[(J+1):(2*J)],
      sbeta = x[2*J+1],
      skappa_X = par$skappa_X,
      skappa_Y = par$skappa_Y,
      sigma_X = x[2*J+2],
      sigma_Y = x[2*J+3]
    )
    scaled_neg_log_gradient(J, N, SS, sigma_G, params, prior_sd)[-c(2*J+2, 2*J+3)]
  }, method = "L-BFGS-B", control = list(maxit = 50000, factr = 1))
  
  list(
    sgamma = optim_MAP$par[1:J],
    salpha = optim_MAP$par[(J+1):(2*J)],
    sbeta = optim_MAP$par[2*J+1],
    skappa_X = par$skappa_X,
    skappa_Y = par$skappa_Y,
    sigma_X = optim_MAP$par[2*J+2],
    sigma_Y = optim_MAP$par[2*J+3]
  )
}


set.seed(1620)

J <- 10
n <- 2
N <- 1000

sd_slab <- 1
sd_spike <- 0.01

true_prior <- list()
true_prior$sgamma <- rep(sd_slab, J)
true_prior$salpha <- rep(sd_slab, J)
true_prior$sbeta <- sd_slab
true_prior$skappa_X <- sd_slab
true_prior$skappa_Y <- sd_spike

sgamma <- rnorm(J, sd = true_prior$sgamma)
salpha <- rnorm(J, sd = true_prior$salpha)
sbeta <- rnorm(1, sd = true_prior$sbeta)
skappa_X <- abs(rnorm(1, sd = true_prior$skappa_X))
skappa_Y <- rnorm(1, sd = true_prior$skappa_Y)

p <- runif(J, 0.1, 0.9)
sigma_G <- sqrt(n * p * (1 - p))
sigma_X <- abs(rnorm(1, sd = 10))
sigma_Y <- abs(rnorm(1, sd = 10))


true_par <- list(
  sgamma = sgamma, 
  salpha = salpha, 
  sbeta = sbeta,
  skappa_X = skappa_X,
  skappa_Y = skappa_Y,
  sigma_X = sigma_X,
  sigma_Y = sigma_Y,
  gamma = sgamma * sigma_X / sigma_G,
  alpha = salpha * sigma_Y / sigma_G,
  beta = sbeta * sigma_Y / sigma_X,
  kappa_X = skappa_X * sigma_X,
  kappa_Y = skappa_Y * sigma_Y
)


data <- generate_data_MASSIVE_model(N, n, p, true_par, 1130)
SS <- data$SS
est_sigma_G <- sqrt(SS[1, 2:(J+1)] * (1 - SS[1, 2:(J+1)] / 2))

hyperparameters <- alternative_determine_hyperparameters(J, N, SS, est_sigma_G)
multimodal_prior <- decode_IV_model("1111111111|0001011100|1|1|1", hyperparameters$sd_slab, hyperparameters$sd_spike)

# Create grid for confounding coefficients
grid_skappa_X <- seq(-6, 6, 0.1)
grid_skappa_Y <- seq(-6, 6, 0.1)
grid_cc <- expand.grid(grid_skappa_X, grid_skappa_Y)

# Compute negative log-posterior on that grid
MASSIVE_posterior <- apply(grid_cc, 1, function(cc) {

  MAP <- optimize_posterior_log_cc(J, N, SS, est_sigma_G, get_ML_solution(SS, cc[1], cc[2], est_sigma_G), multimodal_prior)
  scaled_neg_log_posterior(J, N, SS, est_sigma_G, MAP, multimodal_prior)
})

# Create list of results
five_optima_posterior_surface <- list()

# Derive the five optima using appropriate initialization points
five_optima_posterior_surface$optima <- list(
  optimize_posterior_log(J, N, SS, est_sigma_G, get_ML_solution(SS, 0, 0, est_sigma_G), multimodal_prior),
  optimize_posterior_log(J, N, SS, est_sigma_G, get_ML_solution(SS, 3, 3, est_sigma_G), multimodal_prior),
  optimize_posterior_log(J, N, SS, est_sigma_G, get_ML_solution(SS, -3, -3, est_sigma_G), multimodal_prior),
  optimize_posterior_log(J, N, SS, est_sigma_G, get_ML_solution(SS, 3, -3, est_sigma_G), multimodal_prior),
  optimize_posterior_log(J, N, SS, est_sigma_G, get_ML_solution(SS, -3, 3, est_sigma_G), multimodal_prior)
)

# Save grid and posterior values
five_optima_posterior_surface$grid_values <- 
  data.frame(x = grid_cc[, 1], y = grid_cc[, 2], z = MASSIVE_posterior)

# Save some extra information
five_optima_posterior_surface$simulation_details <- list(
  num_instruments = J,
  num_samples = N,
  true_parameters = true_par,
  SS = SS,
  sigma_G = est_sigma_G,
  true_sd_slab = 1,
  true_sd_spike = sd_spike,
  sd_slab_hyperparameter = hyperparameters$sd_slab,
  sd_spike_hyperparameter = hyperparameters$sd_spike
)

save(five_optima_posterior_surface, file = 'data/five_optima_posterior_surface_reproduced.RData')

# 5. Empirical Results - Figure 4 -----------------------------------------

set.seed(2020)

J <- 50
K <- 5
n <- 2
N <- 100000

sd_slab <- 1
sd_spike <- 0.01

true_prior <- list()
true_prior$sgamma <- rep(sd_slab, J)
true_prior$salpha <- c(rep(sd_spike, K), rep(sd_slab, J-K))
true_prior$sbeta <- sd_slab

true_prior$skappa_X <- sd_slab
true_prior$skappa_Y <- sd_slab

sgamma <- rnorm(J, sd = true_prior$sgamma)
salpha <- rnorm(J, sd = true_prior$salpha)
sbeta <- rnorm(1, sd = true_prior$sbeta)
skappa_X <- abs(rnorm(1, sd = true_prior$skappa_X))
skappa_Y <- rnorm(1, sd = true_prior$skappa_Y)

p <- runif(J, 0.1, 0.9)
sigma_X <- abs(rnorm(1, sd = 10))
sigma_Y <- abs(rnorm(1, sd = 10))
sigma_G <- sqrt(n * p * (1 - p))

gamma <- sgamma * sigma_X / sigma_G
alpha <- salpha * sigma_Y / sigma_G
beta <- sbeta * sigma_Y / sigma_X
kappa_X <- skappa_X * sigma_X
kappa_Y <- skappa_Y * sigma_Y

true_par <- list(
  sgamma = sgamma, 
  salpha = salpha, 
  sbeta = sbeta,
  skappa_X = skappa_X,
  skappa_Y = skappa_Y,
  sigma_X = sigma_X,
  sigma_Y = sigma_Y,
  gamma = gamma,
  alpha = alpha,
  beta = beta,
  kappa_X = kappa_X,
  kappa_Y = kappa_Y
)


data <- generate_data_MASSIVE_model(N, n, p, true_par, 1130)
SS <- data$SS
est_sigma_G <- sqrt(SS[1, 2:(J+1)] * (1 - SS[1, 2:(J+1)] / 2))

# Compile MASSIVE model using Stan
stan_MASSIVE <- rstan::stan_model("inst/stan/MASSIVE.stan")

# 1. Gaussian Prior

samples_Gaussian <- rstan::sampling(
  stan_MASSIVE,
  data = list(
    N = N,
    J = J,
    n = 2,
    SS = SS,
    sigma_G = sigma_G,
    mu_X = 0,
    mu_Y = 0,
    sigma_gamma = rep(1, J),
    sigma_alpha = rep(0.01, J)
  ), iter = 5000, cores = 2, control = list(adapt_delta = 0.99, max_treedepth = 20)
)

# 2. Oracle (Spike-and-slab) Prior 

samples_spike_and_slab <- rstan::sampling(
  stan_MASSIVE,
  data = list(
    N = N,
    J = J,
    n = 2,
    SS = SS,
    sigma_G = sigma_G,
    mu_X = 0,
    mu_Y = 0,
    sigma_gamma = true_prior$sgamma,
    sigma_alpha = true_prior$salpha
  ), iter = 5000, cores = 2, control = list(adapt_delta = 0.99, max_treedepth = 20)
)


# 3. MASSIVE Approach

samples_MASSIVE <- MASSIVE(J, N, SS, est_sigma_G, sd_slab = 1, sd_spike = 0.01, max_iter = 5000)


# Collect and save data
prior_comparison_samples <- list(
  Gaussian = as.matrix(samples_Gaussian),
  Oracle = as.matrix(samples_spike_and_slab),
  MASSIVE = samples_MASSIVE,
  simulation_details = list(
    num_instruments = J,
    num_samples = N,
    true_parameters = true_par,
    SS = SS,
    sigma_G = est_sigma_G,
    true_sd_slab = sd_slab,
    true_sd_spike = sd_spike,
    sd_slab_hyperparameter = sd_slab,
    sd_spike_hyperparameter = sd_spike
  )
)

save(prior_comparison_samples, file = 'data/MASSIVE_prior_comparison_reproduced.RData')


# 5. Empirical Results - Figure 5 -----------------------------------------

simulate_JAMMR_data <- function(N, J, K, beta, sigma = 4, seed = 2020, summary = TRUE) {
  set.seed(seed)
  
  n <- 2
  
  p <- runif(J, 0.1, 0.9); sigma_G <- sqrt(n * p * (1 - p))
  
  gamma <- 0.5 + abs(rnorm(J, 0, sd = 0.5))
  alpha <- sample(c(-1,1), J, replace = T) * c(rep(0, K), rnorm(J-K, 0.7, sd = 0.2))
  
  sigma_X <- sigma_Y <- sigma_U <- sigma
  kappa_X <- kappa_Y <- sigma 
  
  sbeta <- beta * sigma_X / sigma_Y
  sgamma <- gamma * sigma_G / sigma_X
  salpha <- alpha * sigma_G / sigma_Y
  skappa_X <- kappa_X / sigma_X
  skappa_Y <- kappa_Y / sigma_Y
  
  true_par <- list(
    sgamma = sgamma, 
    salpha = salpha, 
    sbeta = sbeta,
    skappa_X = skappa_X,
    skappa_Y = skappa_Y,
    sigma_X = sigma_X,
    sigma_Y = sigma_Y,
    gamma = gamma,
    alpha = alpha,
    beta = beta,
    kappa_X = kappa_X,
    kappa_Y = kappa_Y
  )
  
  data <- generate_data_MASSIVE_model(N, n, p, true_par, seed)
  SS <- data$SS
  
  list(
    true_par = true_par,
    SS = data$SS,
    seed = seed,
    N = N,
    J = J,
    K = K,
    beta = beta,
    sigma = sigma,
    Z = data$Z
  )
}

run_JAMMR_simulation <- function(data) {
  
  J <- nrow(data$SS) - 3
  
  df <- data.frame(data$Z[,-1])
  
  names(df) <- c(paste0("G", 1:J), "X", "Y")
  
  EAF <- unlist(
    df %>% 
      dplyr::select(-X, -Y) %>%
      dplyr::summarise_all(base::mean)
  ) / 2
  
  lm_XG <- df %>% 
    dplyr::select(-X, -Y) %>% 
    purrr::map(~lm(df$X ~ .x)) %>% 
    purrr::map(summary) %>%
    purrr::map(c("coefficients"))
  lm_YG <- df %>% 
    dplyr::select(-X, -Y) %>% 
    purrr::map(~lm(df$Y ~ .x)) %>% 
    purrr::map(summary) %>%
    purrr::map(c("coefficients"))
  
  beta_XG <- lm_XG %>% purrr::map_dbl(2)
  sigma_XG <- lm_XG %>% purrr::map_dbl(4)
  beta_YG <- lm_YG %>% purrr::map_dbl(2)
  sigma_YG <- lm_YG %>% purrr::map_dbl(4)
  
  
  N <- nrow(df)
  
  R2BGLiMS::JAMMR(beta_XG, sigma_XG, beta_YG, sigma_YG, N1 = N, eafs = EAF, 
        w = c(0.03 * N, 0.4 * N, 15 * N), jam.seed = 2020)
}

run_MASSIVE_simulation <- function(data) {
  
  J <- nrow(data$SS) - 3
  N <- nrow(data$Z)
  
  hyperparameters <- determine_hyperparameters(J, N, data$SS, binomial_sigma_G(data$SS))
  MASSIVE(J, N, data$SS, binomial_sigma_G(data$SS), hyperparameters$sd_slab, hyperparameters$sd_spike)
}


J <- 10
num_repetitions <- 100

configurations <- expand.grid(
  N = c(1000, 100000),
  sigma = c(1, 4),
  beta = c(0.0, 0.3),
  K = 0:J
)

results <- list()

for (i in 1:nrow(configurations)) {
  N <- configurations$N[i]
  sigma <- configurations$sigma[i]
  beta <- configurations$beta[i]
  K <- configurations$K[i]
  
  results[[i]] <- replicate(num_repetitions, {
    JAMMR_data <- simulate_JAMMR_data(N, J, K, beta, sigma)
    list(
      JAMMR = run_JAMMR_simulation(JAMMR_data),
      MASSIVE = run_MASSIVE_simulation(JAMMR_data)
    )
  })
}

# For JAM-MR the estimate is the field causal
# For MASSIVE we take the median of the posterior samples as the estimate
# The bootstrap RMSE is computed as follows
# rmse <- function(data, indices, beta) sqrt(mean((data[indices] - beta)^2))
# boot_rmse <- boot::boot(estimated_causal_effects, statistic = rmse, R = 1000, beta = beta)
# rmse <- boot_rmse$t0
# rmse_se <- sd(boot_rmse$t)
# rmse_lower_ci <- quantile(boot_rmse$t, probs = 0.025, names = FALSE)
# rmse_upper_ci <- quantile(boot_rmse$t, probs = 0.975, names = FALSE)

# 6.1 Determinants of Macroeconomic Growth - Figure 6 ---------------------

data(growth)

N <- length(growth$Y) # number of samples
J <- ncol(growth$Z) + ncol(growth$W) # number of instruments

# Intercept has to be slightly variable for MASSIVE, not identical to one
Intercept <- growth$W[, 1] + rnorm(nrow(growth$W), 0, 1e-3)

ivbma_growth <- ivbma::ivbma(growth$Y, growth$X, growth$Z, growth$W, s = 20000, run.diagnostics = TRUE)
colnames(ivbma_growth$rho) <- c(names(growth$X), names(growth$W))

# 1. Effect of institutions (rule of law) on growth
institutions_data <- as.matrix(cbind(1, growth$Z, Intercept, growth$W[, -1], growth$X[1], growth$Y))
institutions_SS <- t(institutions_data) %*% institutions_data / N
institutions_sigma_G <- apply(institutions_data[, 2:(J+1)], 2, sd)
institutions_hyperparameters <- determine_hyperparameters(J, N, institutions_SS, institutions_sigma_G)
institutions_MASSIVE_samples <- MASSIVE(
  J, N, institutions_SS, institutions_sigma_G, 
  institutions_hyperparameters$sd_slab, institutions_hyperparameters$sd_spike, 
  max_iter = 5000, greedy_search = parallel_greedy_search
)

institutions_macroeconomic_growth_samples <- list(
  MASSIVE = institutions_MASSIVE_samples,
  ivbma = ivbma_growth$rho
)

# 2. Effect of (economic) integration on growth
integration_data <- as.matrix(cbind(1, growth$Z, Intercept, growth$W[, -1], growth$X[2], growth$Y))
integration_SS <- t(integration_data) %*% integration_data / N
integration_sigma_G <- apply(integration_data[, 2:(J+1)], 2, sd)
integration_hyperparameters <- determine_hyperparameters(J, N, integration_SS, integration_sigma_G)
integration_MASSIVE_samples <- MASSIVE(
  J, N, integration_SS, integration_sigma_G, 
  integration_hyperparameters$sd_slab, integration_hyperparameters$sd_spike,
  max_iter = 5000, greedy_search = parallel_greedy_search)

integration_macroeconomic_growth_samples <- list(
  MASSIVE = integration_MASSIVE_samples,
  ivbma = ivbma_growth$rho
)

save(institutions_macroeconomic_growth_samples, integration_macroeconomic_growth_samples,
     file = 'data/determinants_macroeconomic_growth_reproduced.RData')


# 6.2 Investigating the Relationship between BMI and Psoriasis - Figure 7 ----

# Information in this experiment was collected from Budu-Aggrey et al. (2019)

# NOTE: we removed rs9581854*, because we could not find the EAF
# We got the EAFs from Locke et al. (2015) - European combined (Table 4 in Supplementary Data)
gwas_eaf <- read.table('inst/extdata/GWAS_EAF.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
gwas_bmi <- read.table('inst/extdata/GWAS_BMI.txt', header = TRUE, stringsAsFactors = FALSE) %>% 
  filter(!grepl("\\*", SNP)) %>%
  inner_join(gwas_eaf, by = "SNP")
gwas_psoriasis <- read.table('inst/extdata/GWAS_psoriasis.txt', header = TRUE, stringsAsFactors = FALSE) %>% 
  filter(!grepl("\\*", SNP))

stopifnot(all(gwas_bmi$SNP == gwas_psoriasis$SNP))

J <- nrow(gwas_bmi) # number of instruments
n <- 2 # number of alleles
EAF <- gwas_bmi$EAF # effect allele frequencies
sigma_G <- sqrt(2 * EAF * (1 - EAF)) # standard deviations for genetic variants

# Genetic associations with exposure (BMI)
beta_XG <- gwas_bmi$beta_biobank
sigma_XG <- (gwas_bmi$upr_biobank - gwas_bmi$lwr_biobank) / (qnorm(0.975) - qnorm(0.025))
nobs_XG <- 378274

# Genetic associations with outcome (risk of psoriasis)
beta_YG <- log(gwas_psoriasis$beta_biobank)
sigma_YG <- (log(gwas_psoriasis$upr_biobank) - log(gwas_psoriasis$lwr_biobank)) / (qnorm(0.975) - qnorm(0.025))
nobs_YG <- 378274

# Observational association between exposure (BMI) and outcome (risk of psoriasis)
nobs_XY <- 378274
beta_XY <- log(1.04)

# derive first- and second-order statistics
SS <- MR_regression_coefficients_to_moments(J, beta_XG, sigma_XG, nobs_XG, beta_YG, sigma_YG, nobs_YG, beta_XY, EAF, n)
N <- min(nobs_XG, nobs_YG, nobs_XY) # choose minimum number of samples from all three data sets

# Run MASSIVE algorithm
hyperparameters <- determine_hyperparameters(J, N, SS, sigma_G)
BMI_psoriasis_samples <- MASSIVE(
  J, N, SS, sigma_G, hyperparameters$sd_slab, hyperparameters$sd_spike, 
  max_iter = 5000, greedy_search = parallel_greedy_search)

saveRDS(BMI_psoriasis_samples, file = "results/final_robust_MASSIVE_BMI_psoriasis.rds")

