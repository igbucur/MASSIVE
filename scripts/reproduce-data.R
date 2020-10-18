

# 0. Setup ----------------------------------------------------------------


set.seed(1620)

source('R/utils.R')
Rcpp::sourceCpp('src/MASSIVE_model.cpp', embeddedR = FALSE)
source('R/MASSIVE_model.R')

optimize_posterior_log <- function(J, N, SS, sigma_G, par, prior_sd, nl_post = scaled_nl_posterior_log, nl_grad = scaled_nl_gradient_log) {
  
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
    nl_post(J, N, SS, sigma_G, params, prior_sd)
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
    nl_grad(J, N, SS, sigma_G, params, prior_sd)
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

optimize_posterior_log_cc <- function(J, N, SS, sigma_G, par, prior_sd, nl_post = scaled_nl_posterior_log, nl_grad = scaled_nl_gradient_log) {
  
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
    nl_post(J, N, SS, sigma_G, params, prior_sd)
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
    nl_grad(J, N, SS, sigma_G, params, prior_sd)[-c(2*J+2, 2*J+3)]
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


analyze_MASSIVE_output <- function(MASSIVE_result, J, N, pruning_threshold = 3e-3) {
  MASSIVE_models <- compare_and_select_models(MASSIVE_result, prob_threshold = pruning_threshold)
  beta_samples <- get_LA_posterior_samples(J, N, MASSIVE_models, 1e5)
  plot(stats::density(beta_samples$betas, bw = 0.01))
  analyze_model_approximations(MASSIVE_models)
}


# 4.2 Computing the Approximation - Figure 3 ------------------------------


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

dh <- old_determine_hyperparameters(J, N, SS, est_sigma_G)
mm_prior <- decode_model("1111111111|0001011100|1|1|1", dh$sd_slab, dh$sd_spike)



# Contour posterior -------------------------------------------------------


grid_skappa_X <- seq(-6, 6, 0.1)
grid_skappa_Y <- seq(-6, 6, 0.1)

# Grid confounding coefficients
grid_cc <- expand.grid(grid_skappa_X, grid_skappa_Y)

MASSIVE_posterior <- apply(grid_cc, 1, function(cc) {

  MAP <- optimize_posterior_log_cc(J, N, SS, est_sigma_G, get_ML_solution(SS, cc[1], cc[2], est_sigma_G), mm_prior)
  scaled_nl_posterior_log(J, N, SS, est_sigma_G, MAP, mm_prior)
})


five_optima_posterior_surface <- list()
five_optima_posterior_surface$optima <- list(
  optimize_posterior_log(J, N, SS, est_sigma_G, get_ML_solution(SS, 0, 0, est_sigma_G), mm_prior),
  optimize_posterior_log(J, N, SS, est_sigma_G, get_ML_solution(SS, 3, 3, est_sigma_G), mm_prior),
  optimize_posterior_log(J, N, SS, est_sigma_G, get_ML_solution(SS, -3, -3, est_sigma_G), mm_prior),
  optimize_posterior_log(J, N, SS, est_sigma_G, get_ML_solution(SS, 3, -3, est_sigma_G), mm_prior),
  optimize_posterior_log(J, N, SS, est_sigma_G, get_ML_solution(SS, -3, 3, est_sigma_G), mm_prior)
)
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
  sd_slab_hyperparameter = dh$sd_slab,
  sd_spike_hyperparameter = dh$sd_spike
)

save(five_optima_posterior_surface, file = 'data/five_optima_posterior_surface.RData')

# 5. Empirical Results - Figure 4 -----------------------------------------


prior_comparison_plot

prior_comparison_samples$num_instruments <- NULL
prior_comparison_samples$num_samples <- NULL
prior_comparison_samples$true_par <- NULL
prior_comparison_samples$SS <- NULL
prior_comparison_samples$estimated_sigma_G <- NULL
prior_comparison_samples$true_sd_slab <- NULL
prior_comparison_samples$true_sd_spike <- NULL
prior_comparison_samples$MASSIVE_sd_slab <- NULL
prior_comparison_samples$MASSIVE_sd_spike <- NULL


# 6.1 Determinants of Macroeconomic Growth - Figure 6 ---------------------




library(ivbma)
data(growth)

N <- length(growth$Y)
J <- ncol(growth$Z) + ncol(growth$W) # remove intercept
Jr <- ncol(growth$Z) + 1

Intercept <- growth$W[, 1] + rnorm(nrow(growth$W), 0, 1e-3)

library(ivbma)
data(growth)
ivbma_growth <- ivbma::ivbma(growth$Y, growth$X, growth$Z, growth$W, s = 20000, run.diagnostics = TRUE)
colnames(ivbma_growth$rho) <- c(names(growth$X), names(growth$W))

institutions <- readRDS('../massive/results/final_robust_MASSIVE_growth_institutions_1000.rds')
institutions_models <- compare_and_select_models(institutions, 3e-3)
institutions_samp <- get_LA_posterior_samples(J, N, institutions_models, nrow(ivbma_growth$rho))

integration <- readRDS('../massive/results/final_robust_MASSIVE_growth_integration_1000.rds')
integration_models <- compare_and_select_models(integration, 3e-3)
integration_samp <- get_LA_posterior_samples(J, N, integration_models, nrow(ivbma_growth$rho))

institutions_macroeconomic_growth_samples <- list(
  MASSIVE = get_LA_posterior_samples(J, N, institutions_models, nrow(ivbma_growth$rho)),
  ivbma = ivbma_growth$rho
)

integration_macroeconomic_growth_samples <- list(
  MASSIVE = get_LA_posterior_samples(J, N, integration_models, nrow(ivbma_growth$rho)),
  ivbma = ivbma_growth$rho
)

save(institutions_macroeconomic_growth_samples, integration_macroeconomic_growth_samples,
     file = 'data/determinants_macroeconomic_growth.RData')


# 6.2 Investigating the Relationship between BMI and Psoriasis - Figure 7 ----

library(dplyr)

source('R/utils.R')

# Add suffix to different scenarios
suffix <- "_bmi_psoriasis"

gwas_eaf <- read.table('inst/extdata/GWAS_BMI_EAF.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
gwas_bmi <- read.table('inst/extdata/GWAS_BMI.txt', header = TRUE, stringsAsFactors = FALSE) %>% filter(!grepl("\\*", SNP)) %>%
  inner_join(gwas_eaf, by = "SNP")
# NOTE: we removed rs9581854*, because we could not find the EAF
# We got the EAFs from Locke et al. (European combined)
gwas_psoriasis <- read.table('data/GWAS_psoriasis.txt', header = TRUE, stringsAsFactors = FALSE) %>% filter(!grepl("\\*", SNP))

# Select a smaller number of genetic variants
selected_snps <- which(gwas_bmi$beta_biobank > 0.0)
print(paste(length(selected_snps), "SNPS selected"))

gwas_eaf <- gwas_eaf[selected_snps, ]
gwas_bmi <- gwas_bmi[selected_snps, ]
gwas_psoriasis <- gwas_psoriasis[selected_snps, ]

stopifnot(all(gwas_bmi$SNP == gwas_psoriasis$SNP))


J <- nrow(gwas_bmi)
n <- 2

beta_XG <- gwas_bmi$beta_biobank
sigma_XG <- (gwas_bmi$upr_biobank - gwas_bmi$lwr_biobank) / (qnorm(0.975) - qnorm(0.025))
nobs_XG <- 378274

beta_YG <- log(gwas_psoriasis$beta_biobank)
sigma_YG <- (log(gwas_psoriasis$upr_biobank) - log(gwas_psoriasis$lwr_biobank)) / (qnorm(0.975) - qnorm(0.025))
nobs_YG <- 378274

nobs_XY <- 378274
beta_XY <- log(1.04)
# sigma_XY <- (log(1.53) - log(1.05))/ (qnorm(0.975) - qnorm(0.025))

EAF <- gwas_bmi$EAF
SiG <- sqrt(2 * EAF * (1 - EAF))
# write.table(SiG, file = paste0('~/surfdrive/BayesMR/data/MR_model_SiG', suffix, '.txt'), row.names = FALSE, col.names = FALSE)

beta_IV <- beta_YG / beta_XG
beta_IVW <- weighted.mean(beta_IV, (beta_XG / sigma_YG)^2)

# library(MendelianRandomization) 
# mr_object <- mr_input(bx = beta_XG, bxse = sigma_XG, by = beta_YG, byse = sigma_YG)
# MendelianRandomization::mr_allmethods(mr_object)

# SS_old <- MR_regression_coefficients_to_moments_old(J, beta_XG, sigma_XG, nobs_XG, beta_YG, sigma_YG, nobs_YG, beta_XY, EAF, n)
SS <- MR_regression_coefficients_to_moments(J, beta_XG, sigma_XG, nobs_XG, beta_YG, sigma_YG, nobs_YG, beta_XY, EAF, n)
N <- min(nobs_XG, nobs_YG, nobs_XY)

source('R/MR_MH_sampler.R')
source('R/MR_model.R')
source('R/MR_Laplace_approximation.R')
source('R/MASSIVE.R')

hyp <- determine_hyperparameters(J, N, SS, SiG)
MASSIVE_BMI_psoriasis <- MASSIVE(J, N, SS, SiG, hyp$sd_slab, hyp$sd_spike, max_iter = 5000,
                                 greedy_search = parallel_greedy_search, Laplace_approximation = robust_safe_smart_LA_log)
# 
saveRDS(MASSIVE_BMI_psoriasis, file = "results/final_robust_MASSIVE_BMI_psoriasis.rds")

# prior_sd <- decode_model(get_random_model(J), sd_slab, sd_spike)

#multi_start_MR_Laplace_approximation(J, N, SS, decode_model(get_null_model(J), sd_slab, sd_spike), random_Gaussian_parameters(J))

#print(range(apply(replicate(30, fast_MR_Laplace_approximation(J, N, SS, random_Gaussian_parameters(J), prior_sd)), 2, '[[', 'value')), digits = 20)
#print(range(apply(replicate(10, multi_start_MR_Laplace_approximation(J, N, SS, prior_sd)), 2, '[[', 'value')), digits = 20)

# Rprof(tmp <- tempfile())
# system.time(result <- MR_Laplace_approximation(J, N, SS, random_Gaussian_parameters(J), decode_model(get_null_model(J), sd_slab, sd_spike)))
# Rprof()

# get CI of sbeta
#result$MAP$sbeta 

#LA_cov <- solve(-result$hessian)
#sbeta_var <- LA_cov[2*J+1, 2*J+1]
#c(result$MAP$sbeta + qnorm(0.025) * sqrt(sbeta_var), result$MAP$sbeta, result$MAP$sbeta + qnorm(0.975) * sqrt(sbeta_var))

## Generate data with the same statistics?
# data <- generate_data_MASSIVE_model(N, n, p, par, 1130)
# est_sigma_G <- sqrt(SS[1, 2:(J+1)] * (1 - SS[1, 2:(J+1)] / 2))

# saveRDS(data, file = paste0('results/MR_Laplace_data', suffix, '.rds'))

# SEM <- get_matrix_SEM(J, n, p, gamma, alpha, beta, 0, 0, kappa_X, kappa_Y, sigma_X, sigma_Y)
# true_SS <- get_moments_matrix(J, SEM$mu, SEM$B, SEM$C)

#      user    system   elapsed 
# 33683.543  8671.359 12504.477 
