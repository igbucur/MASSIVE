

# 0. Setup and data preprocessing -----------------------------------------

library(dplyr)

source('R/utils.R')

# Add suffix to different scenarios
suffix <- "_bmi_psoriasis"

gwas_eaf <- read.table('data/GWAS_BMI_EAF.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
gwas_bmi <- read.table('data/GWAS_BMI.txt', header = TRUE, stringsAsFactors = FALSE) %>% filter(!grepl("\\*", SNP)) %>%
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



# ICML plot ---------------------------------------------------------------

library(ggplot2)
library(bayesplot)

BMI_psoriasis_result <- readRDS('results/final_robust_MASSIVE_BMI_psoriasis.rds')
mod <- compare_and_select_models(BMI_psoriasis_result, 3e-3)
samples <- get_LA_posterior_samples(J, N, mod, 10000)

ggplot(data.frame(beta = samples$betas), aes(x = beta)) + stat_density() + xlim(c(-0.15, 0.15)) +
  ggplot2::geom_vline(xintercept = log(c(1.08, 1.04, 1.13)), linetype = c(1, 2, 2))
  

bayesplot::mcmc_areas(data.frame(beta = samples$betas[samples$betas > -10]), pars = "beta") +
  ggplot2::geom_vline(xintercept = log(c(1.08, 1.04, 1.13)), linetype = c(1, 2, 2)) +
  ggplot2::xlim(c(-0.15, 0.15))

ggplot2::ggsave("../tex/MASSIVE/Camera-Ready Version/plot_beta_BMI_psoriasis.pdf", width = 5, height = 4)

# Use MR_MH_Sampler -------------------------------------------------------
  
source('R/MR_MH_sampler.R')
source('R/write_MR_Gaussian_configuration_file.R')

smart_MR_Laplace_approximation(J, N, SS, decode_model(get_null_model(J), sd_slab, sd_spike))$evidence
smart_MR_Laplace_approximation(J, N, SS, decode_model(get_full_model(J), sd_slab, sd_spike))$evidence
smart_MR_Laplace_approximation(J, N, SS, decode_model(get_ivar_model(J), sd_slab, sd_spike))$evidence

#greedy_start <- stochastic_greedy_search(J, N, SS, LA_function = fast_MR_Laplace_approximation, init_model = get_full_model(J), sd_slab = sd_slab, sd_spike = sd_spike)
#greedy_start_ivar <- stochastic_greedy_search(J, N, SS, LA_function = fast_MR_Laplace_approximation, init_model = get_ivar_model(J), sd_slab = sd_slab, sd_spike = sd_spike)
# greedy_start_smart <- stochastic_greedy_search(J, N, SS, init_model = get_ivar_model(J), sd_slab = sd_slab, sd_spike = sd_spike)
# saveRDS(greedy_start_smart, file = 'data/BMI_psoriasis_stochastic_greedy_start_smart_ivar_slab_1e-1_spike_1e-3.rds')
# greedy_start <- readRDS('data/BMI_psoriasis_stochastic_greedy_start_null.rds')

greedy_start_smart <- readRDS('data/BMI_psoriasis_stochastic_greedy_start_smart_ivar_slab_1e-1_spike_1e-3.rds')

# TODO: change greedy start function to f
system.time(MH_result_BMI_2 <- find_causal_models(J, N, SS, greedy_start = greedy_start_smart, keep_greedy_approximations = TRUE,
                                                sd_slab = sd_slab, sd_spike = sd_spike, stop_after = 2000))


x <- compare_and_select_models(MH_result_BMI)
y <- get_LA_posterior_samples(J, N, x)
plot(density(y$betas))

# load('results/BMI_psoriasis_SGS_null_MH_1000.rdata')

save(greedy_start_smart, MH_result_BMI, sd_spike, sd_slab, file = 'results/BMI_psoriasis_SGS_smart_ivar_MH_1000_slab_1e-1_spike_1e-3.rdata')

  rank_models(MH_result$samples)
analyze_model_samples(MH_result$samples)

x

system.time(multi_MH_result <- multi_start_MR_MH_sampler(J, N, SS, propose_model = simple_proposal, sd_slab = sd_slab, sd_spike = sd_spike, iter = 100000, init_model = get_null_model(J), TRUE))


Rprof(tmp <- tempfile())
system.time(smart_MH_result <- MR_MH_sampler(J, N, SS, init_model = get_null_model(J), propose_model = smart_proposal, sd_slab = sd_slab, sd_spike = sd_spike, iter = 100000))
system.time(sparse_MH_result <- MR_MH_sampler(J, N, SS, init_model = get_null_model(J), propose_model = sparse_proposal, sd_slab = sd_slab, sd_spike = sd_spike, iter = 100000))

Rprof()
summaryRprof(tmp)$by.self

system.time(MH_result_7 <- MR_MH_sampler(J, N, SS, sd_slab = sd_slab, propose_model = simple_proposal, sd_spike = sd_spike, iter = 1e7, init_model = get_null_model(J), TRUE))
system.time(multi_MH_result_7 <- multi_start_MR_MH_sampler(J, N, SS, propose_model = simple_proposal, sd_slab = sd_slab, sd_spike = sd_spike, iter = 1e7, init_model = get_null_model(J), TRUE))
system.time(smart_MH_result_7 <- MR_MH_sampler(J, N, SS, propose_model = smart_proposal, sd_slab = sd_slab, sd_spike = sd_spike, iter = 100000, get_null_model(J), TRUE))

# most_likely_model <- sort(table(MH_result$samples), decreasing = T)[1]

#saveRDS(MH_result, file = paste0('results/MR_Laplace_MH_result', suffix, '.rds'))

model_probs <- sort(prop.table(table(MH_result$samples)), decreasing = T)
models <- names(model_probs)

head(model_probs)

print(paste("There are", sum(model_probs > 0.01), "models with probability greater than one percent."))



# JAMMR -------------------------------------------------------------------

library(R2BGLiMS)

JAMMR_results <- JAMMR(beta_XG, sigma_XG, beta_YG, sigma_YG, nobs_XG, EAF, w = seq(0, 5*nobs_XG, nobs_XG), jam.seed = 2020)
JAMMR_SensitivityPlot(JAMMR_results)

# PolyChordLite -----------------------------------------------------------
source('utils/PolyChord_interface.R')

PC_data_path <- 'dev/PolyChordLite/data/'
SS_file_name <- paste0('MR_model_SSc', suffix, '.txt')
write.table(SS, file = paste0(PC_data_path, '/', SS_file_name), row.names = FALSE, col.names = FALSE)

PC_samples <- list()

for (i in 1:length(models)) {
  if (model_probs[i] > 0.01) {
    
    model_prior_sd <- decode_model(models[i], sd_slab, sd_spike)
    config_filename <- paste0('dev/PolyChordLite/ini/MR_MH_model',  suffix, '_', i, '.ini')
    
    # Write PolyChord configuration file
    write_MR_Gaussian_configuration_file(
      config_file_name = config_filename,
      SS_filename = paste0(PC_data_path, '/', SS_file_name),
      prior_sd = model_prior_sd,
      num_instruments = J,
      num_observations = N,
      slab_precision = 1 / sd_slab / sd_slab,
      spike_precision = 1 / sd_spike / sd_spike,
      PolyChord_control = list(num_live_points = 200)
    )
    
    print(system.time(
      system(paste('dev/PolyChordLite/bin/polychord_MR_Gaussian', config_filename))
    ))
    
    PC_output_file <- paste0('chains/MR_MH_model_', i, '_equal_weights.txt')
    
    posterior_samples <- read_PolyChord_MR_Gaussian_samples(
      paste0('chains/MR_MH_model', suffix, '_', i, '_equal_weights.txt'),
      J, est_sigma_G, model_prior_sd
    )
    
    PC_log_evidence <- read_polychord_logevidence(paste0('chains/MR_MH_model', suffix, '_', i, '.stats'))
    LA_log_evidence <- MR_Laplace_approximation(J, N, SS, random_Gaussian_parameters(J), model_prior_sd)
    
    print(paste("Model:", models[i]))
    print(paste("PolyChord log-evidence:", PC_log_evidence$lev))
    print(paste("Approximated log-evidence:", LA_log_evidence))
    
    PC_samples[[i]] <- list(
      model = models[i],
      PC_log_evidence = PC_log_evidence,
      LA_log_evidence = LA_log_evidence,
      prob = model_probs[i],
      posterior_samples = posterior_samples
    )
  }
}

saveRDS(PC_samples, file = paste0('results/MR_Laplace_PC_samples', suffix, '.rds'))

PC_beta <- do.call('rbind', lapply(PC_samples, function(PC_output) {
  cbind(beta = PC_output$posterior_samples$`beta`, w = PC_output$prob)
}))

plot(density(PC_beta[, 'beta'], weights = PC_beta[, 'w'] / sum(PC_beta[, 'w'])))

# 
# saveRDS(PC_beta, file = paste0('results/MR_Laplace', suffix, '_PC_beta.rds'))



# Analyze multi-chain MAP results -----------------------------------------

MH_results <- readRDS('MR_Laplace_MH_results_bmi_psoriasis.rds')

betas <- sapply(MH_results, '[[', 'betas')
plot(density(betas), main = "Density of causal effect obtained via Laplace approximations")
abline(v = c(0, beta_XY, beta_IVW), col = c('black', 'red', 'blue'))
legend("top", legend = c('zero', 'obs', 'IVW'), lty = 1, col = c('black', 'red', 'blue'))

samples <- sapply(MH_results, '[[', 'samples')
probs <- sort(prop.table(table(samples)), decreasing = TRUE)



# Test causal-confounder hypothesis ---------------------------------------

null_model <- get_null_model(J)

model_1 <- null_model
model_1_LA <- multi_start_MR_Laplace_approximation(J, N, SS, decode_model(model_1, 1, 0.01), num_starting_points = 10)

replicate(10, MR_Laplace_approximation(J, N, SS, random_Gaussian_parameters(J), prior_sd = decode_model(model_7, 1, 0.01)))

model_2 <- model_1
substr(model_2, 2*J+3, 2*J+3) <- "1"
model_2_LA <- multi_start_MR_Laplace_approximation(J, N, SS, decode_model(model_2, 1, 0.01), num_starting_points = 10)

model_3 <- model_1
substr(model_3, 2*J+5, 2*J+5) <- "1"
model_3_LA <- multi_start_MR_Laplace_approximation(J, N, SS, decode_model(model_3, 1, 0.01), num_starting_points = 10)

model_4 <- model_1
substr(model_4, 2*J+7, 2*J+7) <- "1"
model_4_LA <- multi_start_MR_Laplace_approximation(J, N, SS, decode_model(model_4, 1, 0.01), num_starting_points = 10)

model_5 <- model_1
substr(model_5, 2*J+5, 2*J+5) <- "1"
substr(model_5, 2*J+7, 2*J+7) <- "1"
model_5_LA <- multi_start_MR_Laplace_approximation(J, N, SS, decode_model(model_5, 1, 0.01), num_starting_points = 100)

model_7 <- model_1
substr(model_7, 2*J+3, 2*J+3) <- "1"
substr(model_7, 2*J+5, 2*J+5) <- "1"
substr(model_7, 2*J+7, 2*J+7) <- "1"
model_7_LA <- multi_start_MR_Laplace_approximation(J, N, SS, decode_model(model_7, 1, 0.01), num_starting_points = 10)

model_8 <- model_1
substr(model_8, 1, 1) <- "1"
substr(model_8, J+2, J+2) <- "1"
model_8_LA <- multi_start_MR_Laplace_approximation(J, N, SS, decode_model(model_8, 1, 0.01), num_starting_points = 10)

# Fix ---------------------------------------------------------------------

grid_skappa_X <- seq(-0.5, 0.5, 0.05)
grid_skappa_Y <- seq(-0.5, 0.5, 0.05)

grid_cc <- expand.grid(grid_skappa_X, grid_skappa_Y)

post_z_log <- apply(grid_cc, 1, function(cc) {
  
  print(cc)
  
  MAP <- optimize_posterior_log_cc(J, N, SS, SiG, get_ML_solution(SS, cc[1], cc[2], SiG), bad_prior)
  scaled_nl_posterior_log(J, N, SS, SiG, MAP, bad_prior)
})

contour(grid_skappa_X, grid_skappa_Y, matrix(post_z_log, length(grid_skappa_X), length(grid_skappa_Y)), nlevels = 500)
points(LA$optima[[1]]$par[2*J+2], LA$optima[[1]]$par[2*J+3], col = 'red', pch = 3)
points(LA$optima[[2]]$par[2*J+2], LA$optima[[2]]$par[2*J+3], col = 'blue', pch = 3)


