
# 0. Setup and Data Preprocessing -----------------------------------------

library(ivbma)
data(growth)

N <- length(growth$Y)
J <- ncol(growth$Z) + ncol(growth$W) # remove intercept
Jr <- ncol(growth$Z) + 1

Intercept <- growth$W[, 1] + rnorm(nrow(growth$W), 0, 1e-3)

data_1 <- as.matrix(cbind(1, growth$Z, Intercept, growth$W[, -1], growth$X[1], growth$Y))
data_2 <- as.matrix(cbind(1, growth$Z, Intercept, growth$W[, -1], growth$X[2], growth$Y))


  
SS_1 <- t(data_1) %*% data_1 / N
sigma_G_1 <- apply(data_1[, 2:(J+1)], 2, sd)
hyp_1 <- determine_hyperparameters(J, N, SS_1, sigma_G_1)
result_1 <- MASSIVE(J, N, SS_1, sigma_G_1, hyp_1$sd_slab, hyp_1$sd_spike, max_iter = 5000,
                    greedy_search = parallel_greedy_search, Laplace_approximation = robust_safe_smart_LA_log)
saveRDS(result_1, file = 'results/final_robust_MASSIVE_growth_institutions.rds')

oh_1 <- old_determine_hyperparameters(J, N, SS_1, sigma_G_1)
result_oh_1 <- MASSIVE(J, N, SS_1, sigma_G_1, oh_1$sd_slab, oh_1$sd_spike, max_iter = 1000,
                    greedy_search = parallel_greedy_search, Laplace_approximation = safe_smart_LA_log)
saveRDS(result_oh_1, file = 'results/initial_MASSIVE_growth_institutions.rds')

# 
SS_2 <- t(data_2) %*% data_2 / N
sigma_G_2 <- apply(data_2[, 2:(J+1)], 2, sd)
hyp_2 <- determine_hyperparameters(J, N, SS_2, sigma_G_2)
result_2 <- MASSIVE(J, N, SS_2, sigma_G_2, hyp_2$sd_slab, hyp_2$sd_spike, max_iter = 5000,
                    greedy_search = parallel_greedy_search, Laplace_approximation = robust_safe_smart_LA_log)
saveRDS(result_2, file = 'results/final_robust_MASSIVE_growth_integration_500.rds')

oh_2 <- old_determine_hyperparameters(J, N, SS_2, sigma_G_2)
result_oh_2 <- MASSIVE(J, N, SS_2, sigma_G_2, oh_2$sd_slab, oh_2$sd_spike, max_iter = 1000,
                       greedy_search = parallel_greedy_search, Laplace_approximation = safe_smart_LA_log)
saveRDS(result_oh_2, file = 'results/initial_MASSIVE_growth_integration.rds')


# result_1_alt <- MASSIVE(J, N, SS_1, sigma_G_1, sd_slab = 0.1, sd_spike = 0.001)
# incl_prob_1 <- analyze_MASSIVE_output(result_1_alt, J, N)
# 
# result_1_sl0_sp2 <- MASSIVE(J, N, SS_1, sigma_G_1, sd_slab = 1.0, sd_spike = 0.01)
# incl_prob_2 <- analyze_MASSIVE_output(result_1_sl0_sp2, J, N)

# result_1_fin <- MASSIVE(J, N, SS_1, sigma_G_1, sd_slab = 10, sd_spike = 0.001)
# incl_prob_1 <- analyze_MASSIVE_output(result_1_alt, J, N)

data_r1 <- as.matrix(cbind(1, growth$Z, Intercept, growth$X[1], Growth = growth$Y))
data_r2 <- as.matrix(cbind(1, growth$Z, Intercept, growth$X[2], Growth = growth$Y))

SS_r1 <- t(data_r1) %*% data_r1 / N
sigma_G_r1 <- apply(data_r1[, 2:(Jr+1)], 2, sd)
hyp_r1 <- determine_hyperparameters(Jr, N, SS_r1, sigma_G_r1)
result_r1 <- MASSIVE(Jr, N, SS_r1, sigma_G_r1, hyp_r1$sd_slab, hyp_r1$sd_spike, max_iter = 5000)
analyze_MASSIVE_output(result_r1, Jr, N)

saveRDS(result_r1, file = 'results/final_MASSIVE_growth_institutions_Z.rds')

SS_r2 <- t(data_r2) %*% data_r2 / N
sigma_G_r2 <- apply(data_r2[, 2:(Jr+1)], 2, sd)
result_r2 <- MASSIVE(Jr, N, SS_r2, sigma_G_r2, sd_slab = 1, sd_spike = 0.01, max_iter = 5000)
analyze_MASSIVE_output(result_r2, Jr, N)

saveRDS(result_r2, file = 'results/MASSIVE_growth_institutions_dh_Z.rds')

# This is possibly a hint that when you have few covariates, we really should set the hyperparameters instead of determining them empirically.

ivbma_output <- ivbma(growth$Y, growth$X, growth$Z, growth$W, s = 5000)
ivbma_noz <- ivbma(growth$Y, growth$X, matrix(0, 54, 0), growth$W)
summary(ivbma_output, nms.U = c(names(growth$Z), names(growth$W)),nms.V = c(names(growth$X), names(growth$W)))

df <- as.data.frame(data_1)

lm_XG <- df %>% 
  dplyr::select(-`Rule`, -`1`, -`growth$Y`) %>% 
  purrr::map(~lm(df$Rule ~ .x)) %>% 
  # purrr::map(summary) %>%
  purrr::map(c("coefficients"))
beta_XG <- lm_XG %>% purrr::map_dbl(2)
sigma_XG <- lm_XG %>% purrr::map_dbl(4)

lm_YG <- df %>% 
  dplyr::select(-`Rule`, -`1`, -`growth$Y`) %>% 
  purrr::map(~lm(df$`growth$Y` ~ .x)) %>% 
  purrr::map(summary) %>%
  purrr::map(c("coefficients"))
beta_YG <- lm_YG %>% purrr::map_dbl(2)
sigma_YG <- lm_YG %>% purrr::map_dbl(4)

incl_prob_1 <- analyze_model_approximations(result_1)
sel <- which(incl_prob_1$gammas[2,] > 0.5)
beta_XG <- beta_XG[sel]
beta_YG <- beta_YG[sel]
sigma_XG <- sigma_XG[sel]
sigma_YG <- sigma_YG[sel]

MendelianRandomization::mr_ivw(MendelianRandomization::mr_input(beta_XG, sigma_XG, beta_YG, sigma_YG))

lm(Rule ~ ., data = as.data.frame(data_1) %>% select(-c(`1`, `growth$Y`)))
build_

run_JAMMR(list(Z = data_1))


# Plot results for growth data --------------------------------------------

ivbma_growth <- ivbma::ivbma(growth$Y, growth$X, growth$Z, growth$W, s = 20000, run.diagnostics = TRUE)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

institutions <- readRDS('results/final_robust_MASSIVE_growth_institutions_1000.rds')
institutions_models <- compare_and_select_models(institutions, 3e-3)
institutions_samp <- get_LA_posterior_samples(J, N, institutions_models, nrow(ivbma_growth$rho))

df <- data.frame(ivbma = ivbma_growth$rho[,1], MASSIVE = institutions_samp$betas)
library(ggplot2);library(reshape2)
institutions_plot <- ggplot(melt(df),aes(x=value, fill=variable)) + geom_density(alpha=0.25) +
  ggthemes::theme_tufte() + xlim(c(-1, 3.5)) + guides(fill=guide_legend(title="Algorithm")) +
  theme(legend.position = "top") +   scale_fill_manual(values = cbPalette) + theme(text = element_text(size=18)) +
  xlab(expression(beta)) + ylab("Posterior density")
# ggplot2::ggsave("../tex/MASSIVE/plot_beta_institutions_growth.png", width = 5, height = 4)
ggplot2::ggsave("../tex/MASSIVE/Camera-Ready Version/plot_beta_institutions_growth.pdf", institutions_plot, width = 5, height = 4)



#ICML_institutions_growth_plot <- bayesplot::mcmc_areas(data.frame(beta = samples$betas), pars = "beta") + 
 #ggplot2::geom_vline(xintercept = log(c(1.08, 1.04, 1.13)), linetype = c(1, 2, 2)) +

# integration <- readRDS('results/MASSIVE_growth_integration_dh_intercept.rds')
integration <- readRDS('results/final_robust_MASSIVE_growth_integration_1000.rds')
integration_models <- compare_and_select_models(integration, 3e-3)
integration_samp <- get_LA_posterior_samples(J, N, integration_models, nrow(ivbma_growth$rho))

df <- data.frame(ivbma = ivbma_growth$rho[,2], MASSIVE = integration_samp$betas)
library(ggplot2);library(reshape2)
integration_plot <- ggplot(melt(df), aes(x=value, fill=variable)) + geom_density(alpha=0.25) +
  ggthemes::theme_tufte() + xlim(c(-1, 3.5)) + guides(fill=guide_legend(title="Algorithm")) +
  theme(legend.position = "top") +   scale_fill_manual(values = cbPalette) + theme(text = element_text(size=18)) +
  xlab(expression(beta)) + ylab("Posterior density")
  #ggtitle("Estimated effect of economic integration on macroeconomic growth")
ggplot2::ggsave("../tex/MASSIVE/Camera-Ready Version/plot_beta_integration_growth.pdf", integration_plot, width = 5, height = 4)

# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
growth_plot <- ggpubr::ggarrange(institutions_plot, integration_plot, common.legend = TRUE)

ggplot2::ggsave("../tex/MASSIVE/Camera-Ready Version/plot_determinants_growth.pdf", growth_plot, width = 5, height = 3)


# Margarine trial ---------------------------------------------------------

library(ivbma)
data(margarine)

N <- length(margarine$Y)
J <- ncol(margarine$Z) + ncol(margarine$W)
Jr <- ncol(margarine$Z) + 1

Intercept <- margarine$W[, 3] + rnorm(nrow(margarine$W), 0, 1e-3)

data_m <- as.matrix(cbind(1, margarine$Z, Intercept, margarine$W[, -3], margarine$X[1], margarine$Y))
SS_m <- t(data_m) %*% data_m / N
sigma_G_m <- apply(data_m[, 2:(J+1)], 2, sd)
result_m <- MASSIVE(J, N, SS_m, sigma_G_m)

data_rm <- as.matrix(cbind(1, margarine$Z, Intercept, margarine$X[1], margarine$Y))
SS_rm <- t(data_rm) %*% data_rm / N
sigma_G_rm <- apply(data_rm[, 2:(Jr+1)], 2, sd)
result_rm <- MASSIVE(Jr, N, SS_rm, sigma_G_rm)

# saveRDS(result_m, "results/MASSIVE_margarine_price_dh_intercept.rds")
# saveRDS(result_rm, "results/MASSIVE_margarine_price_dh_intercept_Z.rds")

samp <- get_LA_posterior_samples(Jr, N, result_rm, 100000)

ivbma_margarine <- ivbma::ivbma(margarine$Y, margarine$X, margarine$Z, margarine$W, s = 10000, run.diagnostics = TRUE)

price <- readRDS('results/MASSIVE_margarine_price_dh_intercept_new.rds')
price_samp <- get_LA_posterior_samples(J, N, price, nrow(ivbma_margarine$rho))

df <- data.frame(ivbma = ivbma_margarine$rho[, 1], MASSIVE = price_samp$betas)
library(ggplot2);library(reshape2)
ggplot(melt(df), aes(x=value, fill=variable)) + geom_density(alpha=0.25) +
  ggthemes::theme_tufte() + xlim(c(-5, 5)) + guides(fill=guide_legend(title="Algorithm")) +
  theme(legend.position = "top") +   scale_fill_manual(values = cbPalette) + theme(text = element_text(size=20)) +
  xlab(expression(beta)) + ylab("Posterior density")
#ggtitle("Estimated effect of economic integration on macroeconomic growth")
ggplot2::ggsave("../tex/MASSIVE/plot_beta_integration_growth.png", width = 5, height = 4)
