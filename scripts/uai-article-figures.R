

# 0. Setup ----------------------------------------------------------------

library(ggplot2)
library(reshape2)
library(tidyverse)

figures_dir <- "figures/"
if (!dir.exists(figures_dir)) dir.create(figures_dir)

color_blind_palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# 4.2 Computing the Approximation - Figure 3 ------------------------------

load('data/five_optima_posterior_surface.RData')

five_optima_posterior_surface_plot <- ggplot(five_optima_posterior_surface$grid_values) +
  aes(x = x, y = y, z = -z, fill = - five_optima_posterior_surface$simulation_details$num_samples * z) +
  geom_raster(interpolate = TRUE) +
  coord_equal() +
  geom_contour(color = "white", alpha = 0.5, binwidth = 0.005) +
  scale_fill_distiller(palette = "Spectral") +
  ggthemes::theme_tufte() + theme(legend.position="none") +
  xlab(expression(tilde(kappa)[X])) + ylab(expression(tilde(kappa[Y]))) +
  geom_point(aes(x = five_optima_posterior_surface$optima[[1]]$skappa_X, 
                 y = five_optima_posterior_surface$optima[[1]]$skappa_Y), shape = 4) +
  geom_point(aes(x = five_optima_posterior_surface$optima[[2]]$skappa_X, 
                 y = five_optima_posterior_surface$optima[[2]]$skappa_Y), shape = 4) +
  geom_point(aes(x = five_optima_posterior_surface$optima[[3]]$skappa_X, 
                 y = five_optima_posterior_surface$optima[[3]]$skappa_Y), shape = 4) +
  geom_point(aes(x = five_optima_posterior_surface$optima[[4]]$skappa_X, 
                 y = five_optima_posterior_surface$optima[[4]]$skappa_Y), shape = 4) +
  geom_point(aes(x = five_optima_posterior_surface$optima[[5]]$skappa_X, 
                 y = five_optima_posterior_surface$optima[[5]]$skappa_Y), shape = 4)

ggplot2::ggsave(paste0(figures_dir, "Figure_3.pdf"), 
                five_optima_posterior_surface_plot, width = 7, height = 7)

# 5. Empirical Results - Figure 4 -----------------------------------------

load('data/MASSIVE_prior_comparison.RData')

# Prepare data set for plotting by combining beta samples
prior_comparison_betas <- reshape2::melt(
  data.frame(Gaussian = prior_comparison_samples$Gaussian[, 'beta'], 
             Oracle = prior_comparison_samples$Oracle[, 'beta'], 
             MASSIVE = prior_comparison_samples$MASSIVE$betas)
)

prior_comparison_plot <- ggplot(prior_comparison_betas, aes(x = value, fill = variable)) + 
  geom_density(alpha = 0.25) + ggthemes::theme_tufte() + 
  xlim(c(-1.35, -0.85)) + guides(fill = guide_legend(title = "Prior")) +
  theme(legend.position = "top", text = element_text(size = 20)) +
  scale_fill_manual(values = color_blind_palette) +
  xlab(expression(beta)) + ylab("Posterior density") + 
  geom_vline(xintercept = prior_comparison_samples$simulation_details$true_par$beta, linetype = "dashed")

ggplot2::ggsave(paste0(figures_dir, "Figure_4.pdf"), 
                prior_comparison_plot, width = 6, height = 4)


# 5. Empirical Results - Figure 5 -----------------------------------------

load('data/JAMMR_MASSIVE_comparison.RData')

JAMMR_MASSIVE_comparison_plot <- ggplot(
  JAMMR_MASSIVE_comparison,
  aes(x = valid, y = rmse, group = Algorithm, color = Algorithm)) +
  geom_line(size = 1) +
  facet_grid(rows = vars(`N:sigma`), cols = vars(beta), scales = "free_y", 
             labeller = labeller(beta = label_both)) +
  ggthemes::theme_tufte() + 
  theme(strip.text.x = element_text(size = 18), strip.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        axis.title = element_text(size = 20, face = "bold"), legend.position = "top") +
  xlab("Number of valid instruments") + ylab("") + 
  scale_x_continuous(breaks = 0:10) +
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin = rmse_lower_ci, ymax = rmse_upper_ci),
                position = position_dodge(width = 0.1))

ggsave(paste0(figures_dir, "Figure_5.pdf"), 
       JAMMR_MASSIVE_comparison_plot, width = 6, height = 6)


# 6.1 Determinants of Macroeconomic Growth - Figure 6 ---------------------

load('data/determinants_macroeconomic_growth.RData')


plot_determinants_macroeconomic_growth_beta <- function(data_frame) {
  ggplot(data_frame, aes(x = value, fill = variable)) + 
    geom_density(alpha = 0.25) + xlim(c(-1, 3.5)) +
    ggthemes::theme_tufte() + guides(fill=guide_legend(title = "Algorithm")) +
    theme(legend.position = "top", text = element_text(size=18)) + 
    scale_fill_manual(values = color_blind_palette) +
    xlab(expression(beta)) + ylab("Posterior density")
}

determinants_macroeconomic_growth_plot <- ggpubr::ggarrange(
  plot_determinants_macroeconomic_growth_beta(
    data_frame = reshape2::melt(data.frame(
      ivbma = institutions_macroeconomic_growth_samples$ivbma[, 'Rule'], 
      MASSIVE = institutions_macroeconomic_growth_samples$MASSIVE$betas
    ))
  ),
  plot_determinants_macroeconomic_growth_beta(
    data_frame = reshape2::melt(data.frame(
      ivbma = integration_macroeconomic_growth_samples$ivbma[, 'Trade'], 
      MASSIVE = integration_macroeconomic_growth_samples$MASSIVE$betas
    ))
  )
)

ggsave(paste0(figures_dir, "Figure_6.pdf"), 
       determinants_macroeconomic_growth_plot, width = 10, height = 4)


# 6.2 Investigating the Relationship between BMI and Psoriasis - Figure 7 ----

load('data/BMI_psoriasis.RData')

BMI_psoriasis_plot <- ggplot2::ggplot(
    data.frame(beta = BMI_psoriasis_samples$betas), 
    aes(x = beta, fill = color_blind_palette[3])
  ) + geom_density(alpha=0.25) + xlim(c(-0.15, 0.15)) + 
  ggplot2::geom_vline(xintercept = log(c(1.08, 1.04, 1.13)), linetype = c(1, 2, 2)) +
  ggthemes::theme_tufte() + xlab(expression(beta)) +  ylab("Posterior density") +
  scale_fill_manual(values = color_blind_palette[3]) + theme(legend.position = "none")

ggsave(paste0(figures_dir, "Figure_7.pdf"), BMI_psoriasis_plot, width = 5, height = 4)


