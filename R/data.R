#' A dataset containing results for reproducing Figure 3 in Section 4 of the UAI article.
#'
#' @format A list containing all the necessary elements for reproducing the 
#' posterior contour surface in Figure 3.
#' \describe{
#'   \item{optima}{List containing the parameters found at each of the five
#'   local optima in list format.}
#'   \item{grid_values}{Data frame containing maximum posterior values (z) 
#'   computed on a two-dimensional grid (x, y) of values for the scale-free
#'   confounding coefficients (skappa_X, skappa_Y).}
#'   \item{simulation_details}{Details of the data and generating model in
#'   the simulation.}
#' }
"five_optima_posterior_surface"


#' A dataset containing results for reproducing Figure 4 in Section 5 of the UAI article.
#'
#' @format A list containing the posterior samples obtained with the three
#' priors under comparison together with the data simulation details.
#' \describe{
#'   \item{Gaussian}{Posterior samples obtained with rstan when using Gaussian 
#'   prior for the MASSIVE model parameters.}
#'   \item{Oracle}{Posterior samples obtained with rstan when using the Oracle
#'   spike-and-slab prior for the MASSIVE model parameters.}
#'   \item{MASSIVE}{Posterior samples output by MASSIVE algorithm.}
#'   \item{simulation_details}{Details of the data and generating model in
#'   the simulation.}
#' }
"prior_comparison_samples"


#' A dataset containing results for reproducing Figure 5 in Section 5 of the UAI article.
#' 
#' @format A data frame containing evaluations of JAM-MR and MASSIVE in various
#' simulated scenarios.
#' \describe{
#'   \item{Algorithm}{Name of algorithm, either JAM-MR or MASSIVE.}
#'   \item{valid}{Integer number of valid instruments in the simulation.}
#'   \item{beta}{Numeric value of beta (causal effect from X to Y) in the simulation.}
#'   \item{rmse}{Root-mean-square error (RMSE) achieved by the algorithm.}
#'   \item{rmse_se}{Standard error for RMSE metric across replications.}
#'   \item{rmse_lower_ci}{Bootstrapped confidence interval 2.5% quantile for RMSE across replications.}
#'   \item{rmse_upper_ci}{Bootstrapped confidence interval 97.5% quantile for RMSE across replications.}
#'   \item{N:sigma}{Simulation configuration (number of observations):(noise magnitude).}
#' }
"JAMMR_MASSIVE_comparison"


#' A dataset containing results for reproducing Figure 6 in Section 6 of the UAI article.
#' 
#' @format List containing posterior samples for the causal effect of rule-of-law
#' (institutions) on growth, as generated by two competing algorithms, MASSIVE and ivbma.
"institutions_macroeconomic_growth_samples"

#' A dataset containing results for reproducing Figure 6 in Section 6 of the UAI article.
#' 
#' @format List containing posterior samples for the causal effect of (economic)
#' integration on growth, as generated by two competing algorithms, MASSIVE and ivbma.
"integration_macroeconomic_growth_samples"


#' A dataset containing results for reproducing Figure 7 in Section 6 of the UAI article.
"BMI_psoriasis_samples"
#'
#' @format The data set follows the output of the MASSIVE algorithm.
#' \describe{
#'   \item{beta}{Massive posterior samples for the causal effect of the exposure
#'   on the outcome.}
#'   \item{originating_model}{Approximate IV model from which the posterior
#'   sample was generated.}
#' }