% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MASSIVE.R
\name{MASSIVE}
\alias{MASSIVE}
\title{MASSIVE algorithm (Model Assessment and Stochastic Search for Instrumental Variable Estimation)}
\usage{
MASSIVE(
  J,
  N,
  SS,
  sigma_G,
  sd_slab,
  sd_spike,
  max_iter = 1000,
  greedy_search = parallel_greedy_search,
  Laplace_approximation = safe_Laplace_approximation,
  pruning_threshold = 0.003,
  posterior_samples = 10000,
  ...
)
}
\arguments{
\item{J}{Integer number of candidate instruments.}

\item{N}{Integer number of observations.}

\item{SS}{Numeric matrix containing first- and second-order statistics.}

\item{sigma_G}{Numeric vector of instrument standard deviations.}

\item{sd_slab}{Numeric scale parameter of slab component.}

\item{sd_spike}{Numeric scale parameter of spike component.}

\item{max_iter}{Maximum number of stochastic search steps.}

\item{greedy_search}{Function for initial greedy search.}

\item{Laplace_approximation}{Function used to compute Laplace approximation of IV model.}

\item{pruning_threshold}{Numeric threshold for pruning approximated IV models.
Models with probability less that threshold are pruned out.}

\item{posterior_samples}{Integer number of posterior samples to generate.}

\item{...}{Extra arguments for greedy_search, Laplace_approximations and find_causal_models.}
}
\value{
List of explored IV model approximations and their evidences
}
\description{
MASSIVE algorithm (Model Assessment and Stochastic Search for Instrumental Variable Estimation)
}
\examples{
set.seed(2020)
J <- 10
N <- 10000
G <- matrix(rbinom(N * J, 2, 0.3), N, J)
U <- rnorm(N)
X <- G \%*\% runif(J, 0.3, 0.5) + U + rnorm(N)
Y <- G[, 1:5] \%*\% runif(5, 0.1, 0.3) + X + U + rnorm(N)

Z <- cbind(1, G, X, Y)
SS <- t(Z) \%*\% Z / N
sigma_G <- apply(G, 2, sd)

samples <- MASSIVE::MASSIVE(J, N, SS, sigma_G, sd_slab = 1, sd_spike = 0.01, max_iter = 1000)
plot(density(samples$betas))
median(samples$betas)
}
