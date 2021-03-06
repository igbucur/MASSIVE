% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Laplace_approximation.R
\name{derive_smart_starting_points}
\alias{derive_smart_starting_points}
\title{Routine for computing smart starting points from the sufficient statistics.}
\usage{
derive_smart_starting_points(SS, sigma_G = binomial_sigma_G(SS))
}
\arguments{
\item{SS}{Numeric matrix containing first- and second-order statistics.}

\item{sigma_G}{Numeric vector of instrument standard deviations.}
}
\value{
List containing three smart starting points on the ML manifold for
the posterior optimization.
}
\description{
Routine for computing smart starting points from the sufficient statistics.
}
\examples{
J <- 5 # number of instruments
N <- 1000 # number of samples
parameters <- random_Gaussian_parameters(J) 
EAF <- runif(J, 0.1, 0.9) # EAF random values
dat <- generate_data_MASSIVE_model(N, 2, EAF, parameters)
derive_smart_starting_points(dat$SS)
}
