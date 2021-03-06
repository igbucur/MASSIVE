% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posterior.R
\name{scaled_neg_log_gradient}
\alias{scaled_neg_log_gradient}
\title{Routine for computing the negative log-posterior gradient for the MASSIVE model}
\usage{
scaled_neg_log_gradient(
  J,
  N,
  SS,
  sigma_G,
  param_list,
  prior_sd,
  n = 2,
  log_scale = TRUE
)
}
\arguments{
\item{J}{Integer number of genetic instrumental variables.}

\item{N}{Integer number of observations.}

\item{SS}{Numeric matrix containing first- and second-order statistics.}

\item{sigma_G}{Numeric vector of genetic IV standard deviations.}

\item{param_list}{List of IV model parameter values.}

\item{prior_sd}{List of standard deviations for the parameter Gaussian priors.}

\item{n}{Integer number of alleles (trials) for the binomial genetic variants.}

\item{log_scale}{Logical flag indicating whether scale parameters (sigma_X, sigma_Y)
are given on the log-scale or on the normal scale.}
}
\value{
numeric vector; negative log-posterior gradient for the MASSIVE model.
}
\description{
Routine for computing the negative log-posterior gradient for the MASSIVE model
}
