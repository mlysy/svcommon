% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/svcommon-package.R
\docType{package}
\name{svcommon-package}
\alias{svcommon}
\alias{svcommon-package}
\title{svcommon: Fast Inference for Common-Factor Stochastic Volatility Models}
\description{
Provides various tools for estimating the parameters of the common-factor multivariate stochastic volatility model of Fang et al (2020) \doi{10.1002/cjs.11536} and extensions. In particular, the complete-data likelihood implementation scales linearly in the number of assets, and latent volatilities are efficiently marginalized using the Laplace approximation in the R package 'TMB' with very high accuracy. Combined with a carefully initialized block coordinate descent algorithm, maximum likelihood estimation can be conducted two orders of magnitude faster than with alternative parameter inference algorithms.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/mlysy/svcommon}
  \item \url{https://mlysy.github.io/svcommon/}
  \item Report bugs at \url{https://github.com/mlysy/svcommon/issues}
}

}
\author{
\strong{Maintainer}: Martin Lysy \email{mlysy@uwaterloo.ca}

}
\keyword{internal}
