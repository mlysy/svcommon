# svcommon: Fast Inference for Common-Factor Stochastic Volatility Models

*Martin Lysy*

---

### Description

Provides various tools for estimating the parameters of the common-factor multivariate stochastic volatility model of Fang et al (2020) <doi:10.1002/cjs.11536> and extensions.  In particular, the complete-data likelihood implementation scales linearly in the number of assets, and latent volatilities are efficiently marginalized using the Laplace approximation in the R package **TMB** with very high accuracy.  Combined with a carefully initialized block coordinate descent algorithm, maximum likelihood estimation can be conducted two orders of magnitude faster than with alternative parameter inference algorithms.

### Installation

Install the R package [**devtools**](https://CRAN.R-project.org/package=devtools) and run
```r
devtools::install_github("mlysy/svcommon")
```

### Usage

Please see package vignette: `vignette("svcommon")`.
