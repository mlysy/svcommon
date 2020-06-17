
context("eou_adfun")

source("svcommon-testfunctions.R")

test_that("R and TMB likelihoods agree.", {
  ntest <- 10
  for(ii in 1:ntest) {
    nobs <- sample(2:20,1) # number of observations
    alpha <- rnorm(1)
    gamma <- rexp(1)
    mu <- rnorm(1)
    sigma <- rexp(1)
    rho <- runif(1, -1, 1)
    # parameter transformations
    log_gamma <- log(gamma)
    log_sigma <- log(sigma)
    logit_rho <- logit(rho, -1, 1)
    # data
    log_Vt <- rnorm(nobs)
    Xt <- rnorm(nobs)
    dt <- runif(1)
    # TMB object
    eou_ad <- TMB::MakeADFun(data = list(model = "sv_eou",
                                         Xt = Xt, dt = dt),
                             parameters = list(log_Vt = rep(0, nobs),
                                               alpha = 0,
                                               log_gamma = 0,
                                               mu = 0,
                                               log_sigma = 0,
                                               logit_rho = 0),
                             DLL = "svcommon_TMBExports", silent = TRUE)
    expect_equal(
      eou_loglik(alpha = alpha, gamma = gamma, mu = mu, sigma = sigma,
                 rho = rho,
                 log_Vt = log_Vt, Xt = Xt, dt = dt),
      -eou_ad$fn(c(log_Vt, alpha,
                   log_gamma, mu, log_sigma, logit_rho))
    )
  }
})

