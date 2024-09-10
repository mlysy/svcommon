
context("svc_adfun")

## source("svcommon-testfunctions.R")

test_that("Variance matrix factorization is correct.", {
  ntest <- 10
  for(ii in 1:ntest) {
    nq <- 3 # number of regular stocks
    tau <- runif(nq + 1, -1, 1) # cor(VP, V)
    rho <- runif(nq + 1, -1, 1) # cor(V, X)
    omega <- runif(nq, -1, 1) # cor(X - V, X0 - V0)
    VP <- rnorm(1)
    V <- rnorm(nq+1)
    X <- rnorm(nq+1)
    dt <- runif(1)
    expect_equal(svc_lmvn(VP = VP, V = V, X = X, dt = dt,
                          tau = tau, rho = rho, omega = omega),
                 lmvn(c(VP, V, X), Sigma = dt * svc_var(tau, rho, omega)))
  }
})

test_that("Likelihood factorization is correct.", {
  ntest <- 15
  for(ii in 1:ntest) {
    nq <- sample(5, 1) # number of regular stocks
    nobs <- sample(2:20,1) # number of observations
    alpha <- rnorm(nq+1)
    gamma <- rexp(nq+2)
    mu <- rnorm(nq+2)
    sigma <- rexp(nq+2)
    tau <- runif(nq + 1, -1, 1) # cor(VP, V)
    rho <- runif(nq + 1, -1, 1) # cor(V, X)
    omega <- runif(nq, -1, 1) # cor(X - V, X0 - V0)
    log_VPt <- rnorm(nobs)
    log_Vt <- matrix(rnorm(nobs*(nq+1)), nobs, nq+1)
    Xt <- matrix(rnorm(nobs*(nq+1)), nobs, nq+1)
    dt <- runif(1)
    expect_equal(
      svc_loglik(alpha = alpha, gamma = gamma, mu = mu, sigma = sigma,
                 tau = tau, rho = rho, omega = omega,
                 log_VPt = log_VPt, log_Vt = log_Vt, Xt = Xt, dt = dt),
      svc_loglik_fac(alpha = alpha, gamma = gamma, mu = mu, sigma = sigma,
                     tau = tau, rho = rho, omega = omega,
                     log_VPt = log_VPt, log_Vt = log_Vt, Xt = Xt, dt = dt)
    )
  }
})

test_that("R and TMB likelihoods agree.", {
  ntest <- 10
  for(ii in 1:ntest) {
    nq <- sample(5, 1) # number of regular stocks
    nobs <- sample(2:20,1) # number of observations
    alpha <- rnorm(nq+1)
    gamma <- rexp(nq+2)
    mu <- rnorm(nq+2)
    sigma <- rexp(nq+2)
    tau <- runif(nq + 1, -1, 1) # cor(VP, V)
    rho <- runif(nq + 1, -1, 1) # cor(V, X)
    omega <- runif(nq, -1, 1) # cor(X - V, X0 - V0)
    # parameter transformations
    log_gamma <- log(gamma)
    log_sigma <- log(sigma)
    logit_rho <- logit(rho, -1, 1)
    logit_tau <- logit(tau, -1, 1)
    logit_omega <- logit(omega, -1, 1)
    # data
    log_VPt <- rnorm(nobs)
    log_Vt <- matrix(rnorm(nobs*(nq+1)), nobs, nq+1)
    Xt <- matrix(rnorm(nobs*(nq+1)), nobs, nq+1)
    dt <- runif(1)
    # TMB object
    svc_ad <- TMB::MakeADFun(data = list(model = "sv_common",
                                         Xt = Xt, log_VPt = log_VPt, dt = dt),
                             parameters = list(log_Vt = matrix(0, nobs, nq+1),
                                               alpha = rep(0, nq+1),
                                               log_gamma = rep(0, nq+2),
                                               mu = rep(0, nq+2),
                                               log_sigma = rep(0, nq+2),
                                               logit_rho = rep(0, nq+1),
                                               logit_tau = rep(0, nq+1),
                                               logit_omega = rep(0, nq)),
                             DLL = "svcommon_TMBExports", silent = TRUE)
    expect_equal(
      svc_loglik(alpha = alpha, gamma = gamma, mu = mu, sigma = sigma,
                 tau = tau, rho = rho, omega = omega,
                 log_VPt = log_VPt, log_Vt = log_Vt, Xt = Xt, dt = dt) +
      sum(cor_lprior(logit_rho)) + sum(cor_lprior(logit_tau)) +
      sum(cor_lprior(logit_omega)),
      -svc_ad$fn(c(log_Vt, alpha,
                   log_gamma, mu, log_sigma,
                   logit_rho, logit_tau, logit_omega))
    )
  }
})
