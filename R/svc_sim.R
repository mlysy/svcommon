#' Simulate.
#'
#' @param nobs Length of time series.
#' @template param_dt
#' @param X0 Vector of length `nasset + 1` or matrix of size `(nasset + 1) x nseries` of asset log prices at time `t = 0`.
#' @param log_VP0 Scalar or vector of length `nseries` of volatility proxy values at time `t = 0` on the log standard deviation scale.
#' @param log_V0 Vector of length `nasset + 1` or matrix of size `(nasset + 1) x nseries` of volatilities at time `t = 0` on the log standard deviation scale.
#' @param alpha Vector of length `(nasset + 1)` or matrix of size `(nasset + 1) x nseries` of asset growth rate parameters.
#' @param log_gamma Vector of length `nasset + 2` or matrix of size `(nasset + 2) x nseries` log-volatility mean reversion parameters on the log scale.  The first two correspond to the volatility proxy and the common-factor asset's volatility, respectively.
#' @param mu Vector of length `nasset + 2` or matrix of size `(nasset + 2) x nseries` of log-volatility mean parameters.
#' @param log_sigma Vector of length `nasset + 2` or matrix of size `(nasset + 2) x nseries` or log-volatility diffusion parameters on the log scale.
#' @param logit_rho Vector of length `nasset + 1` or matrix of size `(nasset + 1) x nseries` correlation parameters between asset and volatility innovations, on the logit scale.  The first one is that of the common-factor asset proxy.  See 'Details'.
#' @param logit_tau Vector of length `nasset + 1` or matrix of size `(nasset + 1) x nseries` correlation parameters between the latent volatilities and the volatility proxy.
#' @param logit_omega Vector of length `nasset` or matrix of size `nasset x nseries` correlation parameters between the residual asset price of the common-factor proxy and the other residual asset prices.
#' @return A list containing arrays `Xt` and `log_Vt` of size `nobs x (nasset + 1) x nseries`, and the matrix `log_VPt` of size `nobs x nseries` containing the simulations of `nseries` SVC processes observed at times `t = dt, 2dt, ..., nobs*dt`.
svc_sim <- function(nobs, dt, X0, log_VP0, log_V0,
                    alpha, log_gamma, mu, log_sigma,
                    logit_rho, logit_tau, logit_omega) {
  # format inputs
  # transpose matrices so that each asset is a column
  X0 <- t(check_matrix(X0, promote = TRUE))
  log_V0 <- t(check_matrix(log_V0, promote = TRUE))
  alpha <- t(check_matrix(alpha, promote = TRUE))
  log_gamma <- t(check_matrix(log_gamma, promote = TRUE))
  mu <- t(check_matrix(log_gamma, promote = TRUE))
  log_sigma <- t(check_matrix(log_sigma, promote = TRUE))
  logit_rho <- t(check_matrix(logit_rho, promote = TRUE))
  logit_tau <- t(check_matrix(logit_tau, promote = TRUE))
  logit_omega <- t(check_matrix(logit_omega, promote = TRUE))
  # allocate memory
  nseries <- get_max(X0, log_VP0, log_V0, alpha, log_gamma, mu, log_sigma,
                     logit_rho, logit_tau, logit_omega)
  nasset <- ncol(X0) - 1
  Xt <- array(NA, dim = c(nobs, nasset+1, nseries))
  log_VPt <- matrix(NA, nobs, nasset+1)
  log_Vt <- array(NA, dim = c(nobs, nasset+1, nseries))
  # parameter transformations
  gamma <- exp(gamma)
  sigma <- exp(sigma)
  rho <- rho_itrans(logit_rho)
  rho_sqm <- sqrt(1 - rho^2)
  tau <- rho_itrans(logit_tau)
  tau_sqm <- sqrt(1 - tau^2)
  omega <- rho_itrans(logit_omega)
  omega_sqm <- sqrt(1- omega^2)
  # initialize recurrence
  Xcurr <- X0
  log_Vcurr <- log_V0
  log_VPcurr <- log_VP0
  for(ii in 1:nobs) {
    # volatility proxy
    dB_VP <- sqrt(dt) * rnorm(nseries)
    mean_VP <- log_VPcurr - gamma[,1] * (log_VPcurr - mu[,1]) * dt
    sd_VP <- sigma[,1]
    log_VPcurr <- mean_VP + sd_VP * dB_VP
    # latent volatilites
    dB_V <- sqrt(dt) * matrix(rnorm((nasset+1)*nseries), nseries, nasset+1)
    dB_V <- tau * dB_VP + tau_sqm * dB_V
    mean_V <- log_Vcurr - gamma[,-1] * (log_Vcurr - mu[,-1]) * dt
    sd_V <- sigma[,-1]
    log_Vcurr <- mean_V + sd_V * dB_V
    # assets
    sd_X <- exp(log_Vcurr)
    mean_X <- Xcurr + (alpha - .5 * sd_X^2) * dt
    dB_Z <- sqrt(dt) * matrix(rnorm((nasset+1)*nseries), nseries, nasset+1)
    # correlation between asset innovations and common factor
    dB_Z[,-1] <- omega * dB_Z[,1] + omega_sqm * dB_Z[,-1]
    Xcurr <- mean_X + sd_X * (rho * dB_V + rho_sqm * dB_Z)
    # storage
    Xt[ii,,] <- Xcurr
    log_VPt[ii,] <- log_VPcurr
    log_Vt[ii,,] <- log_Vcurr
  }
  list(Xt = Xt, log_Vt = log_Vt, log_VPt = log_VPt)
}
