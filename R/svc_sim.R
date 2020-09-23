#' Simulate time series from the SVC stochastic volatility model.
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
#' @param dBt An optional list of Brownian innovations:
#' \describe{
#'   \item{`VP`}{An `nobs x nseries` matrix.}
#'   \item{`V`}{An `nobs x (nasset+1) x nseries` array.}
#'   \item{`Z`}{An `nobs x (nasset+1) x nseries` array.}
#' }
#' If missing these are all iid draws from an `N(0, dt)` distribution.
#' @return A list containing arrays `Xt` and `log_Vt` of size `nobs x (nasset + 1) x nseries`, and the matrix `log_VPt` of size `nobs x nseries` containing the simulations of `nseries` SVC processes observed at times `t = dt, 2dt, ..., nobs*dt`.
#' @export
svc_sim <- function(nobs, dt, X0, log_VP0, log_V0,
                    alpha, log_gamma, mu, log_sigma,
                    logit_rho, logit_tau, logit_omega, dBt) {
  # problem dimensions
  nseries <- get_max(X0, t(log_VP0), log_V0, alpha, log_gamma, mu, log_sigma,
                     logit_rho, logit_tau, logit_omega, tvec = FALSE)
  nasset <- if(is.vector(X0)) length(X0) - 1 else nrow(X0) - 1
  # format inputs
  # transpose matrices so that each asset is a column
  X0 <- t(check_matrix(X0 = X0, dim = c(nasset+1, nseries), promote = TRUE))
  log_VP0 <- check_vector(log_VP0 = log_VP0, len = nseries, promote = TRUE)
  log_V0 <- t(check_matrix(log_V0 = log_V0,
                           dim = c(nasset+1, nseries), promote = TRUE))
  alpha <- t(check_matrix(alpha = alpha, dim = c(nasset+1, nseries),
                          promote = TRUE))
  log_gamma <- t(check_matrix(log_gamma = log_gamma,
                              dim = c(nasset+2, nseries), promote = TRUE))
  mu <- t(check_matrix(mu = mu,
                       dim = c(nasset+2, nseries), promote = TRUE))
  log_sigma <- t(check_matrix(log_sigma = log_sigma,
                              dim = c(nasset+2, nseries), promote = TRUE))
  logit_rho <- t(check_matrix(logit_rho = logit_rho,
                              dim = c(nasset+1, nseries), promote = TRUE))
  logit_tau <- t(check_matrix(logit_tau = logit_tau,
                              dim = c(nasset+1, nseries), promote = TRUE))
  logit_omega <- t(check_matrix(logit_omega = logit_omega,
                                dim = c(nasset, nseries), promote = TRUE))
  # allocate memory
  # storage is consistent with assets as last dimension.
  # will aperm Xt and log_Vt later
  Xt <- array(NA, dim = c(nobs, nseries, nasset+1))
  log_VPt <- matrix(NA, nobs, nseries)
  log_Vt <- array(NA, dim = c(nobs, nseries, nasset+1))
  # parameter transformations
  gamma <- exp(log_gamma)
  sigma <- exp(log_sigma)
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
    if(!missing(dBt)) {
      dB_VP <- dBt$VP[ii,]
    } else {
      dB_VP <- sqrt(dt) * rnorm(nseries)
    }
    mean_VP <- log_VPcurr - gamma[,1] * (log_VPcurr - mu[,1]) * dt
    sd_VP <- sigma[,1]
    # latent volatilites
    if(!missing(dBt)) {
      dB_V <- t(dBt$V[ii,,])
      dB_Z <- t(dBt$Z[ii,,])
    } else {
      dB_V <- sqrt(dt) * matrix(rnorm((nasset+1)*nseries), nseries, nasset+1)
      dB_Z <- sqrt(dt) * matrix(rnorm((nasset+1)*nseries), nseries, nasset+1)
    }
    dB_V <- tau * dB_VP + tau_sqm * dB_V
    mean_V <- log_Vcurr - gamma[,-1] * (log_Vcurr - mu[,-1]) * dt
    sd_V <- sigma[,-1]
    # assets
    sd_X <- exp(log_Vcurr)
    mean_X <- Xcurr + (alpha - .5 * sd_X^2) * dt
    # correlation between asset innovations and common factor
    dB_Z[,-1] <- omega * dB_Z[,1] + omega_sqm * dB_Z[,-1]
    # recurrence update
    log_VPcurr <- mean_VP + sd_VP * dB_VP
    log_Vcurr <- mean_V + sd_V * dB_V
    Xcurr <- mean_X + sd_X * (rho * dB_V + rho_sqm * dB_Z)
    # storage
    log_VPt[ii,] <- log_VPcurr
    log_Vt[ii,,] <- log_Vcurr
    Xt[ii,,] <- Xcurr
  }
  list(Xt = aperm(Xt, c(1,3,2)),
       log_Vt = aperm(log_Vt, c(1,3,2)),
       log_VPt = log_VPt)
}
