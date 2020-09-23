#' Simulate time series from the exponential Ornstein-Uhlenbeck stochastic volatility model.
#'
#' @param nobs Length of time series.
#' @template param_dt
#' @param X0 Scalar or vector of `nseries` asset log prices at time `t = 0`.
#' @param log_V0 Scalar or vector of `nseries` volatilities at time `t = 0`, on the log standard deviation scale.
#' @param alpha Scalar or vector of `nseries` growth rate parameters.
#' @param log_gamma Scalar or vector of `nseries` log-volatility mean reversion parameters on the log scale.
#' @param mu Scalar or vector of `nseries` log-volatility mean parameters.
#' @param log_sigma Scalar or vector of `nseries` log-volatility diffusion parameters on the log scale.
#' @param logit_rho Scalar or vector of `nseries` correlation parameters between asset and volatility innovations, on the logit scale.
#' @param dBt An optional list with elements `V` and `Z` corresponding to matrices of size `nobs x nseries` of pre-specified Brownian innovations.  If missing these consist of iid draws from an `N(0, dt)` distribution.
#' @return A list containing matrices `Xt` and `log_Vt` of `nobs x nseries` of eOU observations, where each column corresponds to a process observed at times `t = dt, 2dt, ..., nobs*dt`.
#' @export
eou_sim <- function(nobs, dt, X0, log_V0,
                    alpha, log_gamma, mu, log_sigma, logit_rho, dBt) {
  # allocate memory
  nseries <- get_max(X0, log_V0, alpha, log_gamma, mu, log_sigma, logit_rho,
                     tvec = TRUE)
  Xt <- matrix(NA, nobs, nseries)
  log_Vt <- matrix(NA, nobs, nseries)
  # parameter transformations
  gamma <- exp(log_gamma)
  sigma <- exp(log_sigma)
  rho <- rho_itrans(logit_rho)
  rho_sqm <- sqrt(1 - rho^2)
  # initialize recurrence
  Xcurr <- X0
  log_Vcurr <- log_V0
  for(ii in 1:nobs) {
    if(!missing(dB)) {
      dB_V <- dB$V[ii,]
      dB_Z <- dB$Z[ii,]
    } else {
      dB_V <- sqrt(dt) * rnorm(nseries)
      dB_Z <- sqrt(dt) * rnorm(nseries)
    }
    mean_V <- log_Vcurr - gamma * (log_Vcurr - mu) * dt
    sd_V <- sigma
    sd_X <- exp(log_Vcurr)
    mean_X <- Xcurr + (alpha - .5 * sd_X^2) * dt
    log_Vcurr <- mean_V + sd_V * dB_V
    Xcurr <- mean_X + sd_X * (rho * dB_V + rho_sqm * dB_Z)
    # store
    Xt[ii,] <- Xcurr
    log_Vt[ii,] <- log_Vcurr
  }
  list(Xt = Xt, log_Vt = log_Vt)
}

#--- helper functions ----------------------------------------------------------

#' Returns the maximum size of the inputs, where "size" is the trailing dimension of arrays.
#'
#' @param tvec Logical; whether vectors should be transposed into one-row matrices (`tvec = TRUE`) or promoted to one-column matrices (`tvec = FALSE`).
#' @noRd
get_max <- function(..., tvec) {
  N <- sapply(list(...), function(x) {
    if(!is.array(x)) {
      x <- if(tvec) t(x) else as.matrix(x)
    }
    n <- dim(x)[length(dim(x))]
  })
  max(N)
}

#' Convert logit_rho to rho.
#' @noRd
rho_itrans <- function(logit_rho) {
  2/(1 + exp(-logit_rho)) - 1
}
