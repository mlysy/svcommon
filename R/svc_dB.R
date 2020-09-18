#' Extract the Brownian increments from the SVC model.
#'
#' @template param_Xt
#' @template param_log_VPt
#' @template param_log_Vt
#' @template param_dt
#' @template param_alpha
#' @template param_log_gamma
#' @template param_mu
#' @template param_log_sigma
#' @template param_logit_rho
#' @return A list with elements:
#' \describe{
#'   \item{`V`}{An `nobs x (nasset+2)` matrix of the log-volatility innovations.}
#'   \item{`X`}{An `nobs x (nasset+1)` matrix of log-asset innovations.}
#'   \item{`Z`}{An `nobs x (nasset+1)` matrix of residual log-asset innovations, after account for the log-volatility innovations.  That is,
#'     ```
#'     dB_Z = (dB_X - rho dB_V) / sqrt(1 - rho^2).
#'     ```
#'   }
#' }
#' @export
svc_dB <- function(Xt, log_VPt, log_Vt, dt,
                   alpha, log_gamma, mu, log_sigma, logit_rho) {
  nobs <- nrow(Xt)
  nasset <- ncol(Xt)-1
  ind0 <- 1:(nobs-1)
  # do all the volatilities at once
  dB_V <- cbind(log_VPt, log_Vt)
  dimnames(dB_V) <- NULL
  dB_V <- t(matdiff(dB_V)) +
    exp(log_gamma) * (t(dB_V[ind0,,drop=FALSE]) - mu) * dt
  dB_V <- dB_V / exp(log_sigma)
  # innovations for nassets and asset common factor
  Vt <- exp(t(log_Vt[ind0,,drop=FALSE]))
  dB_X <- t(matdiff(Xt)) - (alpha - .5 * Vt^2) * dt
  dB_X <- dB_X / Vt
  # residual part after removing dB_V
  rho <- 2/(1 + exp(-logit_rho)) - 1
  dB_Z <- (dB_X - rho * dB_V[-1,,drop=FALSE]) / sqrt(1 - rho^2)
  list(V = t(dB_V), X = t(dB_X), Z = t(dB_Z))
}

#--- helper functions ----------------------------------------------------------

#' Apply `diff` to each column of a matrix.
#'
#' Output is always a matrix.
#' @noRd
matdiff <- function(X) {
  nobs <- nrow(X)
  dX <- apply(X, 2, diff)
  if(nobs == 2) dX <- t(dX)
  dX
}
