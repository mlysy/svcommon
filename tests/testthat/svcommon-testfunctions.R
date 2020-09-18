#--- a common factor msv model -------------------------------------------------

# correlation structure is:
#
# B_Vi = tau_i B_VP + sqrt(1-tau_i^2) B_epsi
# B_Xi = rho_i B_Vi + sqrt(1-rho_i^2) B_Zi
# B_Zi = omega_i B_Z0 + sqrt(1-omega_i^2) B_etai

#' Full variance matrix.
#'
#' @param tau Vector of length `nq+1` specifying `cor(VP, V)`.
#' @param rho Vector of length `nq+1` specifying `cor(V, X)`.
#' @param omega Vector of length `nq` specifying `cor(X - V, X0 - V0)`
#' @return A correlation matrix of size `2*(nq+1)+1` with the variables ordered as `VP`, `V`, `X`.
svc_var <- function(tau, rho, omega) {
  nq <- length(omega) # number of regular stocks
  Sigma <- diag(2*(nq + 1) + 1)
  Vnm <- paste0("V", 0:nq)
  Xnm <- paste0("X", 0:nq)
  rownames(Sigma) <- c("VP", Vnm, Xnm)
  colnames(Sigma) <- rownames(Sigma)
  Sigma[Vnm, Vnm] <- tau %o% tau + diag(1 - tau^2)
  Sigma[Xnm, Vnm] <- rho * Sigma[Vnm, Vnm]
  Sigma[Vnm, Xnm] <- t(Sigma[Xnm, Vnm])
  Sigma[Xnm, Xnm] <- (rho %o% rho) * (tau %o% tau) +
    sqrt((1-rho^2) %o% (1-rho^2)) * (c(1, omega) %o% c(1, omega))
  Sigma[cbind(Xnm, Xnm)] <- 1
  Sigma["VP", c(Vnm, Xnm)] <- c(tau, rho * tau)
  Sigma[c(Vnm, Xnm), "VP"] <- Sigma["VP", c(Vnm, Xnm)]
  Sigma
}

#' Factorized density evaluation.
#'
svc_lmvn <- function(VP, V, X, dt, tau, rho, omega) {
  nq <- length(omega) # number of regular stocks
  ld_VP <- dnorm(VP, mean = 0, sd = sqrt(dt), log = TRUE)
  ld_V <- sum(dnorm(V, mean = tau * VP,
                    sd = sqrt(dt * (1-tau^2)), log = TRUE))
  ld_X0 <- dnorm(X[1], mean = rho[1] * V[1],
                 sd = sqrt(dt * (1-rho[1]^2)), log = TRUE)
  Z0 <- (X[1] - rho[1] * V[1])/sqrt(1-rho[1]^2)
  ld_X <- sum(dnorm(X[-1],
                    mean = rho[-1] * V[-1] + omega * sqrt(1-rho[-1]^2) * Z0,
                    sd = sqrt(dt * (1-rho[-1]^2) * (1-omega^2)), log = TRUE))
  ld_VP + ld_V + ld_X0 + ld_X
}

#' SVC model loglikelihood.
#'
#' @param alpha Vector of length `nq+1`.
#' @param gamma Vector of length `nq+2`.
#' @param mu Vector of length `nq+2`.
#' @param sigma Vector of length `nq+2`.
#' @param log_VPt Vector of length `nobs`.
#' @param log_Vt Matrix of size `nobs x (nq+1)`.
#' @param Xt Matrix of size `nobs x (nq+1)`.
svc_loglik <- function(alpha, gamma, mu, sigma, tau, rho, omega,
                       log_VPt, log_Vt, Xt, dt) {
  nobs <- length(log_VPt)
  nq <- length(omega)
  Rho <- svc_var(tau, rho, omega) * dt
  ll <- 0
  for(ii in 1:(nobs-1)) {
    log_VPV <- c(log_VPt[ii], log_Vt[ii,])
    mean_VPV <- log_VPV - gamma * (log_VPV - mu) * dt
    sd_VPV <- sigma
    sd_X <- exp(log_VPV[-1])
    mean_X <- Xt[ii,] + (alpha - .5 * sd_X^2) * dt
    Y <- c(log_VPt[ii+1], log_Vt[ii+1,], Xt[ii+1,])
    mean_Y <- c(mean_VPV, mean_X)
    sd_Y <- c(sd_VPV, sd_X)
    var_Y <- t(Rho * sd_Y) * sd_Y
    ll <- ll + lmvn(Y, mu = mean_Y, Sigma = var_Y)
  }
  ll
}

#' Same thing but `O(nq)`.
svc_loglik_fac <- function(alpha, gamma, mu, sigma, tau, rho, omega,
                           log_VPt, log_Vt, Xt, dt) {
  nobs <- length(log_VPt)-1
  nq <- length(omega)
  sqdt <- sqrt(dt)
  mean_VP <- log_VPt[1:nobs] - gamma[1] * (log_VPt[1:nobs] - mu[1]) * dt
  sd_VP <- sigma[1]
  dB_VP <- (log_VPt[1+1:nobs] - mean_VP) / sd_VP
  ll_VP <- sum(dnorm(dB_VP, sd = sqdt, log = TRUE) - log(sd_VP))
  ll_V <- 0
  ll_X <- 0
  for(ii in 1:(nq+1)) {
    mean_V <- log_Vt[1:nobs,ii] - gamma[ii+1] * (log_Vt[1:nobs,ii] - mu[ii+1]) * dt
    sd_V <- sigma[ii+1]
    dB_V <- (log_Vt[1+1:nobs,ii] - mean_V) / sd_V
    ll_V <- ll_V + sum(dnorm(dB_V - tau[ii] * dB_VP, log = TRUE,
                             sd = sqdt * sqrt(1-tau[ii]^2)) - log(sd_V))
    sd_X <- exp(log_Vt[1:nobs,ii])
    mean_X <- Xt[1:nobs,ii] + (alpha[ii] - .5 * sd_X^2) * dt
    if(ii == 1) {
      mean_X <- mean_X + sd_X * (rho[ii] * dB_V)
      sd_X <- sd_X * sqrt(1-rho[ii]^2)
      dB_Z0 <- (Xt[1+1:nobs,ii] - mean_X) / sd_X
      ll_X <- ll_X + sum(dnorm(dB_Z0, sd = sqdt, log = TRUE) - log(sd_X))
    } else {
      mean_X <- mean_X + sd_X * (rho[ii] * dB_V + sqrt(1-rho[ii]^2) * omega[ii-1] * dB_Z0)
      sd_X <- sd_X * sqrt(1-rho[ii]^2) * sqrt(1-omega[ii-1]^2)
      dB_Z <- (Xt[1+1:nobs,ii] - mean_X) / sd_X
      ll_X <- ll_X + sum(dnorm(dB_Z, sd = sqdt, log = TRUE) - log(sd_X))
    }
  }
  ll_VP + ll_V + ll_X
}

#' eOU model likelihood
eou_loglik <- function(alpha, gamma, mu, sigma, rho, log_Vt, Xt, dt) {
  nobs <- length(log_Vt)
  Rho <- cbind(c(1, rho), c(rho, 1)) * dt
  ll <- 0
  for(ii in 1:(nobs-1)) {
    mean_V <- log_Vt[ii] - gamma * (log_Vt[ii] - mu) * dt
    sd_V <- sigma
    sd_X <- exp(log_Vt[ii])
    mean_X <- Xt[ii] + (alpha - .5 * sd_X^2) * dt
    Y <- c(log_Vt[ii+1], Xt[ii+1])
    mean_Y <- c(mean_V, mean_X)
    sd_Y <- c(sd_V, sd_X)
    var_Y <- t(Rho * sd_Y) * sd_Y
    ll <- ll + lmvn(Y, mu = mean_Y, Sigma = var_Y)
  }
  ll
}

#' Log-density of the multivariate normal distribution.
#'
#' @param x Observation vector.
#' @param mu Mean vector.
#' @param Sigma Variance matrix.
#' @return Vector of density evaluations.
lmvn <- function(x, mu, Sigma) {
  if(missing(mu)) mu <- rep(0, length(x))
  C <- chol(Sigma)
  z <- backsolve(C, x-mu, transpose = TRUE)
  -.5 * (sum(z^2) + length(x) * log(2*pi)) - sum(log(diag(C)))
}

# generalized logit and inverse logit
logit <- function(x, min = 0, max = 1) {
  x <- (x - min)/(max - min)
  log(x) - log(1 - x)
}
ilogit <- function(x, min = 0, max = 1) {
  1/(1 + exp(-x)) * (max - min) + min
}

# uniform correlation log-prior on the logit scale
cor_lprior <- function(nu) {
  -nu - 2 * log(1 + exp(-nu))
}

#' Generate a random matrix of size `n x p`.
rmat <- function(n, p) {
  if(missing(p)) p <- n
  matrix(rnorm(n*p), n, p)
}

#' Generate a random array with `dim = ...`.
rarr <- function(...) {
  dim <- unlist(list(...))
  array(rnorm(prod(dim)), dim = dim)
}

#' Remove common factor.
#'
#' @param dB Matrix of size `nobs x nasset`.
#' @param dB0 Vector of length `nobs`.
#' @param logit_rho Vector of length `nasset`.
#' @return Matrix of size `nobs x nasset` corresponding to
#' ```
#' (dB - rho * dB0)/sqrt(1-rho^2)
#' ```.
dB_res <- function(dB, dB0, logit_rho) {
  rho <- svcommon:::rho_itrans(logit_rho)
  rho_sqm <- sqrt(1-rho^2)
  sweep(dB - dB0 %o% rho, 2, rho_sqm, FUN = "/")
}

#' Add common factor.
#'
#' Inverse of [dB_res()].
#'
#' @param dB Matrix of size `nobs x nasset`.
#' @param dB0 Vector of length `nobs`.
#' @param logit_rho Vector of length `nasset`.
#' @return Matrix of size `nobs x nasset` corresponding to
#' ```
#' rho * dB0 + sqrt(1-rho^2) * dB
#' ```.
dB_obs <- function(dB, dB0, logit_rho) {
  rho <- svcommon:::rho_itrans(logit_rho)
  rho_sqm <- sqrt(1-rho^2)
  dB0 %o% rho + sweep(dB, 2, rho_sqm, FUN = "*")
}
