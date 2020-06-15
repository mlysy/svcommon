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

nq <- 3 # number of regular stocks
tau <- runif(nq + 1, -1, 1) # cor(VP, V)
rho <- runif(nq + 1, -1, 1) # cor(V, X)
omega <- runif(nq, -1, 1) # cor(X - V, X0 - V0)

VP <- rnorm(1)
V <- rnorm(nq+1)
X <- rnorm(nq+1)
dt <- runif(1)

ind <- 1:(2*(nq + 1) + 1)
svc_lmvn(VP, V, X, dt, tau, rho, omega)
mvtnorm::dmvnorm(c(VP, V, X)[ind], log = TRUE,
                 sigma = dt * svc_var(tau, rho, omega)[ind,ind,drop=FALSE])

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
    ind <- 1:(2*(nq+1)+1)
    ll <- ll + mvtnorm::dmvnorm(Y[ind], mean = mean_Y[ind], sigma = var_Y[ind,ind,drop=FALSE], log = TRUE)
  }
  ll
}

#' Same thing but `O(nq)`.
svc_loglik2 <- function(alpha, gamma, mu, sigma, tau, rho, omega,
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

svc_loglik(alpha, gamma, mu, sigma, tau, rho, omega,
           log_VPt, log_Vt, Xt, dt)
svc_loglik2(alpha, gamma, mu, sigma, tau, rho, omega,
            log_VPt, log_Vt, Xt, dt)


#--- check TMB code ------------------------------------------------------------

require(TMB)

logit <- function(x, min = 0, max = 1) {
  x <- (x - min)/(max - min)
  log(x) - log(1 - x)
}
ilogit <- function(x, min = 0, max = 1) {
  1/(1 + exp(-x)) * (max - min) + min
}

tmb_flags <- "-std=c++11 -O3 -ffast-math -mtune=native -march=native -Wno-unused-variable -Wno-unused-function -Wno-unknown-pragmas"
model <- "sv_common"
compile(paste0(model, ".cpp"), flags = tmb_flags)
dyn.load(dynlib(model))

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

svc_ad <- MakeADFun(data = list(Xt = Xt, log_VPt = log_VPt, dt = dt),
                    parameters = list(log_Vt = matrix(0, nobs, nq+1),
                                      alpha = rep(0, nq+1),
                                      log_gamma = rep(0, nq+2),
                                      mu = rep(0, nq+2),
                                      log_sigma = rep(0, nq+2),
                                      logit_rho = rep(0, nq+1),
                                      logit_tau = rep(0, nq+1),
                                      logit_omega = rep(0, nq)))

log_gamma <- log(gamma)
log_sigma <- log(sigma)
logit_rho <- logit(rho, -1, 1)
logit_tau <- logit(tau, -1, 1)
logit_omega <- logit(omega, -1, 1)

svc_loglik2(alpha, gamma, mu, sigma, tau, rho, omega,
            log_VPt, log_Vt, Xt, dt)
svc_ad$fn(c(log_Vt, alpha,
            log_gamma, mu, log_sigma,
            logit_rho, logit_tau, logit_omega))

#--- check laplace approximation -----------------------------------------------

require(svcommon)
require(TMB)

nobs <- 1e3 # number of days
nstocks <- 2 # number of stocks, excluding spx

# construct the data
dt <- 1/252
Xt <- as.matrix(cbind(GSPC = snp500[1:nobs, "GSPC"],
                      snp500[1:nobs, 1+1:nstocks]))
Xt <- log(Xt)
log_VPt <- log(as.numeric(snp500[1:nobs, "VIX"]))
# initialize latent variables
log_Vt <- .5 * log(apply(Xt, 2, sv_init, dt = dt, block_size = 10))
# initialize parameters
alpha <- rep(0, nstocks+1)
log_gamma <- rep(0, nstocks+2)
mu <- rep(0, nstocks+2)
log_sigma <- rep(0, nstocks+2)
logit_tau <- rep(0, nstocks+1)
logit_rho <- rep(0, nstocks+1)
logit_omega <- rep(0, nstocks+1)

svc_ad <- MakeADFun(data = list(model = "sv_common",
                                Xt = Xt, log_VPt = log_VPt, dt = dt),
                    parameters = list(log_Vt = log_Vt,
                                      alpha = alpha,
                                      log_gamma = log_gamma,
                                      mu = mu,
                                      log_sigma = log_sigma,
                                      logit_rho = logit_rho,
                                      logit_tau = logit_tau,
                                      logit_omega = logit_omega),
                    random = "log_Vt",
                    DLL = "svcommon_TMBExports", silent = TRUE)

system.time({
  svc_ad$fn(c(alpha, log_gamma, mu, log_sigma, logit_rho, logit_tau, logit_omega))
  svc_ad$gr(c(alpha, log_gamma, mu, log_sigma, logit_rho, logit_tau, logit_omega))
})

# ok let's try optimizing!
system.time({
  opt <- optim(par = svc_ad$par,
               fn = svc_ad$fn,
               gr = svc_ad$gr,
               method = "BFGS", control = list(trace = 1))
})



#--- scratch -------------------------------------------------------------------

sde_sd <- rep(sigma[1], nobs-1)
sde_mean <- log_VPt[1:(nobs-1)] - gamma[1] * (log_VPt[1:(nobs-1)] - mu[1]) * dt
dB_VP <- (log_VPt[2:nobs] - sde_mean)/sde_sd
svc_ad$report(c(log_Vt, alpha,
                log_gamma, mu, log_sigma,
                logit_rho, logit_tau, logit_omega))
list(sde_sd = sde_sd, sde_mean = sde_mean, dB_VP = dB_VP)
