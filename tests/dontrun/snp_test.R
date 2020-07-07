#--- test with real data -------------------------------------------------------

require(svcommon)

nobs <- 1000 # number of days
nasset <- 5 # number of assets, excluding spx

# construct the data
dt <- 1/252
Xt <- as.matrix(cbind(GSPC = snp500[1:nobs, "GSPC"],
                      snp500[1:nobs, 1+1:nasset]))
Xt <- log(Xt)
log_VPt <- log(as.numeric(snp500[1:nobs, "VIX"]))
# initialize latent variables
log_Vt <- log(apply(Xt, 2, sv_init, dt = dt, block_size = 10))
# initialize parameters
alpha <- rep(0, nasset+1)
log_gamma <- rep(0, nasset+2)
mu <- rep(0, nasset+2)
log_sigma <- rep(0, nasset+2)
logit_rho <- rep(0, nasset+1)
logit_tau <- rep(0, nasset+1)
logit_omega <- rep(0, nasset)

curr_par <- list(log_Vt = log_Vt,
                 alpha = alpha,
                 log_gamma = log_gamma,
                 mu = mu,
                 log_sigma = log_sigma,
                 logit_rho = logit_rho,
                 logit_tau = logit_tau,
                 logit_omega = logit_omega)



## #' @param iasset Which parameter to update
## #' @param new_par New parameter vector.
## #' @param old_par Old parameter list.
## svc_update2 <- function(iasset, nasset, new_par, old_par) {
##   map <- svcommon:::svc_map(iasset, nasset)
##   for(nm in names(new_par)) {
##     old_par[[nm]][!is.na(map[[nm]])] <- new_par[nm]
##   }
##   old_par
## }

# initialize parameters of log_Vt

#' Estimate parameters of an Ornstein-Uhlenbeck process.
#'
#' Uses least-squares on the Euler approximation.
ou_fit <- function(Xt, dt) {
  nobs <- length(Xt)
  M <- lm(diff(Xt) ~ Xt[-nobs])
  beta <- as.numeric(coef(M))
  gamma <- -beta[2]/dt
  mu <- beta[1]/(gamma*dt)
  sigma <- sqrt(mean(resid(M)^2)/dt)
  c(gamma = gamma, mu = mu, sigma = sigma)
}

VPt_par <- ou_fit(log_VPt, dt)
curr_par$log_gamma[1] <- log(VPt_par["gamma"])
curr_par$mu[1] <- VPt_par["mu"]
curr_par$log_sigma[1] <- log(VPt_par["sigma"])

# initialize parameters via individual asset
for(iasset in 0:nasset) {
  message("asset = ", iasset)
  eou_ad <- eou_MakeADFun(Xt = Xt[,iasset+1], dt = dt,
                          alpha = curr_par$alpha[iasset+1],
                          log_Vt = curr_par$log_Vt[,iasset+1],
                          log_gamma = curr_par$log_gamma[iasset+2],
                          mu = curr_par$mu[iasset+2],
                          log_sigma = curr_par$log_sigma[iasset+2],
                          logit_rho = curr_par$logit_rho[iasset+1])
  ## eou_ad2 <- MakeADFun(data = list(model = "sv_eou",
  ##                                 Xt = Xt[,iasset+1], dt = dt),
  ##                     parameters = list(log_Vt = log_Vt[,iasset+1],
  ##                                       alpha = alpha[iasset+1],
  ##                                       log_gamma = log_gamma[iasset+2],
  ##                                       mu = mu[iasset+2],
  ##                                       log_sigma = log_sigma[iasset+2],
  ##                                       logit_rho = logit_rho[iasset+1]),
  ##                     random = "log_Vt",
  ##                     DLL = "svcommon_TMBExports", silent = TRUE)
  ## all.equal(eou_ad, eou_ad2)
  ## tm <- system.time({
  ##   opt <- optim(par = eou_ad$par,
  ##                fn = eou_ad$fn,
  ##                gr = eou_ad$gr,
  ##                method = "BFGS", control = list(trace = 1))
  ## })
  ## tm <- system.time({
  ##   opt <- nlminb(start = eou_ad$par,
  ##                 objective = eou_ad$fn,
  ##                 gradient = eou_ad$gr, control = list(trace = 10))
  ## })
  tm <- system.time({
    opt <- nlm(p = eou_ad$par,
               f = function(par) {
                 out <- eou_ad$fn(par)
                 attr(out,"gradient") <- eou_ad$gr()
                 out
               }, print.level = 1)
  })
  message("Time: ", round(tm[3], 2), " seconds")
  curr_par <- svc_update(eou_ad, old_par = curr_par, iasset = iasset)
  ## curr_par2 <- svc_update2(iasset, nasset = nasset,
  ##                         new_par = opt$par, old_par = curr_par)
  ## curr_par2$log_Vt[,iasset+1] <- eou_ad$env$last.par.best[1:nobs]
  ## message("svc_update == svc_update2: ", identical(curr_par, curr_par2))
}

# initialize the remaining parameters

#' Extract the Browian increments from the SVC model.
#'
#' @return A list with elements:
#' \describe{
#'   \item{`V`}{An `nobs x (nasset+2)` matrix of the log-volatility innovations.}
#'   \item{`X`}{An `nobs x (nasset+1)` matrix of log-asset innovations.}
#'   \item{`Z`}{An `nobs x (naxsset+1)` matrix of residual log-asset innovations, after account for the log-volatility innovations.  That is,
#'     ```
#'     dB_Z = (dB_X - rho dB_V) / sqrt(1 - rho^2).
#'     ```
#'   }
#' }
svc_dB <- function(Xt, log_VPt, log_Vt, dt,
                   alpha, log_gamma, mu, log_sigma, logit_rho) {
  nobs <- nrow(Xt)
  nasset <- ncol(Xt)-1
  ind0 <- 1:(nobs-1)
  # do all the volatilities at once
  dB_V <- cbind(log_VPt, log_Vt)
  dB_V <- t(apply(dB_V, 2, diff)) + exp(log_gamma) * (t(dB_V[ind0,]) - mu) * dt
  dB_V <- dB_V / exp(log_sigma)
  # innovations for nassets and asset common factor
  Vt <- exp(t(log_Vt[ind0,]))
  dB_X <- t(apply(Xt, 2, diff)) - (alpha - .5 * Vt^2) * dt
  dB_X <- dB_X / Vt
  # residual part after removing dB_V
  rho <- 2/(1 + exp(-logit_rho)) - 1
  dB_Z <- (dB_X - rho * dB_V[-1,]) / sqrt(1 - rho^2)
  list(V = t(dB_V), X = t(dB_X), Z = t(dB_Z))
}

#' Convert correlation parameter to unconstrained scale
cor_trans <- function(rho) {
  nu <- .5 * (rho + 1)
  log(nu) - log(1-nu)
}

dB <- svc_dB(Xt, log_VPt = log_VPt, log_Vt = curr_par$log_Vt, dt = dt,
             alpha = curr_par$alpha, log_gamma = curr_par$log_gamma,
             mu = curr_par$mu, log_sigma = curr_par$log_sigma,
             logit_rho = curr_par$logit_rho)

curr_par$logit_tau[] <- as.numeric(cor_trans(cor(dB$V[,1], dB$V[,-1])))
curr_par$logit_omega[] <- as.numeric(cor_trans(cor(dB$Z[,1], dB$Z[,-1])))

# blockwise coordinate descent
# checked: fix_Vt = T/F gives same result, former much faster :)

curr_par2 <- curr_par
curr_par3 <- curr_par

for(iasset in -1:nasset) {
  message("asset = ", iasset)
  svc_ad <- svc_MakeADFun(Xt = Xt, log_VPt = log_VPt, dt = dt,
                          par_list = curr_par,
                          iasset = iasset, fix_Vt = TRUE)
  ## svc_ad2 <- MakeADFun(data = list(model = "sv_common",
  ##                                  Xt = Xt, log_VPt = log_VPt, dt = dt),
  ##                      parameters = curr_par,
  ##                      random = "log_Vt",
  ##                      map = svc_map(iasset, nasset = nasset),
  ##                      DLL = "svcommon_TMBExports", silent = TRUE)
  ## all.equal(svc_ad, svc_ad2)
  ## tm <- system.time({
  ##   opt <- optim(par = svc_ad$par,
  ##                fn = svc_ad$fn,
  ##                gr = svc_ad$gr,
  ##                method = "BFGS", control = list(trace = 1))
  ## })
  ## tm <- system.time({
  ##   opt <- nlminb(start = svc_ad$par,
  ##                 objective = svc_ad$fn,
  ##                 gradient = svc_ad$gr, control = list(trace = 10))
  ## })
  tm <- system.time({
    opt <- nlm(p = svc_ad$par,
               f = function(par) {
                 out <- svc_ad$fn(par)
                 attr(out,"gradient") <- svc_ad$gr()
                 out
               }, print.level = 1)
  })
  message("Time: ", round(tm[3], 2), " seconds")
  # update parameters
  curr_par <- svc_update(svc_ad, old_par = curr_par, iasset = iasset)
  ## curr_par2 <- svc_update(svc_ad, old_par = curr_par2, iasset = iasset)
  ## curr_par3 <- svc_update(svc_ad, old_par = curr_par3, iasset = iasset)
  ## curr_par2 <- svc_update2(iasset, nasset = nasset,
  ##                         new_par = opt$par, old_par = curr_par)
  ## curr_par2$log_Vt <- matrix(svc_ad$env$last.par.best[1:(nobs * (nasset+1))],
  ##                           nobs, nasset+1)
  ## message("svc_update == svc_update2: ", identical(curr_par, curr_par2))
}

# joint opt
iasset <- "all"
svc_ad <- svc_MakeADFun(Xt = Xt, log_VPt = log_VPt, dt = dt,
                        par_list = curr_par,
                        iasset = iasset)

tm <- system.time({
  opt <- optim(par = svc_ad$par,
               fn = svc_ad$fn,
               gr = svc_ad$gr,
               method = "BFGS", control = list(trace = 1))
})

curr_par <- svc_update(svc_ad, old_par = curr_par, iasset = iasset)

# check mode
finfo <- numDeriv::hessian(svc_ad$fn, x = svc_ad$par)

optimCheck::optim_proj(xsol = svc_ad$par,
                       fun = svc_ad$fn,
                       maximize = FALSE)

#--- scratch -------------------------------------------------------------------

svc_update(svc_ad, curr_par, iasset = iasset)

svc_update(2, opt$par, old_par = curr_par)

svc_ad <- MakeADFun(data = list(model = "sv_common",
                                Xt = Xt, log_VPt = log_VPt, dt = dt),
                    parameters = curr_par,
                    random = "log_Vt",
                    map = svc_map(2),
                    DLL = "svcommon_TMBExports", silent = TRUE)

system.time({
  svc_ad$fn(c(alpha, log_gamma[-1], mu, log_sigma, logit_rho, logit_tau, logit_omega))
  svc_ad$gr(c(alpha, log_gamma, mu, log_sigma, logit_rho, logit_tau, logit_omega))
})

# ok let's try optimizing!
system.time({
  opt <- optim(par = svc_ad$par,
               fn = svc_ad$fn,
               gr = svc_ad$gr,
               method = "BFGS", control = list(trace = 1))
})



#--- convert to package functions ----------------------------------------------

#' Need:
#'
#' - Wrapped `ADFun` constructors for the models.  This for argument checking, default inputs, and mappings.  Can also input parameters as a par_list.
#' - `svc_update()`: update par_list.
