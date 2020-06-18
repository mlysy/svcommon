#--- test with real data -------------------------------------------------------

require(svcommon)

nobs <- 1e3 # number of days
nstocks <- 10 # number of stocks, excluding spx

# construct the data
dt <- 1/252
Xt <- as.matrix(cbind(GSPC = snp500[1:nobs, "GSPC"],
                      snp500[1:nobs, 1+1:nstocks]))
Xt <- log(Xt)
log_VPt <- log(as.numeric(snp500[1:nobs, "VIX"]))
# initialize latent variables
log_Vt <- log(apply(Xt, 2, sv_init, dt = dt, block_size = 10))
# initialize parameters
alpha <- rep(0, nstocks+1)
log_gamma <- rep(0, nstocks+2)
mu <- rep(0, nstocks+2)
log_sigma <- rep(0, nstocks+2)
logit_rho <- rep(0, nstocks+1)
logit_tau <- rep(0, nstocks+1)
logit_omega <- rep(0, nstocks)

curr_par <- list(log_Vt = log_Vt,
                 alpha = alpha,
                 log_gamma = log_gamma,
                 mu = mu,
                 log_sigma = log_sigma,
                 logit_rho = logit_rho,
                 logit_tau = logit_tau,
                 logit_omega = logit_omega)

#' @param stock Which stock to optimize over.
#' @param univ If `TRUE` applies only to parameters in univariate model.
svc_map <- function(stock, univ = FALSE) {
  map <- list(alpha = 0:nstocks,
              log_gamma = -1:nstocks,
              mu = -1:nstocks,
              log_sigma = -1:nstocks,
              logit_rho = 0:nstocks,
              logit_tau = 0:nstocks,
              logit_omega = 1:nstocks)
  # all NA except `stock`
  lapply(map, function(x) {
    x[x != stock] <- NA
    factor(x, levels = if(all(is.na(x))) NULL else stock)
  })
}

#' @param stock Which parameter to update
#' @param new_par New parameter vector.
#' @param old_par Old parameter list.
svc_update <- function(stock, univ = FALSE, new_par, old_par) {
  map <- svc_map(stock, univ)
  for(nm in names(new_par)) {
    old_par[[nm]][!is.na(map[[nm]])] <- new_par[nm]
  }
  old_par
}

# initialize parameters via individual stocks
for(istock in 0:nstocks) {
  message("parameter = ", istock)
  eou_ad <- eou_MakeADFun(Xt = Xt[,istock+1], dt = dt,
                          alpha = alpha[istock+1],
                          log_Vt = log_Vt[,istock+1],
                          log_gamma = log_gamma[istock+2],
                          mu = mu[istock+2],
                          log_sigma = log_sigma[istock+2],
                          logit_rho = logit_rho[istock+1])
  ## eou_ad2 <- MakeADFun(data = list(model = "sv_eou",
  ##                                 Xt = Xt[,istock+1], dt = dt),
  ##                     parameters = list(log_Vt = log_Vt[,istock+1],
  ##                                       alpha = alpha[istock+1],
  ##                                       log_gamma = log_gamma[istock+2],
  ##                                       mu = mu[istock+2],
  ##                                       log_sigma = log_sigma[istock+2],
  ##                                       logit_rho = logit_rho[istock+1]),
  ##                     random = "log_Vt",
  ##                     DLL = "svcommon_TMBExports", silent = TRUE)
  tm <- system.time({
    opt <- optim(par = eou_ad$par,
                 fn = eou_ad$fn,
                 gr = eou_ad$gr,
                 method = "BFGS", control = list(trace = 1))
  })
  message("Time: ", round(tm[3], 2), " seconds")
  curr_par <- svc_update(istock, univ = TRUE,
                         new_par = opt$par, old_par = curr_par)
  curr_par$log_Vt[,istock+1] <- eou_ad$env$last.par.best[1:nobs]
}

# blockwise coordinate descent

for(istock in -1:nstocks) {
  message("parameter = ", istock)
  svc_ad <- MakeADFun(data = list(model = "sv_common",
                                  Xt = Xt, log_VPt = log_VPt, dt = dt),
                      parameters = curr_par,
                      random = "log_Vt",
                      map = svc_map(istock),
                      DLL = "svcommon_TMBExports", silent = TRUE)
  tm <- system.time({
    opt <- optim(par = svc_ad$par,
                 fn = svc_ad$fn,
                 gr = svc_ad$gr,
                 method = "BFGS", control = list(trace = 1))
  })
  message("Time: ", round(tm[3], 2), " seconds")
  # update parameters
  curr_par <- svc_update(istock, new_par = opt$par, old_par = curr_par)
  curr_par$log_Vt <- matrix(svc_ad$env$last.par.best[1:(nobs * (nstocks+1))],
                            nobs, nstocks+1)
}


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
