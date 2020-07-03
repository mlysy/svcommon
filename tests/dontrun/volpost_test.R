#--- getting the volatility posterior distribution -----------------------------

require(svcommon)

nobs <- 250 # number of days
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

iasset <- 0 # SPX
eou_ad <- eou_MakeADFun(Xt = Xt[,iasset+1], dt = dt,
                        alpha = curr_par$alpha[iasset+1],
                        log_Vt = curr_par$log_Vt[,iasset+1],
                        log_gamma = curr_par$log_gamma[iasset+2],
                        mu = curr_par$mu[iasset+2],
                        log_sigma = curr_par$log_sigma[iasset+2],
                        logit_rho = curr_par$logit_rho[iasset+1])

tm <- system.time({
  opt <- optim(par = eou_ad$par,
               fn = eou_ad$fn,
               gr = eou_ad$gr,
               method = "BFGS", control = list(trace = 1))
})

# ok let's sdreport
sdreport(eou_ad)
sqrt(diag(sdreport(eou_ad)$cov))
sdreport(eou_ad)$value

# ok let's try optimization with part of latent variables fixed
