#--- test with real data -------------------------------------------------------

require(svcommon)

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

#' Convert correlation parameter to unconstrained scale.
cor_trans <- function(rho) {
  nu <- .5 * (rho + 1)
  log(nu) - log(1-nu)
}

#--- data setup ----------------------------------------------------------------

nobs <- 2000 # number of days
nasset <- 10 # number of assets, excluding spx

# construct the data
dt <- 1/252
Xnames <- colnames(snp500)
Xnames <- Xnames[!Xnames %in% c("Date", "GSPC", "VIX")]
Xnames <- sample(Xnames, nasset)
## Xnames <- colnames(snp500)[1+1:nasset]
Xt <- as.matrix(cbind(GSPC = snp500[1:nobs, "GSPC"],
                      snp500[1:nobs, Xnames]))
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


#--- initialize parameters via marginal fits -----------------------------------

# initialize parameters of log_Vt
VPt_par <- ou_fit(log_VPt, dt)
curr_par$log_gamma[1] <- log(VPt_par["gamma"])
curr_par$mu[1] <- VPt_par["mu"]
curr_par$log_sigma[1] <- log(VPt_par["sigma"])

opt_method <- "nlminb"

system.time({
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
    if(opt_method == "optim") {
      tm <- system.time({
        opt <- optim(par = eou_ad$par,
                     fn = eou_ad$fn,
                     gr = eou_ad$gr,
                     method = "BFGS", control = list(trace = 1))
      })
    } else if(opt_method == "nlminb") {
      tm <- system.time({
        opt <- nlminb(start = eou_ad$par,
                      objective = eou_ad$fn,
                      gradient = eou_ad$gr, control = list(trace = 10))
      })
    } else if(opt_method == "nlm") {
      tm <- system.time({
        opt <- nlm(p = eou_ad$par,
                   f = function(par) {
                     out <- eou_ad$fn(par)
                     attr(out,"gradient") <- eou_ad$gr()
                     out
                   }, print.level = 1)
      })
    } else stop("Invalid opt_method.")
    message("Time: ", round(tm[3], 2), " seconds")
    curr_par <- svc_update(eou_ad, old_par = curr_par, iasset = iasset)
    ## curr_par2 <- svc_update2(iasset, nasset = nasset,
    ##                         new_par = opt$par, old_par = curr_par)
    ## curr_par2$log_Vt[,iasset+1] <- eou_ad$env$last.par.best[1:nobs]
    ## message("svc_update == svc_update2: ", identical(curr_par, curr_par2))
  }
})

# initialize the correlation parameters
dB <- svc_dB(Xt, log_VPt = log_VPt, log_Vt = curr_par$log_Vt, dt = dt,
             alpha = curr_par$alpha, log_gamma = curr_par$log_gamma,
             mu = curr_par$mu, log_sigma = curr_par$log_sigma,
             logit_rho = curr_par$logit_rho)
curr_par$logit_tau[] <- as.numeric(cor_trans(cor(dB$V[,1], dB$V[,-1])))
curr_par$logit_omega[] <- as.numeric(cor_trans(cor(dB$Z[,1], dB$Z[,-1])))

#--- blockwise coordinate descent ----------------------------------------------

# checked: fix_Vt = T/F gives same result, former much faster :)

## curr_par2 <- curr_par
## curr_par3 <- curr_par

nepoch <- 3

system.time({
  for(iepoch in 1:nepoch) {
    message("--- epoch: ", iepoch, " ---")
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
    if(opt_method == "optim") {
      tm <- system.time({
        opt <- optim(par = svc_ad$par,
                     fn = svc_ad$fn,
                     gr = svc_ad$gr,
                     method = "BFGS", control = list(trace = 1))
      })
    } else if(opt_method == "nlminb") {
      tm <- system.time({
        opt <- nlminb(start = svc_ad$par,
                      objective = svc_ad$fn,
                      gradient = svc_ad$gr, control = list(trace = 10))
      })
    } else if(opt_method == "nlm") {
      tm <- system.time({
        opt <- nlm(p = svc_ad$par,
                   f = function(par) {
                     out <- svc_ad$fn(par)
                     attr(out,"gradient") <- svc_ad$gr()
                     out
                   }, print.level = 1)
      })
    } else stop("Invalid opt_method.")
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
  }
})

#--- joint optimization --------------------------------------------------------

iasset <- "all"
svc_ad <- svc_MakeADFun(Xt = Xt, log_VPt = log_VPt, dt = dt,
                        par_list = curr_par,
                        iasset = iasset, fix_Vt = FALSE)

system.time({
  if(opt_method == "optim") {
    tm <- system.time({
      opt <- optim(par = svc_ad$par,
                   fn = svc_ad$fn,
                   gr = svc_ad$gr,
                   method = "BFGS", control = list(trace = 1))
    })
  } else if(opt_method == "nlminb") {
    tm <- system.time({
      opt <- nlminb(start = svc_ad$par,
                    objective = svc_ad$fn,
                    gradient = svc_ad$gr, control = list(trace = 10))
    })
  } else if(opt_method == "nlm") {
    tm <- system.time({
      opt <- nlm(p = svc_ad$par,
                 f = function(par) {
                   out <- svc_ad$fn(par)
                   attr(out,"gradient") <- svc_ad$gr()
                   out
                 }, print.level = 1)
    })
  } else stop("Invalid opt_method.")
})

## tm <- system.time({
##   opt <- optim(par = svc_ad$par,
##                fn = svc_ad$fn,
##                gr = svc_ad$gr,
##                method = "BFGS", control = list(trace = 1))
## })

curr_par <- svc_update(svc_ad, old_par = curr_par, iasset = iasset)

system.time({
  est <- sdreport(svc_ad)
})

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
