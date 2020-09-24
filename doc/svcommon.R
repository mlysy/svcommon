params <-
list(load_calcs = TRUE)

## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load_data, include = FALSE, message = FALSE------------------------------
# required packages
library(svcommon)
library(mvtnorm) # multivariate normal draws
library(tidyr); library(dplyr) # data frame manipulation
library(ggplot2) # plotting

dim(snp500) # automatically loaded with svcommon
snp500[1:6, ncol(snp500)-8 + 1:8]

## ----disp_data, ref.label = "load_data"---------------------------------------
# required packages
library(svcommon)
library(mvtnorm) # multivariate normal draws
library(tidyr); library(dplyr) # data frame manipulation
library(ggplot2) # plotting

dim(snp500) # automatically loaded with svcommon
snp500[1:6, ncol(snp500)-8 + 1:8]

## ----init_setup, include = FALSE----------------------------------------------
# problem dimensions
nobs <- 1000 # number of days
nasset <- 5 # number of assets, excluding GSPC

# format the data
dt <- 1/252
Xt <- as.matrix(cbind(GSPC = snp500[1:nobs, "GSPC"],
                      snp500[1:nobs, 1+1:nasset]))
Xt <- log(Xt)
log_VPt <- log(as.numeric(snp500[1:nobs, "VIX"]))

# initialize latent variables with rolling window standard deviations
log_Vt <- log(apply(Xt, 2, sv_init, dt = dt, block_size = 10))
# initialize parameters with default values
alpha <- rep(0, nasset+1)
log_gamma <- rep(0, nasset+2)
mu <- rep(0, nasset+2)
log_sigma <- rep(0, nasset+2)
logit_rho <- rep(0, nasset+1)
logit_tau <- rep(0, nasset+1)
logit_omega <- rep(0, nasset)

# parameter and latent variable list for convenient updating
curr_par <- list(log_Vt = log_Vt,
                 alpha = alpha,
                 log_gamma = log_gamma,
                 mu = mu,
                 log_sigma = log_sigma,
                 logit_rho = logit_rho,
                 logit_tau = logit_tau,
                 logit_omega = logit_omega)

# optimization control parameters
opt_control <- list(eval.max = 1000, iter.max = 1000, trace = 10)


## ----init_opt_load, include = FALSE-------------------------------------------
if(params$load_calcs) curr_par <- readRDS("init_opt.RDS")

## ----init_opt, eval = !params$load_calcs--------------------------------------
#  # problem dimensions
#  nobs <- 1000 # number of days
#  nasset <- 5 # number of assets, excluding GSPC
#  
#  # format the data
#  dt <- 1/252
#  Xt <- as.matrix(cbind(GSPC = snp500[1:nobs, "GSPC"],
#                        snp500[1:nobs, 1+1:nasset]))
#  Xt <- log(Xt)
#  log_VPt <- log(as.numeric(snp500[1:nobs, "VIX"]))
#  
#  # initialize latent variables with rolling window standard deviations
#  log_Vt <- log(apply(Xt, 2, sv_init, dt = dt, block_size = 10))
#  # initialize parameters with default values
#  alpha <- rep(0, nasset+1)
#  log_gamma <- rep(0, nasset+2)
#  mu <- rep(0, nasset+2)
#  log_sigma <- rep(0, nasset+2)
#  logit_rho <- rep(0, nasset+1)
#  logit_tau <- rep(0, nasset+1)
#  logit_omega <- rep(0, nasset)
#  
#  # parameter and latent variable list for convenient updating
#  curr_par <- list(log_Vt = log_Vt,
#                   alpha = alpha,
#                   log_gamma = log_gamma,
#                   mu = mu,
#                   log_sigma = log_sigma,
#                   logit_rho = logit_rho,
#                   logit_tau = logit_tau,
#                   logit_omega = logit_omega)
#  
#  # optimization control parameters
#  opt_control <- list(eval.max = 1000, iter.max = 1000, trace = 10)
#  
#  # initialize parameters of each individual asset
#  tm_tot <- system.time({
#    for(iasset in 0:nasset) {
#      message("asset = ", iasset)
#      # construct model object
#      eou_ad <- eou_MakeADFun(Xt = Xt[,iasset+1], dt = dt,
#                              alpha = curr_par$alpha[iasset+1],
#                              log_Vt = curr_par$log_Vt[,iasset+1],
#                              log_gamma = curr_par$log_gamma[iasset+2],
#                              mu = curr_par$mu[iasset+2],
#                              log_sigma = curr_par$log_sigma[iasset+2],
#                              logit_rho = curr_par$logit_rho[iasset+1])
#      # optimize with quasi-newton method
#      tm <- system.time({
#        opt <- nlminb(start = eou_ad$par,
#                      objective = eou_ad$fn,
#                      gradient = eou_ad$gr,
#                      control = opt_control)
#      })
#      message("Time: ", round(tm[3], 2), " seconds")
#      # update parameters
#      curr_par <- svc_update(eou_ad, old_par = curr_par, iasset = iasset)
#    }
#  })
#  message("Total Time: ", round(tm_tot[3], 2), " seconds.")

## ----init_opt_save, include = FALSE-------------------------------------------
if(!params$load_calcs) saveRDS(curr_par, file = "init_opt.RDS")

## ----block_opt_load, include = FALSE------------------------------------------
if(params$load_calcs) curr_par <- readRDS("block_opt.RDS")

## ----block_opt, eval = !params$load_calcs-------------------------------------
#  nepoch <- 3
#  
#  tm_tot <- system.time({
#    for(iepoch in 1:nepoch) {
#      message("epoch = ", iepoch)
#      for(iasset in -1:nasset) {
#        message("asset = ", iasset)
#        svc_ad <- svc_MakeADFun(Xt = Xt, log_VPt = log_VPt, dt = dt,
#                                par_list = curr_par,
#                                iasset = iasset)
#        tm <- system.time({
#          opt <- nlminb(start = svc_ad$par,
#                        objective = svc_ad$fn,
#                        gradient = svc_ad$gr,
#                        control = opt_control)
#        })
#        message("Time: ", round(tm[3], 2), " seconds")
#        curr_par <- svc_update(svc_ad, old_par = curr_par, iasset = iasset)
#      }
#    }
#  })
#  message("Total Time: ", round(tm_tot[3], 2), " seconds.")

## ----block_opt_save, include = FALSE------------------------------------------
if(!params$load_calcs) saveRDS(curr_par, file = "block_opt.RDS")

## ----joint_opt_load, include = FALSE------------------------------------------
if(params$load_calcs) curr_par <- readRDS("joint_opt.RDS")

## ----joint_opt, eval = !params$load_calcs-------------------------------------
#  # joint parameter optimization
#  iasset <- "all"
#  svc_ad <- svc_MakeADFun(Xt = Xt, log_VPt = log_VPt, dt = dt,
#                          par_list = curr_par,
#                          iasset = iasset)
#  tm <- system.time({
#    opt <- nlminb(start = svc_ad$par,
#                  objective = svc_ad$fn,
#                  gradient = svc_ad$gr,
#                  control = opt_control)
#  })
#  message("Time: ", round(tm[3], 2), " seconds")
#  curr_par <- svc_update(svc_ad, old_par = curr_par, iasset = iasset)

## ----joint_opt_save, include = FALSE------------------------------------------
if(!params$load_calcs) saveRDS(curr_par, file = "joint_opt.RDS")

## ----summary_load, include = FALSE--------------------------------------------
if(params$load_calcs) svc_est <- readRDS("summary.RDS")

## ----summary, eval = !params$load_calcs---------------------------------------
#  system.time({
#    svc_est <- TMB::sdreport(svc_ad)
#  })
#  knitr::kable(summary(svc_est, select = "report"), digits = 2)

## ----summary_save, include = FALSE--------------------------------------------
if(!params$load_calcs) saveRDS(svc_est, file = "summary.RDS")

## ----forecast, fig.width = 7, fig.height = 4----------------------------------
nfwd <- 10 # number of days to forecast
nsim <- 1e4 # number of simulated forecasts

# sample from the mode-quadrature approximation
curr_post <- mvtnorm::rmvnorm(nsim,
                              mean = svc_est$value,
                              sigma = svc_est$cov)

# convert to appropriate inputs to svc_sim
sim_args <- sapply(c("alpha", "log_gamma", "mu", "log_sigma", "logit_rho",
                       "logit_tau", "logit_omega", "log_VT"),
                     function(nm) {
                       t(curr_post[,colnames(curr_post) == nm])
                     }, simplify = FALSE)
names(sim_args)[8] <- "log_V0" # rename log_VT
# remaining arguments
sim_args <- c(sim_args,
              list(X0 = Xt[nobs,],
                   log_VP0 = log_VPt[nobs],
                   nobs = nfwd, dt = dt))
# forward simulation
fwd_post <- do.call(svc_sim, args = sim_args)

# plot posteriors on day t = nobs + nfwd
# format data
X_fwd <- t(fwd_post$Xt[nfwd,,])
colnames(X_fwd) <- colnames(Xt)
X_fwd <- data.frame(type = "X", X_fwd)
V_fwd <- exp(t(fwd_post$log_Vt[nfwd,,]))
colnames(V_fwd) <- colnames(log_Vt)
V_fwd <- data.frame(type = "V", V_fwd)
XV_fwd <- rbind(X_fwd, V_fwd)
# plot
XV_fwd %>%
  as_tibble() %>%
  pivot_longer(GSPC:USB, names_to = "asset", values_to = "value") %>%
  mutate(asset = factor(asset, levels = colnames(Xt))) %>%
  ggplot(aes(x = value)) +
  geom_density() +
  facet_wrap(~ type:asset, scales = "free", nrow = 2) +
  xlab("") + ylab("Posterior Forecast Distribution")

