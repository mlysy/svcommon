context("svc_sim")

## source("svcommon-testfunctions.R")

test_that("Innovations from svc_sim and svc_dB are the same.", {
  ntest <- 20
  for(ii in 1:ntest) {
    # simulation parameters
    nasset <- sample(1:3, 1)
    nobs <- sample(1:10, 1)
    nseries <- sample(1:5, 1)
    dt <- runif(1)/10
    alpha <- rmat(nasset + 1, nseries)
    log_gamma <- rmat(nasset + 2, nseries)
    mu <- rmat(nasset + 2, nseries)
    log_sigma <- rmat(nasset + 2, nseries)
    logit_rho <- rmat(nasset + 1, nseries)
    logit_tau <- rmat(nasset + 1, nseries)
    logit_omega <- rmat(nasset, nseries)
    X0 <- rmat(nasset + 1, nseries)
    log_VP0 <- rnorm(nseries)
    log_V0 <- rmat(nasset + 1, nseries)
    dB_Zt <- sqrt(dt) * rarr(nobs, nasset+1, nseries)
    dB_Vt <- sqrt(dt) * rarr(nobs, nasset+1, nseries)
    dB_VPt <- sqrt(dt) * rarr(nobs, nseries)
    # simulate data
    Yt <- svc_sim(nobs = nobs, dt = dt,
                  X0 = X0, log_VP0 = log_VP0, log_V0 = log_V0,
                  alpha = alpha,
                  log_gamma = log_gamma, mu = mu, log_sigma = log_sigma,
                  logit_rho = logit_rho,
                  logit_tau = logit_tau, logit_omega = logit_omega,
                  dBt = list(Z = dB_Zt, V = dB_Vt, VP = dB_VPt))
    # extract innovations
    for(iseries in 1:nseries) {
      # add back in the initial values
      Xt <- rbind(X0[,iseries], Yt$Xt[,,iseries])
      log_VPt <- c(log_VP0[iseries], Yt$log_VPt[,iseries])
      log_Vt <- rbind(log_V0[,iseries], Yt$log_Vt[,,iseries])
      # calculate innovations
      dB <- svc_dB(Xt = Xt,
                   log_VPt = log_VPt,
                   log_Vt = log_Vt,
                   dt = dt,
                   alpha = alpha[,iseries], log_gamma = log_gamma[,iseries],
                   mu = mu[,iseries], log_sigma = log_sigma[,iseries],
                   logit_rho = logit_rho[,iseries])
      # check VP
      expect_equal(dB$V[,1], dB_VPt[,iseries])
      # check VU
      expect_equal(dB_res(dB = dB$V[,-1], dB0 = dB$V[,1],
                          logit_rho = logit_tau[,iseries]),
                   matrix(dB_Vt[,,iseries], nobs, nasset+1))
      # check ZP
      expect_equal(dB$Z[,1], dB_Zt[,1,iseries])
      # check Zu
      expect_equal(dB_res(dB = dB$Z[,-1], dB0 = dB$Z[,1],
                          logit_rho = logit_omega[,iseries]),
                   matrix(dB_Zt[,-1,iseries], nobs, nasset))
    }
  }
})

#--- scratch -------------------------------------------------------------------

## dB_VPt[-1,iseries] %o% svcommon:::rho_itrans(logit_tau[,iseries]) +
##   sweep(dB_Vt[-1,,iseries], 2, sqrt(1-svcommon:::rho_itrans(logit_tau[,iseries])^2), FUN = "*")
## dB$V[,-1]

## dB_obs(dB = dB_Vt[,,iseries],
##        dB0 = dB_VPt[,iseries],
##        logit_rho = logit_tau[,iseries])
## dB$V[,-1]
