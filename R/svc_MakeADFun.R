#' Construct a [TMB::MakeADFun()] object for the exponential Ornstein-Uhlenbeck stochastic volatility model.
#'
#' @template param_Xt
#' @template param_log_VPt
#' @template param_dt
#' @template param_log_Vt
#' @template param_alpha
#' @template param_log_gamma
#' @template param_mu
#' @template param_log_sigma
#' @template param_logit_rho
#' @template param_logit_tau
#' @template param_logit_omega
#' @param par_list Optional list with named elements consisting of a subset of `log_Vt`, `alpha`, `log_gamma`, `mu`, `log_sigma`, `logit_rho`, `logit_tau`, and `logit_omega`.  Values in `par_list` will supercede those of the corresponding individual argument if both are provided.
#' @param iasset Index of asset for which parameters are to be treated as non-fixed.  Either the character string "all" indicating that no parameters are fixed, or an integer in `-1:nasset`, where `-1` denotes the proxy for the volatility factor, `0` denotes the proxy for the asset common factor, and `1:nasset` denotes the remaining assets.
#' @param ... Additional arguments to [TMB::MakeADFun()].
#'
#' @return The result of a call to [TMB::MakeADFun()].
#'
#' @details The common-factor multivariate stochastic volatility (SVC) model for multiple assets is given by the stochastic differential equation (SDE) ...
#'
#' `svc_MakeADFun` implements the Euler approximation to this SDE...
#'
#' The optional latent variable and parameter inputs `log_Vt`, `alpha`, ..., `logit_rho` can be set to initialize optimization routines.  The default values are for each parameter vector to consist of the zero vector of the appropriate length, and the columns of `log_Vt` to be the log of windowed standard deviation estimates for the corresponding asset as calculated by [sv_init()].
#'
#' `svc_MakeADFun` is a wrapper to [TMB::MakeADFun()].  This function may be called on the underlying C++ template provided by \pkg{svcommon} via
#' ```
#' TMB::MakeADFun(data = list(model = "sv_common", ...),
#'                parameters = list(...),
#'                DLL = "svcommon_TMBExports",
#'                ...)
#' ```
#' @export
svc_MakeADFun <- function(Xt, log_VPt, dt,
                          log_Vt, alpha, log_gamma, mu, log_sigma, logit_rho,
                          logit_tau, logit_omega, par_list,
                          iasset = "all", ...) {
  # extract arguments from par_list
  if(!missing(par_list)) {
    for(arg_name in c("log_Vt", "alpha", "log_gamma", "mu",
                      "log_sigma", "logit_rho", "logit_tau", "logit_omega")) {
      if(!is.null(par_list[[arg_name]])) {
        assign(arg_name, value = par_list[[arg_name]])
      }
    }
  }
  # argument checking
  check_matrix(Xt = Xt)
  nobs <- nrow(Xt)
  nasset <- ncol(Xt)-1
  check_vector(log_VPt = log_VPt, len = nobs)
  check_scalar(dt = dt)
  if(missing(log_Vt)) {
    log_Vt <- log(apply(Xt, 2, sv_init, dt = dt, block_size = 10))
  }
  log_Vt <- check_matrix(log_Vt = log_Vt, dim = c(nobs, nasset+1))
  alpha <- check_vector(alpha = alpha, len = nasset+1, default = 0)
  log_gamma <- check_vector(log_gamma = log_gamma,
                            len = nasset+2, default = 0)
  mu <- check_vector(mu = mu, len = nasset+2, default = 0)
  log_sigma <- check_vector(log_sigma = log_sigma, len = nasset+2, default = 0)
  logit_rho <- check_vector(logit_rho = logit_rho, len = nasset+1, default = 0)
  logit_tau <- check_vector(logit_tau = logit_tau, len = nasset+1, default = 0)
  logit_omega <- check_vector(logit_omega = logit_omega, len = nasset,
                              default = 0)
  # construct ADFun object
  par_list <- list(log_Vt = log_Vt, alpha = alpha, log_gamma = log_gamma,
                   mu = mu, log_sigma = log_sigma, logit_rho = logit_rho,
                   logit_tau = logit_tau, logit_omega = logit_omega)
  # for blockwise coord descent,
  # checked that T and F gives same result but former much faster
  # for joint optimization unfix everything
  fix_Vt <- iasset != "all"
  map_list <- svc_map(iasset = iasset, nasset = nasset, nobs = nobs,
                      fix_Vt = fix_Vt)
  TMB::MakeADFun(data = list(model = "sv_common",
                             Xt = Xt, log_VPt = log_VPt, dt = dt),
                 parameters = par_list,
                 random = "log_Vt",
                 map = map_list,
                 DLL = "svcommon_TMBExports", silent = TRUE, ...)
}

#--- helper fuctions -----------------------------------------------------------

#' Create map to fix certain parameters.
#'
#' @param iasset,nasset,nobs,fix_Vt See [svc_MakeADFun()].
#' @details With `fix_Vt = TRUE`, both `iasset` and the common asset (`iasset = 0`) are the only ones to get updated.
#' @noRd
svc_map <- function(iasset, nasset, nobs, fix_Vt) {
  if(iasset == "all") {
    map_list <- list()
    if(fix_Vt) {
      map_list <- c(map_list, list(log_Vt = factor(rep(NA, nobs*(nasset+1)))))
    }
  } else if(iasset %in% -1:nasset) {
    map_list <- list(alpha = 0:nasset,
                     log_gamma = -1:nasset,
                     mu = -1:nasset,
                     log_sigma = -1:nasset,
                     logit_rho = 0:nasset,
                     logit_tau = 0:nasset,
                     logit_omega = 1:nasset)
    # all NA except `iasset`
    map_list <- lapply(map_list, function(x) {
      x[x != iasset] <- NA
      factor(x, levels = if(all(is.na(x))) NULL else iasset)
    })
    if(fix_Vt) {
      # all but individual asset and common asset are fixed
      log_Vt <- matrix(1:(nobs*(nasset+1)), nobs, nasset+1)
      log_Vt[,!(0:nasset %in% c(0, iasset))] <- NA
      log_Vt <- factor(as.numeric(log_Vt))
      map_list <- c(map_list, list(log_Vt = log_Vt))
    }
  } else {
    stop("Invalid specification of iasset.")
  }
  map_list
}
