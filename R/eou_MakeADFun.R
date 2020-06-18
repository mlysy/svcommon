#' Construct a [TMB::MakeADFun()] object for the exponential Ornstein-Uhlenbeck stochastic volatility model.
#'
#' @param Xt Vector of `nobs` asset log prices.
#' @param dt Interobservation time.
#' @param log_Vt Optional vector of `nobs` volatilities on the log standard deviation scale.  See 'Details'.
#' @param alpha Optional asset growth rate parameter.  See 'Details'.
#' @param log_gamma Optional log-volatility mean reversion parameter on the log scale.  See 'Details'.
#' @param mu Optional log-volatility mean parameter.  See 'Details'.
#' @param log_sigma Optional log-volatility diffusion parameter on the log scale.  See 'Details'.
#' @param logit_rho Optional correlation parameter between asset and volatility innovations, on the logit scale.  See 'Details'.
#' @param par_list Optional list with named elements consisting of a subset of `log_Vt`, `alpha`, `log_gamma`, `mu`, `log_sigma`, and `logit_rho`.  Values in `par_list` will supercede those of the corresponding individual argument if both are provided.
#'
#' @return The result of a call to [TMB::MakeADFun()].
#'
#' @details The exponential Ornstein-Uhlenbeck (eOU) stochastic volatility model for a single asset is given by the stochastic differential equation (SDE) ...
#'
#' `eou_MakeADFun` implements the Euler approximation to this SDE...
#'
#' The optional inputs `log_Vt`, `alpha`, ..., `logit_rho` can be set to initialize optimization routines.  The default values are `alpha = 0`, ..., `logit_rho = 0`, and `log_Vt` as the log of windowed standard deviation estimates returned by [sv_init()].
#'
#' `eou_MakeADFun` is a wrapper to [TMB::MakeADFun()].  This function may be called on the underlying C++ template provided by \pkg{svcommon} via
#' ```
#' TMB::MakeADFun(data = list(model = "sv_eou", ...),
#'                parameters = list(...),
#'                DLL = "svcommon_TMBExports",
#'                ...)
#' ```
#' @export
eou_MakeADFun <- function(Xt, dt,
                          log_Vt, alpha, log_gamma, mu, log_sigma, logit_rho,
                          par_list) {
  # argument formatting
  check_vector(Xt = Xt)
  nobs <- length(Xt)
  check_scalar(dt = dt)
  # extract arguments from par_list
  if(!missing(par_list)) {
    for(arg_name in c("log_Vt", "alpha", "log_gamma",
                      "mu", "log_sigma", "logit_rho")) {
      if(!is.null(par_list[[arg_name]])) {
        assign(arg_name, value = par_list[[arg_name]])
      }
    }
  }
  # check arguments
  log_Vt <- check_vector(log_Vt = log_Vt,
                         default = log(sv_init(Xt, dt, block_size = 10)))
  check_dims(log_Vt = log_Vt, Xt = Xt)
  alpha <- check_scalar(alpha = alpha, default = 0)
  log_gamma <- check_scalar(log_gamma = log_gamma, default = 0)
  mu <- check_scalar(mu = mu, default = 0)
  log_sigma <- check_scalar(log_sigma = log_sigma, default = 0)
  logit_rho <- check_scalar(logit_rho = logit_rho, default = 0)
  # construct ADFun object
  par_list <- list(log_Vt = log_Vt, alpha = alpha, log_gamma = log_gamma,
                   mu = mu, log_sigma = log_sigma, logit_rho = logit_rho)
  MakeADFun(data = list(model = "sv_eou", Xt = Xt, dt = dt),
            parameters = par_list,
            random = "log_Vt",
            DLL = "svcommon_TMBExports", silent = TRUE)
}


#--- helper functions ----------------------------------------------------------

check_scalar <- function(..., default) {
  ## if(!missing(default) && missing(...)) return(default)
  ## if(!missing(default)) return(default)
  dots <- list(...)
  if(is.na(dots)) return(default)
  xname <- names(dots)
  ## x <- if(!missing(default)) default else dots[[1]]
  x <- dots[[1]]
  if(!is.numeric(x) || !is.vector(x) || length(x) != 1) {
    stop("`", xname, "` must be a numeric scalar.")
  }
  x
}


# behavior is as follows
check_vector <- function(..., default) {
  ## if(!missing(default) && missing(...)) return(default)
  ## if(!missing(default)) return(default)
  dots <- tryCatch(list(...), error = function(e) NA)
  if(is.na(dots)) return(default)
  xname <- names(dots)
  x <- dots[[1]]
  ## x <- if(!missing(default)) default else dots[[1]]
  if(!is.numeric(x) || !is.vector(x)) {
    stop("`", xname, "` must be a numeric vector.")
  }
  x
}

check_dims <- function(...) {
  dots <- list(...)
  xn1 <- names(dots)[1]
  xn2 <- names(dots)[2]
  x1 <- dots[[1]]
  x2 <- dots[[2]]
  if(!all(dim(as.array(x1)) == dim(as.array(x2)))) {
    stop(xn1, " and ", xn2, " have incompatible dimensions.")
  }
}
