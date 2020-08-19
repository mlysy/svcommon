#' Construct a [TMB::MakeADFun()] object for the exponential Ornstein-Uhlenbeck stochastic volatility model.
#'
#' @param Xt Vector of `nobs` asset log prices.
#' @template param_dt
#' @param log_Vt Optional vector of `nobs` volatilities on the log standard deviation scale.  See 'Details'.
#' @param alpha Optional asset growth rate parameter.  See 'Details'.
#' @param log_gamma Optional log-volatility mean reversion parameter on the log scale.  See 'Details'.
#' @param mu Optional log-volatility mean parameter.  See 'Details'.
#' @param log_sigma Optional log-volatility diffusion parameter on the log scale.  See 'Details'.
#' @param logit_rho Optional correlation parameter between asset and volatility innovations, on the logit scale.  See 'Details'.
#' @param par_list Optional list with named elements consisting of a subset of `log_Vt`, `alpha`, `log_gamma`, `mu`, `log_sigma`, and `logit_rho`.  Values in `par_list` will supercede those of the corresponding individual argument if both are provided.
#' @param ... Additional arguments to [TMB::MakeADFun()].
#'
#' @return The result of a call to [TMB::MakeADFun()].
#'
#' @details The exponential Ornstein-Uhlenbeck (eOU) stochastic volatility model for a single asset is given by the stochastic differential equation (SDE)
#' ```
#' dlog_Vt = - gamma (log_Vt - mu) dt + sigma dB_Vt
#' dXt = (alpha - .5 Vt^2) dt + Vt (rho dB_Vt + sqrt(1-rho^2) dB_Zt),
#' ```
#' where `B_Vt` and `B_Zt` are independent Brownian motions.
#'
#' `eou_MakeADFun()` implements the Euler approximation to this SDE...
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
                          par_list, ...) {
  # extract arguments from par_list
  if(!missing(par_list)) {
    for(arg_name in c("log_Vt", "alpha", "log_gamma",
                      "mu", "log_sigma", "logit_rho")) {
      if(!is.null(par_list[[arg_name]])) {
        assign(arg_name, value = par_list[[arg_name]])
      }
    }
  }
  # argument checking
  check_vector(Xt = Xt)
  nobs <- length(Xt)
  check_scalar(dt = dt)
  log_Vt <- check_vector(log_Vt = log_Vt, len = nobs,
                         default = log(sv_init(Xt, dt, block_size = 10)))
  alpha <- check_scalar(alpha = alpha, default = 0)
  log_gamma <- check_scalar(log_gamma = log_gamma, default = 0)
  mu <- check_scalar(mu = mu, default = 0)
  log_sigma <- check_scalar(log_sigma = log_sigma, default = 0)
  logit_rho <- check_scalar(logit_rho = logit_rho, default = 0)
  # construct ADFun object
  par_list <- list(log_Vt = log_Vt, alpha = alpha, log_gamma = log_gamma,
                   mu = mu, log_sigma = log_sigma, logit_rho = logit_rho)
  TMB::MakeADFun(data = list(model = "sv_eou", Xt = Xt, dt = dt),
                 parameters = par_list,
                 random = "log_Vt",
                 DLL = "svcommon_TMBExports", silent = TRUE, ...)
}


#--- helper functions ----------------------------------------------------------

#' Check whether argument is a numeric scalar.
#'
#' @param ... Named argument.  The name is used to create an informative error message.
#' @param default An optional default value.  Function should throw an error if `...` is missing and no default is provided.
#' @return The scalar or its default value if the format is correct.
#' @noRd
check_scalar <- function(..., default) {
  ## xname <- names(sapply(match.call(), deparse)[-1])
  ## xname <- xname[!(xname %in% "default")]
  ## dots <- tryCatch(list(...), error = function(e) NA)
  ## if(is.na(dots)) x <- default else x <- dots[[1]]
  ## dots <- list(...)
  ## if(is.na(dots)) return(default)
  ## xname <- names(dots)
  ## x <- dots[[1]]
  xname <- names(sapply(match.call(), deparse)[-1])
  xname <- xname[!(xname %in% "default")]
  x <- tryCatch(list(...)[[1]], error = function(e) NULL)
  if(is.null(x)) {
    if(!missing(default)) {
      x <- default
    } else stop("object '", xname, "' not found.")
  }
  if(!is.numeric(x) || !is.vector(x) || length(x) != 1) {
    stop("'", xname, "' must be a numeric scalar.")
  }
  x
}

## #' @param x A named list of length 1.
## #' @param len An optional named list of length 1.
## #' @param default An optional default value.
## #' @details Should work as follows:
## #' - If `is.null(x)` and no default provided, says `names(x)` is missing with no default.
## #' - Otherwise, if `len` provided and doesn't match, says `names(x)` should have length `names(len)`.
## check_vector <- function(x, len, default) {

## }

#' Checks whether argument is a numeric vector of right length.
#'
#' @param ... Named argument.  The name is used to create an informative error message.
#' @param len Optional required length of vector.
#' @param default An optional default value.
#' @details Function should behave as follows:
#' - Error if `...` missing and no default provided.
#' - Error if `...` not a numeric vector.
#' - Error if `...` not a vector of right length.
#' - If `...` missing and default provided, set to default.
#' - If no errors encountered, return value.
#' @noRd
check_vector <- function(..., len, default) {
  ## if(is.na(dots)) {
  ##   if(length(default) == 1) default <- rep(default, len)
  ##   return(default)
  ## }
  ## xname <- names(dots)
  ## x <- dots[[1]]
  xname <- names(sapply(match.call(), deparse)[-1])
  xname <- xname[!(xname %in% c("len", "default"))]
  x <- tryCatch(list(...)[[1]], error = function(e) NULL)
  if(is.null(x)) {
    if(!missing(default)) {
      if(!missing(len) && length(default) == 1) default <- rep(default, len)
      x <- default
    } else stop("object '", xname, "' not found.")
  }
  if(!is.numeric(x) || !is.vector(x)) {
    stop("'", xname, "' must be a numeric vector.")
  }
  if(!missing(len) && (length(x) != len)) {
    stop("'", xname, "' has incorrect length.")
  }
  x
}

#' Checks whether argument is a numeric matrix of right dimensions.
#'
#' @param ... Named argument.  The name is used to create an informative error message.
#' @param dim Optional required dimensions.
#' @param default An optional default value.  Function should throw an error if `...` is missing and no default is provided.
#' @param promote If `TRUE`, vectors are first promoted to column matrices.
#' @noRd
check_matrix <- function(..., dim, default, promote = FALSE) {
  xname <- names(sapply(match.call(), deparse)[-1])
  xname <- xname[!(xname %in% c("dim", "default"))]
  x <- tryCatch(list(...)[[1]], error = function(e) NULL)
  if(is.null(x)) {
    if(!missing(default)) {
      if(!missing(dim) && !is.matrix(default) && length(default) == 1) {
        default <- array(default, dim = dim)
      }
      x <- default
    } else stop("object '", xname, "' not found.")
  }
  if(promote) x <- as.matrix(x)
  if(!is.numeric(x) || !is.matrix(x)) {
    stop("'", xname, "' must be a numeric matrix.")
  }
  if(!missing(dim) && !identical(base::dim(x), as.integer(dim))) {
    stop("'", xname, "' has incorrect dimensions.")
  }
  x
}

## #' Checks whether two array objects have the same dimensions.
## #'
## #' @param ... Vector of one or two named arguments, let's call them `x1` and `x2`.
## #' @param dims Optional vector of dimensions.
## #' @return Nothing, but throws an informative error if:
## #' - `...` contains two arguments and `dim(x1) != dim(x2)`.
## #' - `...` contains one argument and `dim(x1) != dims`.
## #' @noRd
## check_dims <- function(..., dims) {
##   dots <- list(...)
##   xn1 <- names(dots)[1]
##   x1 <- dots[[1]]
##   if(length(dots) == 1) {
##     if(!all(dim(as.array(x1)) == dims)) {
##       stop(xn1, " has incorrect dimensions.")
##     }
##   } else {
##     xn2 <- names(dots)[2]
##     x2 <- dots[[2]]
##     if(!all(dim(as.array(x1)) == dim(as.array(x2)))) {
##       stop(xn1, " and ", xn2, " have incompatible dimensions.")
##     }
##   }
## }
