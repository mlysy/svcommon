#--- helper functions for checking inputs --------------------------------------

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
#' @param promote If `TRUE` and `len` is provided, scalars are replicated as necessary.
#' @details Function should behave as follows:
#' - Error if `...` missing and no `default` provided.
#' - Error if `...` not a numeric vector.
#' - Error if `len` is provided and `...` is not of right length.
#' - If `...` missing and `default` provided, argument is set to `default` and the above checks are made.
#' - If `default` is a scalar and `len` is provided, `default` is set to `rep(default, len)` and the above checks are made.
#' - If no errors encountered, return value.
#' @noRd
check_vector <- function(..., len, default, promote = FALSE) {
  ## if(is.na(dots)) {
  ##   if(length(default) == 1) default <- rep(default, len)
  ##   return(default)
  ## }
  ## xname <- names(dots)
  ## x <- dots[[1]]
  xname <- names(sapply(match.call(), deparse)[-1])
  xname <- xname[!(xname %in% c("len", "default", "promote"))]
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
  if(promote && !missing(len) && (length(x) == 1)) {
      x <- rep(x, len)
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
#' @param default An optional default value.
#' @param promote If `TRUE`, vectors are first promoted to column matrices.  If `dim` is specified, vectors of length `dim[1]` are replicated as needed.
#' @details Essentially the same behavior as [check_vector()].
#' @noRd
check_matrix <- function(..., dim, default, promote = FALSE) {
  xname <- names(sapply(match.call(), deparse)[-1])
  xname <- xname[!(xname %in% c("dim", "default", "promote"))]
  x <- tryCatch(list(...)[[1]], error = function(e) NULL)
  if(is.null(x)) {
    if(!missing(default)) {
      if(!missing(dim) && !is.matrix(default) && length(default) == 1) {
        default <- array(default, dim = dim)
      }
      x <- default
    } else stop("object '", xname, "' not found.")
  }
  if(promote) {
    if(!missing(dim) && is.vector(x) && (length(x) == dim[1])) {
      dn <- if(is.null(names(x))) NULL else list(names(x), NULL)
      x <- matrix(x, nrow = dim[1], ncol = dim[2], dimnames = dn)
    } else {
      x <- as.matrix(x)
    }
  }
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
