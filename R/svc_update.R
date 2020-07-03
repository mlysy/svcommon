#' Update the parameters and latent variables of an SVC model.
#'
#' @param fit An object returned by [svc_MakeADFun()] or [eou_MakeADFun()], typically after it has been used in a numerical optimization routine.
#' @param old_par A list of SVC parameters to be updated.  Must be in the same format as the `par_list` argument of [svc_MakeADFun()].
#' @param iasset Which parameter(s) to update.  See corresponding argument of [svc_MakeADFun()].
#' @return A list in the same format as `old_par` but with the corresponding parameters and latent variables updated with those of `svc_fit`.
#' @export
svc_update <- function(fit, old_par, iasset) {
  # problem dimensions
  nobs <- nrow(old_par$log_Vt)
  nasset <- ncol(old_par$log_Vt)-1
  # extract parameters and latent variables from svc_fit
  new_par <- fit$env$last.par.best
  new_par <- sapply(unique(names(new_par)), function(pn) {
    setNames(new_par[names(new_par) == pn], nm = NULL)
  }, simplify = FALSE)
  new_Vt <- new_par$log_Vt
  new_par <- new_par[names(new_par) != "log_Vt"]
  # parameter updates
  ## if(iasset == "all") stop('iasset == "all" currently not supported.')
  map_list <- svc_map(iasset, nasset, fix_Vt = FALSE)
  for(nm in names(new_par)) {
    if(iasset == "all") {
      old_par[[nm]] <- new_par[[nm]]
    } else {
      old_par[[nm]][!is.na(map_list[[nm]])] <- new_par[[nm]]
    }
  }
  # volatility update
  if(length(new_Vt) == nobs) {
    old_par$log_Vt[,iasset+1] <- new_Vt
  } else if(length(new_Vt) == nobs * (nasset+1)) {
    old_par$log_Vt <- matrix(new_Vt, nobs, nasset+1)
  }
  old_par
}
