#' Initialize latent volatilities in a stochastic volatility model.
#'
#'@param Xt Vector of log-prices.
#'@param dt Interobservation time.
#'@param block_size Block size over which to calculate variance.
#'@return A vector of volatilities on the standard deviation scale, of the same length as `Xt`.
#'@export
sv_init <- function(Xt, dt, block_size) {
  dX <- diff(Xt) # log-returns
  N <- length(dX)
  ind <- rbind(floor(seq(1, N-block_size+1, len = N)),
               ceiling(seq(block_size, N, len = N)))
  Vt <- apply(ind, 2, function(ii) var(dX[ii[1]:ii[2]]))/dt
  sqrt(c(Vt[1], Vt))
}
