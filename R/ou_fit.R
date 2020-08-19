#' Estimate parameters of an Ornstein-Uhlenbeck process.
#'
#' Calculate the MLE using least-squares on the one-step Euler approximation.
#'
#' @param Xt Vector of log-prices.
#' @param dt Interobservation time.
#' @return A named vector containing the maximum likelihood estimate of the OU parameters.
#'
#' @details The Ornstein-Uhlenbeck process is described by the stochastic differential equation
#' ```
#' dXt = - gamma (Xt - mu) dt + sigma dBt.
#' ```
#' Here, the MLE is computed analytically from the Euler approximation
#' ```
#' X_{t+1} | X_t ~ N(X_t - gamma (X_t - mu) dt, sigma^2 dt).
#' ```
#' @export
ou_fit <- function(Xt, dt) {
  nobs <- length(Xt)
  M <- lm(diff(Xt) ~ Xt[-nobs])
  beta <- as.numeric(coef(M))
  gamma <- -beta[2]/dt
  mu <- beta[1]/(gamma*dt)
  sigma <- sqrt(mean(resid(M)^2)/dt)
  c(gamma = gamma, mu = mu, sigma = sigma)
}
