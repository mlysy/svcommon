#--- jacobian of the logit prior -----------------------------------------------

nu_fun <- function(rho) log((1+rho)/(1-rho))
rho_fun <- function(nu) (exp(nu) - 1)/(1 + exp(nu))
rho_fun2 <- function(nu) 2/(1 + exp(-nu)) - 1
rho_grad <- function(nu) 2 * exp(nu) / (1 + exp(nu))^2
rho_grad2 <- function(nu) 2 * exp(-nu - 2 * log(1 + exp(-nu)))

rho <- runif(10, -1, 1)
rho - rho_fun2(nu_fun(rho))

curve(rho_fun, from = -10, to = 10)

sapply(nu_fun(rho), numDeriv::grad, func = rho_fun) -
  rho_grad2(nu_fun(rho))
