/// @filename sv_eou.hpp
/// @brief The exponential Ornstein-Uhlenbeck stochastic volatility model.

// #include <TMB.hpp>
#include "include/sv_utils.hpp"

/// @param[in] Xt Log-asset vector of size `n_obs`.
/// @param dt Interobservation time.
/// @param log_Vt Log-volatility vector of size `n_obs`.
/// @param alpha Instantaneous growth rate (scalar).
/// @param log_gamma Mean-reversion parameter (log scale).
/// @param mu Mean log-volatility.
/// @param log_sigma Volatility diffusion (log scale).
/// @param logit_rho Correlation parameter between log-asset and log-volatility innovations (logit scale).
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type>
Type sv_eou(objective_function<Type>* obj) {
  DATA_VECTOR(Xt);
  DATA_SCALAR(dt);
  PARAMETER_VECTOR(log_Vt); 
  PARAMETER(alpha);
  PARAMETER(log_gamma);
  PARAMETER(mu);
  PARAMETER(log_sigma);
  PARAMETER(logit_rho);
  // constants
  int n_obs = Xt.size();
  Type sqrt_dt = sqrt(dt);
  // parameter conversions
  Type gamma = exp(log_gamma);
  Type sigma = exp(log_sigma);
  Type rho = Type(2.0) / (Type(1.0) + exp(-logit_rho)) - Type(1.0);
  // temporary storage
  vector<Type> sde_mean(n_obs-1);
  vector<Type> sde_sd(n_obs-1);
  vector<Type> dB_V(n_obs-1);
  vector<Type> dB_Z(n_obs-1);
  Type rho_sqm;
  // output
  Type llik = Type(0.0);
  Type log_VT; // log-volatility at last timepoint
  // volatility
  ou_ms<Type>(sde_mean, sde_sd, log_Vt.segment(0, n_obs-1), dt,
	      gamma, mu, sigma);
  residual<Type>(dB_V, log_Vt.segment(1, n_obs-1), sde_mean, sde_sd);
  llik += (dnorm(dB_V, Type(0.0), sqrt_dt, 1) - log(sde_sd)).sum();
  // log-asset
  bm_ms<Type>(sde_mean, sde_sd,
	      Xt.segment(0, n_obs-1), log_Vt.segment(0, n_obs-1), dt, alpha);
  rho_sqm = corr_sqm(rho);
  sde_mean += sde_sd.array() * dB_V.array() * rho;
  sde_sd.array() *= rho_sqm;
  residual<Type>(dB_Z, Xt.segment(1, n_obs-1), sde_mean, sde_sd);
  llik += (dnorm(dB_Z, Type(0.0), sqrt_dt, 1) - log(sde_sd)).sum();
  // add log_VT to the output
  log_VT = log_Vt(n_obs-1);
  ADREPORT(log_VT);
  // add all parameters to report
  ADREPORT(alpha);
  ADREPORT(log_gamma);
  ADREPORT(mu);
  ADREPORT(log_sigma);
  ADREPORT(logit_rho);
  return -llik;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
