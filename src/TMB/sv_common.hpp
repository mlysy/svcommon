/// @filename sv_common.hpp
/// @brief A scalable common-factor stochastic volatility model.

// #include <TMB.hpp>
#include "include/sv_utils.hpp"

/// @param[in] Xt Log-return matrix of size `n_obs x n_stocks`.  The first row is the stock corresponding to the common factor.
/// @param[in] log_VPt Volatility proxy, on the log scale.  Vector of length `n_obs`.
/// @param dt Interobservation time.
/// @param log_Vt Log-volatility matrix of size `n_obs x n_stocks`.
/// @param alpha Vector of `n_stocks` instantaneous growth rates.
/// @param log_gamma Vector of `n_stocks + 1` mean-reversion parameters (log scale).  The first is that of the common factor, the last is that of the volatility proxy.
/// @param mu Vector of `n_stocks + 1` mean log-volatilities.
/// @param log_sigma Vector of `n_stocks + 1` volatility diffusions (log scale).
/// @param logit_rho Vector of `n_stocks` correlation parameters between log-return and log-volatility innovations (logit scale).
/// @param logit_tau Vector of `n_stocks` correlation parameters between log-volatilty and log-proxy innovations (logit scale).
/// @param logit_omega Vector of `n_stocks - 1` correlation parameters between common factor log-return innovation and other log-return residual innovations.
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type>
Type sv_common(objective_function<Type>* obj) {
  DATA_MATRIX(Xt);
  DATA_VECTOR(log_VPt);
  DATA_SCALAR(dt);
  PARAMETER_MATRIX(log_Vt); 
  PARAMETER_VECTOR(alpha);
  PARAMETER_VECTOR(log_gamma);
  PARAMETER_VECTOR(mu);
  PARAMETER_VECTOR(log_sigma);
  PARAMETER_VECTOR(logit_rho);
  PARAMETER_VECTOR(logit_tau);
  PARAMETER_VECTOR(logit_omega);
  // constants
  int n_obs = Xt.rows();
  int n_stocks = Xt.cols();
  Type sqrt_dt = sqrt(dt);
  // parameter conversions
  vector<Type> gamma = exp(log_gamma);
  vector<Type> sigma = exp(log_sigma);
  vector<Type> rho = Type(2.0) / (Type(1.0) + exp(-logit_rho)) - Type(1.0);
  vector<Type> tau = Type(2.0) / (Type(1.0) + exp(-logit_tau)) - Type(1.0);
  vector<Type> omega = Type(2.0) / (Type(1.0) + exp(-logit_omega)) - Type(1.0);
  // temporary storage
  vector<Type> sde_mean(n_obs-1);
  vector<Type> sde_sd(n_obs-1);
  vector<Type> dB_VP(n_obs-1);
  vector<Type> dB_V(n_obs-1);
  vector<Type> dB_Z0(n_obs-1);
  vector<Type> dB_Z(n_obs-1);
  Type rho_sqm;
  // output
  Type llik = Type(0.0);
  vector<Type> log_VT(n_stocks); // log-volatilities at last timepoint
  // volatility proxy
  ou_ms<Type>(sde_mean, sde_sd, log_VPt.segment(0, n_obs-1), dt,
	      gamma(0), mu(0), sigma(0));
  residual<Type>(dB_VP, log_VPt.segment(1, n_obs-1), sde_mean, sde_sd);
  // dB_VP = (log_VPt.segment(1, n_obs-1) - sde_mean).array() / sde_sd.array();
  REPORT(sde_mean);
  REPORT(sde_sd);
  REPORT(dB_VP);
  llik += (dnorm(dB_VP, Type(0.0), sqrt_dt, 1) - log(sde_sd)).sum();
  for(int ii=0; ii<n_stocks; ii++) {
    // volatilities
    ou_ms<Type>(sde_mean, sde_sd, log_Vt.block(0, ii, n_obs-1, 1), dt,
  		gamma(ii+1), mu(ii+1), sigma(ii+1));
    residual<Type>(dB_V, log_Vt.block(1, ii, n_obs-1, 1), sde_mean, sde_sd);
    sde_mean = tau(ii) * dB_VP.array();
    llik += (dnorm(dB_V, sde_mean,
  		   sqrt_dt * corr_sqm(tau(ii)), 1) - log(sde_sd)).sum();
    // assets
    bm_ms<Type>(sde_mean, sde_sd,
		Xt.block(0, ii, n_obs-1, 1),
		log_Vt.block(0, ii, n_obs-1, 1), dt, alpha(ii));
    rho_sqm = corr_sqm(rho(ii));
    if(ii == 0) {
      // common-factor asset
      sde_mean += sde_sd.array() * dB_V.array() * rho(ii);
      sde_sd.array() *= rho_sqm;
      residual<Type>(dB_Z0, Xt.block(1, ii, n_obs-1, 1), sde_mean, sde_sd);
      llik += (dnorm(dB_Z0, Type(0.0), sqrt_dt, 1) - log(sde_sd)).sum();
    } else {
      // other assets
      sde_mean += sde_sd.array() * (rho(ii) * dB_V.array() + rho_sqm * omega(ii-1) * dB_Z0.array());
      sde_sd.array() *= rho_sqm * corr_sqm(omega(ii-1));
      residual<Type>(dB_Z, Xt.block(1, ii, n_obs-1, 1), sde_mean, sde_sd);
      llik += (dnorm(dB_Z, Type(0.0), sqrt_dt, 1) - log(sde_sd)).sum();
    }
  }
  // add log_VT to the output
  log_VT = log_Vt.row(n_obs-1);
  ADREPORT(log_VT);
  // add all parameters to report
  ADREPORT(alpha);
  ADREPORT(log_gamma);
  ADREPORT(mu);
  ADREPORT(log_sigma);
  ADREPORT(logit_rho);
  ADREPORT(logit_tau);
  ADREPORT(logit_omega);
  return -llik;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
