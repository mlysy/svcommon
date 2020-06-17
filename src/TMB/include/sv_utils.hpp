/// @file sv_utils.hpp
/// @brief helper functions for sv models.

#ifndef svcommon_sv_utils_hpp
#define svcommon_sv_utils_hpp

// for passing arguments to Eigen functions
template <class Type>
using RefMatrix = Eigen::Ref <Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> >;
template <class Type>
using cRefMatrix = const Eigen::Ref <const Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> >;
// template <class Type>
// using RefVector = Eigen::Ref <Eigen::Vector<Type, Eigen::Dynamic> >;
// template <class Type>
// using cRefVector = const Eigen::Ref <const Eigen::Matrix<Type, Eigen::Dynamic> >;

/// @brief Calculate mean and standard deviation of OU model.
///
/// @param[out] mean Vector of `nobs` OU means.
/// @param[out] sd Vector of `nobs` OU standard deviations.
/// @param[in] logv Vector of `nobs` OU observations.
/// @param[in] dt Interobservation time.
/// @param[in] gamma Mean reversion parameter.
/// @param[in] mu Mean parameter.
/// @param[in] sigma Diffusion parameter.
template <class Type>
void ou_ms(RefMatrix<Type> mean, RefMatrix<Type> sd,
	   cRefMatrix<Type>& logv, Type dt,
	   Type gamma, Type mu, Type sigma) {
  mean = logv.array() - gamma * (logv.array() - mu) * dt;
  sd.setConstant(sigma);
  return;
}

/// @brief Calculate mean and standard deviation of Brownian motion model.
///
/// @param[out] mean Vector of `nobs` BM means.
/// @param[out] sd Vector of `nobs` BM standard deviations.
/// @param[in] x Vector of `nobs` BM observations.
/// @param[in] logv Vector of `nobs` OU observations.
/// @param[in] dt Interobservation time.
/// @param[in] alpha Growth rate parameter.
template <class Type>
void bm_ms(RefMatrix<Type> mean, RefMatrix<Type> sd,
	   cRefMatrix<Type>& x, cRefMatrix<Type>& logv, Type dt,
	   Type alpha) {
  sd = logv.array().exp();
  mean = x.array() + (alpha - Type(.5) * sd.array() * sd.array()) * dt;
  return;
}

/// @brief Calculate residual.
///
/// @param[out] z Residual vector of size `n`.
/// @param[in] x Original vector of size `n`.
/// @param[in] mu Mean vector of size `n`.
/// @param[in] sigma Standard deviation vector of size `n`.
template <class Type>
void residual(RefMatrix<Type> z, cRefMatrix<Type>& x,
	      cRefMatrix<Type>& mu, cRefMatrix<Type> sigma) {
  z = (x - mu).array() / sigma.array();
  return;
}

template <class Type>
Type corr_sqm(Type rho) {
  return sqrt(Type(1.0) - rho*rho);
}

#endif
