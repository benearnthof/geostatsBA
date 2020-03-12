#define STAN_MATH_PRIM_FUN_GP_MATERN32_COV_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/divide_columns.hpp>
#include <stan/math/prim/fun/distance.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <cmath>
#include <vector>

/**
 * Returns a Matern 3/2 covariance matrix
 *
 * \f[ k(x, x') = \sigma^2(1 +
 *  \frac{\sqrt{3}d(x, x')}{l})exp(-\frac{\sqrt{3}d(x, x')}{l})
 * \f]
 *
 * where \f$ d(x, x') \f$ is the Euclidean distance.
 *
 * @tparam T_x type for each scalar
 * @tparam T_s type of parameter of sigma
 * @tparam T_l type of parameter length scale
 *
 * @param x std::vector of scalars that can be used in squared distance
 * @param length_scale length scale
 * @param sigma marginal standard deviation or magnitude
 * @throw std::domain error if sigma <= 0, l <= 0, or x is nan or inf
 */
template <typename T_x, typename T_s, typename T_l>
inline typename Eigen::Matrix<return_type_t<T_x, T_s, T_l>, Eigen::Dynamic,
                              Eigen::Dynamic>
gp_matern32_cov(const std::vector<T_x> &x, const T_s &sigma,
                const T_l &length_scale) {
  using std::exp;
  using std::pow;
  
  size_t x_size = stan::math::size(x);
  Eigen::Matrix<return_type_t<T_x, T_s, T_l>, Eigen::Dynamic, Eigen::Dynamic>
    cov(x_size, x_size);
  
  if (x_size == 0) {
    return cov;
  }
  
  const char *function = "gp_matern32_cov";
  size_t x_obs_size = stan::math::size(x[0]);
  for (size_t i = 0; i < x_size; ++i) {
    check_size_match(function, "x row", x_obs_size, "x's other row",
                     stan::math::size(x[i]));
  }
  
  for (size_t n = 0; n < x_size; ++n) {
    check_not_nan(function, "x", x[n]);
  }
  
  check_positive_finite(function, "magnitude", sigma);
  check_positive_finite(function, "length scale", length_scale);
  
  T_s sigma_sq = square(sigma);
  T_l root_3_inv_l = std::sqrt(3.0) / length_scale;
  
  for (size_t i = 0; i < x_size; ++i) {
    cov(i, i) = sigma_sq;
    for (size_t j = i + 1; j < x_size; ++j) {
      return_type_t<T_x> dist = distance(x[i], x[j]);
      cov(i, j)
        = sigma_sq * (1.0 + root_3_inv_l * dist) * exp(-root_3_inv_l * dist);
      cov(j, i) = cov(i, j);
    }
  }
  return cov;
}