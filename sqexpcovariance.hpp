#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/distance.hpp>
#include <stan/math/prim/fun/divide_columns.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <cmath>
#include <vector>


template <typename T_x, typename T_s, typename T_l>
inline typename Eigen::Matrix<return_type_t<T_x, T_s, T_l>, Eigen::Dynamic,
                              Eigen::Dynamic>
gp_exponential_cov(const std::vector<T_x> &x, const T_s &sigma,
                   const T_l &length_scale) {
  using std::exp;
  using std::pow;
  
  size_t x_size = stan::math::size(x);
  Eigen::Matrix<return_type_t<T_x, T_s, T_l>, Eigen::Dynamic, Eigen::Dynamic>
    cov(x_size, x_size);
  if (x_size == 0) {
    return cov;
  }
  
  const char *function = "gp_exponential_cov";
  size_t x_obs_size = stan::math::size(x[0]);
  for (size_t i = 0; i < x_size; ++i) {
    check_size_match(function, "x row", x_obs_size, "x's other row",
                     stan::math::size(x[i]));
  }
  
  for (size_t i = 0; i < x_size; ++i) {
    check_not_nan(function, "x", x[i]);
  }
  
  check_positive_finite(function, "magnitude", sigma);
  check_positive_finite(function, "length scale", length_scale);
  
  T_s sigma_sq = square(sigma);
  T_l neg_inv_l = -1.0 / length_scale;
  
  for (size_t i = 0; i < x_size; ++i) {
    cov(i, i) = sigma_sq;
    for (size_t j = i + 1; j < x_size; ++j) {
      cov(i, j) = sigma_sq * exp(neg_inv_l * distance(x[i], x[j]));
      cov(j, i) = cov(i, j);
    }
  }
  return cov;
}