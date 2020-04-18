// returns a Matern 3/2 covariance matrix 
gp_matern32_cov(const std::vector<Eigen::Matrix<T_x, -1, 1>> &x,
                const T_s &sigma, const std::vector<T_l> &length_scale) {
  using std::exp;

  size_t x_size = stan::math::size(x);
  Eigen::Matrix<return_type_t<T_x, T_s, T_l>, Eigen::Dynamic, Eigen::Dynamic>
      cov(x_size, x_size);

  if (x_size == 0) {
    return cov;
  }
  const char *function = "gp_matern32_cov";
  size_t l_size = length_scale.size();
  for (size_t n = 0; n < x_size; ++n) {
    check_not_nan(function, "x", x[n]);
  }

  T_s sigma_sq = square(sigma);
  double root_3 = std::sqrt(3.0);

  std::vector<Eigen::Matrix<return_type_t<T_x, T_l>, -1, 1>> x_new
      = divide_columns(x, length_scale);

  for (size_t i = 0; i < x_size; ++i) {
    for (size_t j = i; j < x_size; ++j) {
      return_type_t<T_x, T_l> dist = distance(x_new[i], x_new[j]);
      cov(i, j) = sigma_sq * (1.0 + root_3 * dist) * exp(-root_3 * dist);
      cov(j, i) = cov(i, j);
    }
  }
  return cov;
}