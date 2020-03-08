// generated with brms 2.10.0
functions {
}
data {
  int<lower=1> N;  // number of observations
  vector[N] Y;  // response variable
  // data for splines
  int Ks;  // number of linear effects
  matrix[N, Ks] Xs;  // design matrix for the linear effects
  // data for spline s(x2)
  int nb_1;  // number of bases
  int knots_1[nb_1];  // number of knots
  // basis function matrices
  matrix[N, knots_1[1]] Zs_1_1;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  // temporary intercept for centered predictors
  real Intercept;
  // spline coefficients
  vector[Ks] bs;
  // parameters for spline s(x2)
  // standarized spline coefficients
  vector[knots_1[1]] zs_1_1;
  // standard deviations of the coefficients
  real<lower=0> sds_1_1;
  real<lower=0> sigma;  // residual SD
}
transformed parameters {
  // actual spline coefficients
  vector[knots_1[1]] s_1_1 = sds_1_1 * zs_1_1;
}
model {
  // initialize linear predictor term
  vector[N] mu = Intercept + rep_vector(0, N) + Xs * bs + Zs_1_1 * s_1_1;
  // priors including all constants
  target += student_t_lpdf(Intercept | 3, 8, 10);
  target += normal_lpdf(zs_1_1 | 0, 1);
  target += student_t_lpdf(sds_1_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += student_t_lpdf(sigma | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  // likelihood including all constants
  if (!prior_only) {
    target += normal_lpdf(Y | mu, sigma);
  }
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept;
}
