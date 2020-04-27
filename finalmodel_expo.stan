// generated with brms 2.12.0
functions {

  /* compute a latent Gaussian process
   * Args:
   *   x: array of continuous predictor values
   *   sdgp: marginal SD parameter
   *   lscale: length-scale parameter
   *   zgp: vector of independent standard normal variables 
   * Returns:  
   *   a vector to be added to the linear predictor
   */ 
  vector gp(vector[] x, real sdgp, vector lscale, vector zgp) { 
    int Dls = rows(lscale);
    int N = size(x);
    matrix[N, N] cov;
    if (Dls == 1) {
      // one dimensional or isotropic GP
      cov = cov_exp_quad(x, sdgp, lscale[1]);
    } else {
      // multi-dimensional non-isotropic GP
      cov = cov_exp_quad(x[, 1], sdgp, lscale[1]);
      for (d in 2:Dls) {
        cov = cov .* cov_exp_quad(x[, d], 1, lscale[d]);
      }
    }
    for (n in 1:N) {
      // deal with numerical non-positive-definiteness
      cov[n, n] += 1e-8;
    }
    return cholesky_decompose(cov) * zgp;
  }

  /* Spectral density function of a Gaussian process
   * Args:
   *   x: array of numeric values of dimension NB x D
   *   sdgp: marginal SD parameter
   *   lscale: vector of length-scale parameters
   * Returns: 
   *   numeric values of the function evaluated at 'x'
   */
  vector spd_cov_exp_quad(vector[] x, real sdgp, vector lscale) {
    int NB = dims(x)[1];
    int D = dims(x)[2];
    int Dls = rows(lscale);
    vector[NB] out;
    if (Dls == 1) {
      // one dimensional or isotropic GP
      real constant = square(sdgp) * (sqrt(2 * pi()) * lscale[1])^D;
      real neg_half_lscale2 = -0.5 * square(lscale[1]);
      for (m in 1:NB) {
        out[m] = constant * exp(neg_half_lscale2 * dot_self(x[m]));
      }
    } else {
      // multi-dimensional non-isotropic GP
      real constant = square(sdgp) * sqrt(2 * pi())^D * prod(lscale);
      vector[Dls] neg_half_lscale2 = -0.5 * square(lscale);
      for (m in 1:NB) {
        out[m] = constant * exp(dot_product(neg_half_lscale2, square(x[m])));
      }
    }
    return out;
  }
  /* compute an approximate latent Gaussian process
   * Args:
   *   X: Matrix of Laplacian eigen functions at the covariate values
   *   sdgp: marginal SD parameter
   *   lscale: vector of length-scale parameters
   *   zgp: vector of independent standard normal variables 
   *   slambda: square root of the Laplacian eigen values
   * Returns:  
   *   a vector to be added to the linear predictor
   */ 
  vector gpa(matrix X, real sdgp, vector lscale, vector zgp, vector[] slambda) { 
    vector[cols(X)] diag_spd = sqrt(spd_cov_exp_quad(slambda, sdgp, lscale));
    return X * (diag_spd .* zgp);
  }
}
data {
  int<lower=1> N;  // number of observations
  int Y[N];  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // data for splines
  int Ks;  // number of linear effects
  matrix[N, Ks] Xs;  // design matrix for the linear effects
  // data for spline s(dem)
  int nb_1;  // number of bases
  int knots_1[nb_1];  // number of knots
  // basis function matrices
  matrix[N, knots_1[1]] Zs_1_1;
  // data for spline s(rain)
  int nb_2;  // number of bases
  int knots_2[nb_2];  // number of knots
  // basis function matrices
  matrix[N, knots_2[1]] Zs_2_1;
  // data for spline s(tpi)
  int nb_3;  // number of bases
  int knots_3[nb_3];  // number of knots
  // basis function matrices
  matrix[N, knots_3[1]] Zs_3_1;
  // data related to GPs
  int<lower=1> Kgp_1;  // number of sub-GPs (equal to 1 unless 'by' was used)
  int<lower=1> Dgp_1;  // GP dimension
  // number of basis functions of an approximate GP
  int<lower=1> NBgp_1;
  // approximate GP basis matrices
  matrix[N, NBgp_1] Xgp_1;
  // approximate GP eigenvalues
  vector[Dgp_1] slambda_1[NBgp_1];
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  vector[Ks] bs;  // spline coefficients
  // parameters for spline s(dem)
  // standarized spline coefficients
  vector[knots_1[1]] zs_1_1;
  real<lower=0> sds_1_1;  // standard deviations of spline coefficients
  // parameters for spline s(rain)
  // standarized spline coefficients
  vector[knots_2[1]] zs_2_1;
  real<lower=0> sds_2_1;  // standard deviations of spline coefficients
  // parameters for spline s(tpi)
  // standarized spline coefficients
  vector[knots_3[1]] zs_3_1;
  real<lower=0> sds_3_1;  // standard deviations of spline coefficients
  vector<lower=0>[Kgp_1] sdgp_1;  // GP standard deviation parameters
  vector<lower=0>[1] lscale_1[Kgp_1];  // GP length-scale parameters
  vector[NBgp_1] zgp_1;  // latent variables of the GP
}
transformed parameters {
  // actual spline coefficients
  vector[knots_1[1]] s_1_1;
  // actual spline coefficients
  vector[knots_2[1]] s_2_1;
  // actual spline coefficients
  vector[knots_3[1]] s_3_1;
  // compute actual spline coefficients
  s_1_1 = sds_1_1 * zs_1_1;
  // compute actual spline coefficients
  s_2_1 = sds_2_1 * zs_2_1;
  // compute actual spline coefficients
  s_3_1 = sds_3_1 * zs_3_1;
}
model {
  // initialize linear predictor term
  vector[N] mu = Intercept + Xc * b + Xs * bs + Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1 + Zs_3_1 * s_3_1 + gpa(Xgp_1, sdgp_1[1], lscale_1[1], zgp_1, slambda_1);
  // priors including all constants
  target += student_t_lpdf(Intercept | 3, 0, 10);
  target += student_t_lpdf(sds_1_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(zs_1_1 | 0, 1);
  target += student_t_lpdf(sds_2_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(zs_2_1 | 0, 1);
  target += student_t_lpdf(sds_3_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(zs_3_1 | 0, 1);
  target += student_t_lpdf(sdgp_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += inv_gamma_lpdf(lscale_1[1][1] | 0.494077, 7e-05);
  target += normal_lpdf(zgp_1 | 0, 1);
  // likelihood including all constants
  if (!prior_only) {
    target += bernoulli_logit_lpmf(Y | mu);
  }
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}
