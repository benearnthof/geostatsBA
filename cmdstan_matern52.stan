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
      cov = gp_matern52_cov(x, sdgp, lscale[1]);
    } else {
      // multi-dimensional non-isotropic GP
      cov = gp_matern52_cov(x[, 1], sdgp, lscale[1]);
      for (d in 2:Dls) {
        cov = cov .* gp_matern52_cov(x[, d], 1, lscale[d]);
      }
    }
    for (n in 1:N) {
      // deal with numerical non-positive-definiteness
      cov[n, n] += 1e-12;
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
  vector[N] Y;  // response variable
  // data related to GPs
  int<lower=1> Kgp_1;  // number of sub-GPs (equal to 1 unless 'by' was used)
  int<lower=1> Dgp_1;  // GP dimension
  vector[Dgp_1] Xgp_1[N];  // covariates of the GP
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  real Intercept;  // temporary intercept for centered predictors
  vector<lower=0>[Kgp_1] sdgp_1;  // GP standard deviation parameters
  vector<lower=0>[1] lscale_1[Kgp_1];  // GP length-scale parameters
  vector[N] zgp_1;  // latent variables of the GP
  real<lower=0> sigma;  // residual SD
}
transformed parameters {
}
model {
  // initialize linear predictor term
  vector[N] mu = Intercept + rep_vector(0, N) + gp(Xgp_1, sdgp_1[1], lscale_1[1], zgp_1);
  // priors including all constants
  target += student_t_lpdf(Intercept | 3, 7, 10);
  target += student_t_lpdf(sdgp_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += inv_gamma_lpdf(lscale_1[1][1] | 1.149226, 0.019503);
  target += normal_lpdf(zgp_1 | 0, 1);
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
  // to visualize results
  vector[N] mu_sample = Intercept + rep_vector(0, N) + gp(Xgp_1, sdgp_1[1], lscale_1[1], zgp_1);
}
