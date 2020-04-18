// generated with brms 2.10.0
functions {

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
      cov[n, n] += 1e-12;
    }
    return cholesky_decompose(cov) * zgp;
  }

}
data {
  int<lower=1> N;  // number of observations
  vector[N] Y;  // response variable
  // data related to GPs
  // number of sub-GPs (equal to 1 unless 'by' was used)
  int<lower=1> Kgp_1;
  int<lower=1> Dgp_1;  // GP dimension
  // covariates of the GP
  vector[Dgp_1] Xgp_1[N];
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  // temporary intercept for centered predictors
  real Intercept;
  // GP standard deviation parameters
  vector<lower=0>[Kgp_1] sdgp_1;
  // GP length-scale parameters
  vector<lower=0>[1] lscale_1[Kgp_1];
  // latent variables of the GP
  vector[N] zgp_1;
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
  target += normal_lpdf(zgp_1 | 0, 1);
  target += inv_gamma_lpdf(lscale_1[1] | 0.576656, 0.000279);
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
  vector[N] mu_generated = Intercept + rep_vector(0, N) + gp(Xgp_1, sdgp_1[1], lscale_1[1], zgp_1);
}
