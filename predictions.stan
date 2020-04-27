functions {
  matrix gp(vector[] x, real sdgp, real lscale, int D) { 
    int N = size(x);
    matrix[N, N] cov;
    if (D == 1) {
      // one dimensional or isotropic GP
      cov = cov_exp_quad(x, sdgp, lscale);
    } else {
      // multi-dimensional GP
      cov = cov_exp_quad(x[, 1], sdgp, lscale);
      for (d in 2:D) {
        cov = cov .* cov_exp_quad(x[, d], 1, lscale);
      }
    }
    for (n in 1:N) {
      // deal with numerical non-positive-definiteness
      cov[n, n] += 1e-8;
    }
    return cov;
  }
}

data {
  int<lower=1> N;  // number of observations
  // number of sub-GPs (equal to 1 unless 'by' was used)
  int<lower=1> Kgp_1;
  int<lower=1> Dgp_1;  // GP dimension
  real<lower=0> sdgp_1;
  real<lower=0> lscale_1;
  real Intercept;
  // covariates of the GP
  vector[Dgp_1] Xgp_1[N];
}

transformed data {
  matrix[N, N] cov = gp(Xgp_1, sdgp_1, lscale_1, Dgp_1);
  matrix[N, N] L_cov = cholesky_decompose(cov);
}

parameters {}
model {}

generated quantities {
  vector[N] f = multi_normal_cholesky_rng(rep_vector(Intercept, N), L_cov);
  int y[N];
  for (n in 1:N) {
    y[n] = bernoulli_logit_rng(f[n]);
    // y[n] = poisson_log_rng(f[n]); 
  }
   
}
