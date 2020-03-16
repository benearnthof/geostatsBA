// latent variable gaussian process for non normal outcome
// add a small number to diagonal of covariance matrix to ensure positive definiteness
data {
  int<lower=1> N;
  real x[N];
  vector[N] y;
}
transformed data {
  real delta = 1e-9;
}

functions {
  matrix gp_matern32_cov(vector[] x, real alpha, real rho) {
    real sigma_sq = square(alpha);
    real root_3 = sqrt(3.0);
    matrix [size(x), size(x)] cov;
    vector[size(x)] x_new;

    for (i in 1:size(x)) {
       for (j in 1:size(x)) {
         real dist = distance(x[i]/rho, x[j]/rho);
        cov[i, j] = sigma_sq * (1.0 + root_3 * dist) * exp(-root_3 * dist);
        cov[j, i] = cov[i, j];
      }
    }
    return cov;
  }
}

parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
  vector[N] eta;
}
model {
  vector[N] f;
  {
    matrix[N, N] L_K;
    // standard cov_exo_quad
    // brms equivalent: cov_exp_quad(x, sdgp, lscale[1]);
    matrix[N, N] K = gp_matern32_cov(x, alpha, rho);

    // diagonal elements
    for (n in 1:N)
      K[n, n] = K[n, n] + delta;

    L_K = cholesky_decompose(K);
    f = L_K * eta;
  }

  rho ~ inv_gamma(5, 5);
  alpha ~ std_normal();
  sigma ~ std_normal();
  eta ~ std_normal();

  y ~ normal(f, sigma);
}
