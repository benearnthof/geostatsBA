//lets try swapping out the covariance function for an equivalent 
data {
  int<lower=1> N1;
  real x1[N1];
  vector[N1] y1;
  int<lower=1> N2;
  real x2[N2];
}
transformed data {
  real delta = 1e-9;
  int<lower=1> N = N1 + N2;
  real x[N];
  for (n1 in 1:N1) x[n1] = x1[n1];
  for (n2 in 1:N2) x[N1 + n2] = x2[n2];
}
parameters {
  real<lower=0> rho; //scale parameter
  real<lower=0> alpha; //sigma parameter for diagonals
  real<lower=0> sigma;
  vector[N] eta;
}
transformed parameters {
  vector[N] f;
  {
    matrix[N, N] L_K;
    //matrix[N, N] K = cov_exp_quad(x, alpha, rho);
    matrix[N, N] K;
    /*using squared distances like it is implemented in the stan math library
    cov.diagonal().array() = sigma_sq;
    for (size_t j = 0; j < x_size; ++j) {
      for (size_t i = j + 1; i < x_size; ++i) {
        cov(i, j)
          = sigma_sq * exp(squared_distance(x[i], x[j]) * neg_half_inv_l_sq);
      }
    }
    */ 
    //diagonal elements
    for (i in 1:N) {
      K[i, i] = alpha * alpha;
    }
    for (i in 1:N) {
      for (j in (i+1):N) {
        real sq_alpha = alpha * alpha; //square alpha
        real sq_dist = (x[i] - x[j]) * (x[i] - x[j]); //squared euclidian dist
        K[i,j] = sq_alpha * exp(-(sq_dist/(rho * rho))); //squared exponential cov
        K[j,i] = K[i, j]; //symmetric matrix
      }
    }
    
    // making sure matrix is positive definite for cholesky
    for (n in 1:N)
      K[n, n] = K[n, n] + delta;

    L_K = cholesky_decompose(K);
    f = L_K * eta;
  }
}
model {
  rho ~ inv_gamma(5, 5);
  alpha ~ std_normal();
  sigma ~ std_normal();
  eta ~ std_normal();

  y1 ~ normal(f[1:N1], sigma);
}
generated quantities {
  vector[N2] y2;
  for (n2 in 1:N2)
    y2[n2] = normal_rng(f[N1 + n2], sigma);
}
