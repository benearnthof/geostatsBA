//function block must appear before all other blocks
//function calling does not work because of mistyped variables
//lets try calculating the matrix directly in the transformed data block
/*
 functions {
  matrix gp_matern32_cov(real[] x, real alpha, real rho) {
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
*/
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
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
  vector[N] eta;
}

transformed parameters{
    vector[N] f; //The log-intensity process
    {
    matrix[N, N] L_K;
    matrix[N, N] K; // = gp_matern32_cov(x, alpha, rho);
    for (i in 1:N) {
        K[i,i] = alpha*alpha;
        for (j in 1:(i-1)) {
          // dt = distance between 2 measurements
            K[i,j] = alpha*alpha * (1 + sqrt(3)*dt*(i-j)/rho) * exp(-sqrt(3)*dt*(i-j)/rho);
            K[j,i] = K[i,j];
        }
    }
    // diagonal elements
    for (n in 1:N)
      K[n, n] = K[n, n] + delta;

    L_K = cholesky_decompose(K);
    f = L_K * eta;
  }

    

    x = mu + cholesky_decompose(Sigma)*x_std;
}
transformed parameters {
  vector[N] f;
  {
    matrix[N, N] L_K;
    matrix[N, N] K = gp_matern32_cov(x, alpha, rho);

    // diagonal elements
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
