data {
  int K;
  int Y[K];
  int N[K];
}

parameters{
  real<lower=0, upper=1> alpha;
  real<lower=1> kappa;
  vector<lower=0, upper=1>[K] theta;
}

model {
  for (i in 1:K) {
    Y[i] ~ binomial(N[i], theta[i]); // binomial likelihood
  }
  // prior
  theta ~ beta(alpha * kappa, (1 - alpha) * kappa);
  
  // hyper priors
  kappa ~ pareto(1, 0.3);
  alpha ~ beta(5, 5);
}

generated quantities{
  real<lower=0, upper=1> aTheta;
  aTheta = beta_rng(alpha * kappa, (1 - alpha) * kappa);
}
