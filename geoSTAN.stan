data {
  int <lower = 1> N1; //number of data points
  int x1[N1]; //site yes no
  // int n1[N1]; //number of observations
  int <lower = 1> N2;//number of new points
  matrix[N1+N2, N1+N2] dist; //distances between points
}

transformed data {
  int <lower = 1> N;
  N = N1 + N2;
}

parameters {
  vector[N1] y1;
  vector[N2] y2;
  real beta;
  real sigma_sq;
  real phi;
}

transformed parameters {
  vector [N1+N2] mu;
  for(i in 1:N) mu[i] = beta;
}

model {
  vector[N] y;
  matrix[N, N] Sigma; 
  matrix[N, N] L;
  for(i in 1:N1) y[i] = y1[i];
  for(i in 1:N2) y[N1+i] = y2[i];
  // covariance matrix
  for(i in 1:(N-1)){
   for(j in (i+1):N){
     Sigma[i,j] = exp((-1)*phi*dist[i,j]);
     Sigma[j,i] = Sigma[i,j];
  }
 }
 for(i in 1:N) Sigma[i,i] = sigma_sq;
 // using cholesky decomposition to speed up sampling
 L = cholesky_decompose(Sigma);
 sigma_sq ~ normal(0, 5);
 phi ~ normal(0, 5);
 y ~ multi_normal_cholesky(mu,L);
 beta ~ normal(0,5);
 x1 ~ bernoulli(y1);
}

generated quantities {
  vector[N2] y_pred;
  // predicted values at new points
  for(i in 1:N2) y_pred[i] = inv_logit(beta+y2[i]);
}
