// generated with brms 2.10.0
functions {
  
  /* matrix cov_custom(vector x, real sdgp, real lscale) {
    int N = cols(x);
    matrix[N, N] cov;
    //diagonal elements
    for (i in 1:N) {
      cov[i, i] = sdgp * sdgp;
    }
    for (i in 1:N) {
      for (j in (i+1):N) {
        real sq_sigma = sdgp * sdgp; //square alpha
        real sq_dist = (x[i] - x[j]) * (x[i] - x[j]); //squared euclidian dist
        cov[i,j] = sq_sigma * exp(-(sq_dist/(lscale * lscale))); //squared exponential cov
        cov[j,i] = cov[i, j]; //symmetric matrix
      }
    }
    return(cov);
  } */ 
   matrix gp_matern52_cov(vector x, real sigma, real length_scale) {
    matrix[cols(x), cols(x)] cov;

    real sigma_sq = sigma^2;
    real root_5_inv_l = sqrt(5.0) / length_scale;
    real inv_l_sq_5_3 = 5.0 / (3.0 * length_scale^2);

    for (i in 1:cols(x)) {
      cov[i, i] = sigma_sq;
        for (j in (i + 1):cols(x)) {
        real sq_distance = (x[i] - x[j])^2;
        real dist = sqrt(sq_distance);
        cov[i, j] = sigma_sq
                  * (1.0 + root_5_inv_l * dist + inv_l_sq_5_3 * sq_distance)
                  * exp(-root_5_inv_l * dist);
        cov[j, i] = cov[i, j];
      }
    }
    return(cov);
  }  
  
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
      // cov = cov_exp_quad(x, sdgp, lscale[1]);
      // cov = cov_custom(to_vector(x[, 1]) , sdgp, lscale[1]);
      // cannot hand over an array here needs to be vector
      cov = gp_matern52_cov(x, sdgp, lscale[1]);
    } else {
      // multi-dimensional non-isotropic GP
      // cov = cov_exp_quad(x[, 1], sdgp, lscale[1]);
      // cov = cov_custom(to_vector(x[, 1]), sdgp, lscale[1]);
      cov = gp_matern52_cov(x[, 1], sdgp, lscale[1]);
      for (d in 2:Dls) {
        // cov = cov .* cov_exp_quad(x[, d], 1, lscale[d]);
        // cov = cov .* cov_custom(to_vector(x[, d]), 1, lscale[d]);
        cov = cov .* gp_matern52_cov(x[, d], sdgp, lscale[d]);
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
  vector<lower=0>[1] lscale_1[1];
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
  target += inv_gamma_lpdf(lscale_1[1] | 1.089792, 0.015279);
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
