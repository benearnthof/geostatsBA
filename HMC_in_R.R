# Hamiltonian Monte Carlo in R
# https://jonnylaw.netlify.com/2019/07/31/hamiltonian_monte_carlo_in_r/
# https://blogs.rstudio.com/tensorflow/posts/2019-10-03-intro-to-hmc/
# Determining th eposterior distribution for the parameters of a real-world 
# Bayesian model inevitably requires calculating high-dimensional integrals. 
# Often these are tedious or impossible to calculate by hand. MCMC algorithms 
# are popular approoaches, samplers such as the Gibbs sampler can be used to 
# sample from models with conditionally conjugate specifications and the MH-
# algorithm can be used when the conditionally conjugate form is not present. 
# There are downsides to these established methods however: Gibbs sampling 
# is limited to a restrictive form of the prior distribution and, although MH 
# allows any prior distribution from which the probability density function can 
# be evaluated, the proposal distribution for MH is often a multivariate normal 
# distribution with mean corresponding to the previous value of the parameters.
# This form of proposals results in random walk behaviour, which in turn causes
# the MH algorithm to converge very slowly, thus taking a very long time to
# yield satisfactory draws from the posterior distribution. 

# The gradient of the un-normalized log-posterior distribution can be used to
# explore the posterior distribution in a more efficient manner. Hamiltonian 
# Monte Carlo (HMC) is a MCMC method that utilises a discretisation of 
# Hamilton's equations in order to model a physical system where the parameters 
# are represented by the position of a particle in the parameter space. 
# In order to implement HMC, the posterior distribution is augmented with a 
# momentum vector Phi, which is used to propose updates to the position of this
# massless particle. These proposed positions can be very far away from the 
# initial position, allowing us to traverse a wide array of locations in 
# posterior space even in a limited number of steps. 

# HMC in R
# Hamilton's equations are discretised via "leapfrog-integration". These 
# leapfrog steps can be written in the following manner: 

# one leapfrog step
leapfrog_step <- function(gradient, step_size, position, momentum, d) {
  momentum1 <- momentum + gradient(position) * 0.5 * step_size
  position1 <- position + step_size * momentum1
  momentum2 <- momentum1 + gradient(position1) * 0.5 * step_size
  matrix(c(position1, momentum2), ncol = d * 2)
}

# multiple leapfrog steps
# l = number of steps
leapfrogs <- function(gradient, step_size, l, position, momentum, d) {
  for (i in 1:l) {
    pos_mom <- leapfrog_step(gradient, step_size, position, momentum, d)
    position <- pos_mom[seq_len(d)]
    momentum <- pos_mom[-seq_len(d)]
  }
  pos_mom
}

# the log-acceptance can be written as: 
log_acceptance <- function(propPosition, propMomentum, position, momentum, log_posterior) {
  log_posterior(propPosition) + sum(dnorm(propMomentum, log = T)) - log_posterior(position) - sum(dnorm(momentum, log = T))
}

# in order to propose a new set of parameters, a random momentum vector is drawn
# from the normal distribution and used in the leapfrog steps. To ensure 
# detailed balance and that the stationary distribution of the markov chain is 
# equivalent to the target distribution, a metropolis step is used accepting the
# newly proposed parameters with log-probability equal to the log_acceptance 
# defined above. To implement this in R, we compare a log uniform random number 
# to the log acceptance criterion. 

# A single, complete HMC step can be written as: 

hmc_step <- function(log_posterior, gradient, step_size, l, position) {
  d <- length(position)
  momentum <- rnorm(d)
  pos_mom <- leapfrogs(gradient, step_size, l, position, momentum, d)
  propPosition <- pos_mom[seq_len(d)]
  propMomentum <- pos_mom[-seq_len(d)]
  a <- log_acceptance(propPosition, propMomentum, position, momentum, log_posterior)
  if (log(runif(1)) < a) {
    propPosition
  } else {
    position
  }
}

# The complete HMC algorithm can be written as 
hmc <- function(log_posterior, gradient, step_size, l, initP, m) {
  out <- matrix(NA_real_, nrow = m, ncol = length(initP))
  out[1,] <- initP
  for (i in 2:m) {
    out[i, ] <- hmc_step(log_posterior, gradient, step_size, l, out[i - 1, ])
  }
  out
}

# example model: Bivariate Normal Model
# theory added in post

# the gradient can be written like: 
gradient <- function(ys) {
  function(theta) {
    mu <- c(theta[1], theta[3])
    sigma <- c(theta[2], theta[4])
    n <- nrow(ys)
    c(1/sigma[1]^2*sum(ys[,1] - mu[1]) - mu[1]/9,
      -n/sigma[1] + sum((ys[,1] - mu[1])^2) / sigma[1]^3 + 2/sigma[1] - 3,
      1/sigma[2]^2*sum(ys[,2] - mu[2]) - mu[2]/9,
      -n/sigma[2] + sum((ys[,2] - mu[2])^2) / sigma[2]^3 + 2/sigma[2] - 3)
  }
}

# HMC works best when the leapfrog proposal can propose unconstrained values of 
# the parameters which lie on the real number line. 
# We need a transform function for the parameters, which calculates the 
# exponential of the standard deviation parameters. The log-posterior is 
# calculated using the transformed values, the appropriate transformation and 
# inverse transformation can be written as: 

transform <- function(theta) {
  c(theta[1], exp(theta[2]), theta[3], exp(theta[4]))
}

inv_transform <- function(theta) {
  c(theta[1], log(theta[2]), theta[3], log(theta[4]))
}

# The leapfrog step proposal is calculated using the unconstrained parameters, 
# hence the derivative of the log-jacobian of the transformation is required to 
# be added to the value of the gradient of the log-density. Then the derivative 
# of the log-jacobian is calculated to get the value of the gradient 
# corresponding to the unconstrained parameters int he leapfrog step. 
log_jacobian <- function(theta) {
  c(0, theta[2], 0, theta[4])
}

deriv_log_jacobian <- function(theta) {
  c(0, 1, 0, 1)
}

# In our case the derivative of the log-jacobian contributes the value 1 to each 
# of the partial derivatives. 