require(rstan)
K <- 30 # number of polls
aN <- 10 # number of individuals in each poll
theta <- 0.35
Y <- rbinom(K, aN, theta)
N <- rep(aN, K)

aModel <- stan_model("hierarchical_eu.stan")
fit <- sampling(aModel, data = list(N=N, K=K, Y=Y), iter = 400, chains = 4)

require(shinystan)
pairs(fit)
shinystan::launch_shinystan(fit)
