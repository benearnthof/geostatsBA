# icar example code
library(mgcv)
library(brms)

# generate spatial data 
east <- north <- 1:10
Grid <- expand.grid(east, north)
K <- nrow(Grid)

# set up distance and neighbourhood matrices
distance <- as.matrix(dist(Grid))
W <- array(0, c(K, K))
W[distance == 1] <- 1

# generate the covariates and response data
x1 <- rnorm(K)
x2 <- rnorm(K)
theta <- rnorm(K, sd = 0.05)
phi <- rmulti_normal(
  1, mu = rep(0, K), Sigma = 0.4 * exp(-0.1 * distance)
)
eta <- x1 + x2 + phi
prob <- exp(eta) / (1 + exp(eta))
size <- rep(50, K)
y <- rbinom(n = K, size = size, prob = prob)

dat <- data.frame(y, size, x1, x2)

# fit a CAR model to the data 
fit <- brm(y | trials(size) ~ x1 + x2, data = dat, 
           family = binomial(), autocor = cor_car(W),
           chains = 4, cores = 4, iter = 5000, warmup = 2000)

plot(marginal_effects(fit))

fitted_values <- fitted(fit)
data <- as.data.frame(cbind(Y = standata(fit)$Y, fitted_values))
library(ggplot2)
ggplot(data) +
  geom_point(aes(x = Estimate, y = Y))

fit2 <- brm(y | trials(size) ~ x1 + x2, data = dat, 
           family = binomial(),
           chains = 4, cores = 4, iter = 5000, warmup = 2000)

data2 <- as.data.frame(cbind(Y = standata(fit2)$Y, fitted(fit2)))
ggplot(data2) +
  geom_point(aes(x = Estimate, y = Y))

######
make_stancode(rating ~ treat + period + carry + (1|subject), 
              data = inhaler, family = "cumulative")

library(bayesplot)
bayesplot::mcmc_plot(fit)
plot(fit)
