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
pairs(fit)

# trying gp models
dat <- mgcv::gamSim(1, n = 300, scale = 2)

fit1 <- brm(y ~ s(x1, x2, bs = "gp"), dat, chains = 4, cores = 4)

plot(fit1)
me <- marginal_effects(fit1)
plot(me)
ms <- marginal_smooths(fit1)
plot(ms)

fit2 <- gam(y ~ s(x1, x2, bs = "gp"), data = dat)
# the library schoenberg has been taken over by a package for 12 tone music
# the original one has been renamed to "gratia"
# library(schoenberg)
# install.packages('gratia')
remotes::install_github("gavinsimpson/gratia")
draw(fit2)
plot(fit2)

dat <- mgcv::gamSim(1, n = 200, scale = 2)
fit1 <- brm(y ~ s(x1, x2, bs = "gp"), dat, chains = 4, cores = 4)
fit2 <- brm(y ~ gp(x1, x2), dat, chains = 4, cores = 4)
fit3 <- gam(y ~ s(x1, x2, bs = "gp"), data = dat)
fit4 <- brm(y ~ t2(x1, x2), data = dat)

fits <- list(fit1, fit2, fit3, fit4)
saveRDS(fits, file = "testfits.RDS")

ms1 <- marginal_smooths(fit1)
ms2 <- marginal_smooths(fit2)
ms2 <- marginal_effects(fit2)
ms3 <- marginal_effects(fit3)
ms4 <- marginal_smooths(fit4)

plot(ms1)
plot(ms2)
plot(fit3)
plot(ms4)

plot(fit3)
plot(fit2)

# cor_sar 

data(oldcol, package = "spdep")

fit0 <- brm(CRIME ~ INC + HOVAL, data = COL.OLD,
            chains = 2, cores = 2)

fit1 <- brm(CRIME ~ INC + HOVAL, data = COL.OLD, 
            autocor = cor_lagsar(COL.nb), 
            chains = 2, cores = 2)
summary(fit1)
plot(fit1)

fit2 <- brm(CRIME ~ INC + HOVAL, data = COL.OLD, 
            autocor = cor_errorsar(COL.nb), 
            chains = 2, cores = 2)
summary(fit2)
plot(fit2)

plt_fitted <- function(fit) {
  fitted_values <- fitted(fit)
  data <- as.data.frame(cbind(Y = standata(fit)$Y, fitted_values))
  plt <- ggplot(data) +
    geom_point(aes(x = Estimate, y = Y))
  return(plt)
}

plt_fitted(fit0)
plt_fitted(fit1)
plt_fitted(fit2)
