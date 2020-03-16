# custom stan 
dat <- gamSim(1, n = 100)
m1 <- gam(y ~ s(x2), data = dat)
plot(m1)

library(rstan)
b1 <- stan(file = "customstan.stan", 
           data = list(N = length(dat$x2),
                       x = dat$x2,
                       y = dat$y),
           iter = 400, chains = 2, cores = 2)

plot(b1)

library(brms)
dummy <- brm(y ~ gp(x2), data = dat,
             iter = 400, chains = 2, cores = 2)

#marginal_smooths(dummy)
marginal_effects(dummy) # effects is needed for gaussian process
class(dummy$fit)
class(b1)

dummy2 <- dummy
dummy2$fit <- b1

marginal_smooths(dummy2)
marginal_effects(dummy2)

stancode(dummy)

b1_extract <- extract(b1)

# define a model in terms of the joint distribution of the observed and unobserved 
# data
newdata <- seq(from = 0.01, to = 1, by = 0.01)
dat$newdata <- newdata

b2 <- stan(file="customstan_predict.stan", 
           data=list(x1=dat$x2, y1=dat$y, N1=length(dat$x2), 
                    x2=dat$newdata, N2=length(dat$newdata)),
            iter=400, chains=2, cores = 2)

wot <- extract(b2)
y_new <- colMeans(wot$y2)
plot(y_new~newdata)

# predictions can be made through newdata syntax.

# 2 steps left: 
# custom covariance functions
# multidimensional gaussian processes

# trying out custom covariance function
b3 <- stan(file="customstan_matern32_predict.stan", 
           data=list(x1=dat$x2, y1=dat$y, N1=length(dat$x2), 
                     x2=dat$newdata, N2=length(dat$newdata)),
           iter=400, chains=2, cores = 2) # matern 3/2
m1 <- gam(y ~ s(x2, bs = "gp", m = 2), data = dat) # power exponential
plot(m1)
marginal_effects(dummy)
m3 <- gam(y ~ s(x2, bs = "gp", m = 3), data = dat) # matern 3/2
plot(m3)
# asdf <- predict(m1, newdata = data.frame(x2 = newdata))
# plot(asdf~newdata)


