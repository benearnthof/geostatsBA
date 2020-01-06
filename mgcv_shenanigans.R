# MGCV shenanigans
# https://people.maths.bris.ac.uk/~sw15190/mgcv/tampere/mgcv.pdf
source("packageloader.R")
?gamSim
set.seed(24)
eg <- gamSim(2, n = 300, scale = 0.05)
b  <- gam(y ~ s(x, z, bs = "gp", k = 50, m = c(3, 0.175)), data = eg$data, method = "REML") ## Matern spline
b1 <- gam(y ~ s(x, z, bs = "gp", k = 50, m = c(1, 0.175)), data = eg$data, method = "REML") ## spherical
b2 <- gam(y ~ s(x, z, bs = "gp", k = 50, m = c(2, 0.175)), data = eg$data, method = "REML") ## exponential

op <- par(mfrow=c(2,2), mar = c(0.5,0.5,3,0.5))
with(eg$truth, persp(x, z, f, theta = 30, main = "Truth")) ## truth
vis.gam(b, theta=30, main = "Matern")
vis.gam(b1, theta=30, main = "Spherical")
vis.gam(b2, theta=30, main = "Exponential")
par(op)

?adaptive.smooth

# https://www.fromthebottomoftheheap.net/2018/04/21/fitting-gams-with-brms/
# https://m-clark.github.io/workshops/stars/extensions.html

# load the example data mcycle
data(mcycle, package = 'MASS')
head(mcycle)

theme_set(theme_bw())
ggplot(mcycle, aes(x = times, y = accel)) +
  geom_point() +
  labs(x = "Miliseconds post impact", y = "Acceleration (g)",
       title = "Simulated Motorcycle Accident",
       subtitle = "Measurements of head acceleration")

# modeling acceleration as a smooth function of time
# default thin plate regression spline basis
m1 <- gam(accel ~ s(times), data = mcycle, method = "REML")
summary(m1)
plot(m1)
draw(m1)

# the same model can be estimated the bayesian way with the brm function of the
# brms package
# cant use te() or ti() tensor smooths, need to use t2() instead.

m2 <- brm(bf(accel ~ s(times)),
          data = mcycle, family = gaussian(), cores = 4, seed = 17,
          iter = 5000, warmup = 1000, thin = 10, refresh = 1,
          control = list(adapt_delta = 0.90))
plot(m2)
summary(m2)
pairs(m2)
shinystan::launch_shinystan(m2)
# sds(stimes_1) is the variance parameter which controls the wiggliness of the
# smooth. the larger this value the more wiggly the smooth.
# The credible interval does not include 0 so there is evidence that a smooth is
# required over and above a linear parametric effect of times.
# stimes_1 is th efixed effect part of the spline which is the linear function
# that is perfectly smooth.

# comparing both models:
gam.vcomp(m1, rescale = FALSE)
# compares well to the credible interval of the brm() version
# (453 & 1153)

# extracting marginal effect of the spline from the model
ms <- marginal_smooths(m2)
plot(ms)
draw(m1)
# they are almost equivalent

# posterior predictive checks to assess model fit
pp_check(m2)
pp_check(m2, type = "ecdf_overlay")
# both plots show deviations between the posterior simulations and the observed
# data. This is because of non constant variance of the acceleration data
# conditional upon the covariate. Both models assumed that the observations
# are distributed gaussian with constant variance.
# https://peerj.com/preprints/27320.pdf

dat <- gamSim(1, n = 200, scale = 2)
fit <- brm(y ~ s(x0) + s(x1) + s(x2) + s(x3), data = dat)
plot(fit)

# all smooth terms
plot(marginal_smooths(fit), rug = TRUE, ask = FALSE)

# fit and plot a two-dimensional smooth term
fit2 <- brm(y ~ t2(x0, x2), data = dat)
ms <- marginal_smooths(fit2)
plot(ms, stype = "contour")
plot(ms, stype = "raster")

# mgcv kovarianzmatrix wird automatisch erzeugt? 
# kann diese an stan ohne weiteres übergeben werden?
# => welche funktioniert mit kriging am besten?
# 
# zeitliche komponente => epochen als kovariable mit aufnehmen
# aufspalten für kleinere nachbarschaftsstrukturen 
# zeitepochen, modellwahl, evtl stan als vergleichsmodell
# kovarianzmatrizen, mgcv schätzverfahren explizit erklären
# penalisierte likelihood
# fahrmeir STAR modelle => strukturiert additive regression
# bayesian theorie 