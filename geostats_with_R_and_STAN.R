# Geostatistical modeling with R and STAN
# https://aheblog.com/2016/12/07/geostatistical-modelling-with-r-and-stan/

set.seed(1)
evidence <- readRDS("Daten/evidence.csv")
mixed <- dplyr::sample_n(evidence, size = nrow(evidence), replace = FALSE)

testevidence <- head(mixed, n = 500)
# randomly select points where we wish to make predictions
newlocations <- mixed[501:1501,]
newlocations <- dplyr::filter(newlocations, site == 0)
head(testevidence)
head(newlocations)

# create distance matrix
df.loc <- data.frame(lon = c(testevidence$lon, newlocations$lon),
                     lat = c(testevidence$lat, newlocations$lat))

d1 <- as.matrix(dist(df.loc))

# standata
dat <- list(
  N1 = length(testevidence[,1]),
  x1 = testevidence$site,
  N2 = length(newlocations[,1]),
  dist = d1
)

#stanoptions
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 4)

# stanfit
fit <- rstan::stan(
  "geoSTAN.stan",
  data = dat, 
  chains = 4,
  iter = 1000
)
# target variable is missing

## using brms to generate stan code 
library(brms)
set.seed(1)
testevidence <- sample_n(evidence, size = 250, replace = FALSE)
stancode <- brms::make_stancode(formula = "site ~ gp(lon, lat)",
                    data = testevidence, 
                    family = bernoulli(), 
                    chains = 4,
                    cores = 4,
                    iter = 1000,
                    control = list(adapt_delta = 0.8, 
                                   max_treedepth = 13))
# saved as geoSTAN_testevidence.STAN

standata <- brms::make_standata(formula = "site ~ gp(lon, lat)",
                                data = testevidence, 
                                family = bernoulli(), 
                                chains = 4,
                                cores = 4,
                                iter = 1000,
                                control = list(adapt_delta = 0.8, 
                                               max_treedepth = 13))

fit <- rstan::stan(
  file = "geoSTAN_testevidence.stan",
  data = standata,
  cores = 4, 
  chains = 4,
  iter = 1000
)

# this seems to work but takes about 425 seconds for 
# 1000 transitions using 10 leapfrog steps per transition
# lets try with 250 points first

traceplot(fit)

extract(fit)

predictors <- readRDS(file = "Daten/predictors.RDS")
r.pts <- rasterToPoints(predictors, spatial = TRUE)
proj4string(r.pts)

r.pts@coords

r.pts@data <- data.frame(r.pts@data, lon = coordinates(r.pts)[,1],
                         lat = coordinates(r.pts)[,2])

head(r.pts@data)

newdata <- data.frame(lon = r.pts@data$lon, lat = r.pts@data$lat)

test <- brms::posterior_predict(fit, newdata = newdata)

fit2 <- brms::brm(site ~ s(lon, lat),
                  family = bernoulli,
                  data = testevidence,
                  chains = 4, 
                  cores = 4, 
                  iter = 2000,
                  control = list(adapt_delta = 0.8,
                                 max_treedepth = 13)
)

plot(marginal_smooths(fit2))

test <- brms::posterior_predict(fit2, newdata = newdata, nsamples = 10)
averages <- data.frame(predictions = colMeans(test), lon = newdata$lon, lat = newdata$lat)
df_preds <- data.frame(predictions = colMeans(test), lon = newdata$lon, lat = newdata$lat)
colnames(df_preds) <- c("z", "x", "y")
df_preds_new <- data.frame(x = df_preds$x, y = df_preds$y, z = df_preds$z)
coordinates(averages) <- ~lon+lat
# seems to stil resemble bavaria which is a good sign

preds <- SpatialPointsDataFrame(averages@coords, averages@data)
crs(preds) <- crs(predictors)

r <- rasterFromXYZ(df_preds_new)
plot(r)
# thats how we can do predictive mapping with brms objects!

fit3 <- brms::brm(site ~ s(lon, lat , bs="gp"),
                  family = bernoulli,
                  data = testevidence,
                  chains = 4, 
                  cores = 4, 
                  iter = 2000,
                  control = list(adapt_delta = 0.8,
                                 max_treedepth = 13))

plot(marginal_smooths(fit3))

test <- brms::posterior_predict(fit3, newdata = newdata, nsamples = 10)
averages <- data.frame(predictions = colMeans(test), lon = newdata$lon, lat = newdata$lat)
df_preds <- data.frame(predictions = colMeans(test), lon = newdata$lon, lat = newdata$lat)
colnames(df_preds) <- c("z", "x", "y")
df_preds_new <- data.frame(x = df_preds$x, y = df_preds$y, z = df_preds$z)
coordinates(averages) <- ~lon+lat
# seems to stil resemble bavaria which is a good sign

preds <- SpatialPointsDataFrame(averages@coords, averages@data)
crs(preds) <- crs(predictors)

r2 <- rasterFromXYZ(df_preds_new)
plot(r2)

# lets verify if the shortcut over mgcv actually works for real 

data(mcycle, package = 'MASS')
head(mcycle)

m1 <- gam(accel ~ s(times, bs = "gp", m = 1), data = mcycle)
plot(m1)
m2 <- gam(accel ~ s(times, bs = "gp", m = 2), data = mcycle)
plot(m2)
m3 <- gam(accel ~ s(times, bs = "gp", m = 3), data = mcycle)
plot(m3)
m4 <- gam(accel ~ s(times, bs = "gp", m = 4), data = mcycle)
plot(m4)
m5 <- gam(accel ~ s(times, bs = "gp", m = 5), data = mcycle)
plot(m5)

set.seed(1)
b1 <- brms::brm(accel ~ s(times, bs="gp", m = 1), data = mcycle, chains = 4, cores = 4, iter = 1000)
plot(marginal_smooths(b1))
set.seed(1)
b2 <- brms::brm(accel ~ s(times, bs="gp", m = 2), data = mcycle, chains = 4, cores = 4, iter = 1000)
plot(marginal_smooths(b2))
set.seed(1)
b3 <- brms::brm(accel ~ s(times, bs="gp", m = 3), data = mcycle, chains = 4, cores = 4, iter = 1000)
plot(marginal_smooths(b3))
set.seed(1)
b4 <- brms::brm(accel ~ s(times, bs="gp", m = 4), data = mcycle, chains = 4, cores = 4, iter = 1000)
plot(marginal_smooths(b4))
set.seed(1)
b5 <- brms::brm(accel ~ s(times, bs="gp", m = 5), data = mcycle, chains = 4, cores = 4, iter = 1000)
plot(marginal_smooths(b5))
# m = 5 seems to not work

# lets try a 2 dimensional example

dat <- mgcv::gamSim(1, n = 200, scale = 2)
library(mgcViz)
t1 <- gam(y ~ s(x1, x2, bs = "gp", m = 1), data = dat)
plot(getViz(t1))
t4 <- gam(y ~ s(x1, x2, bs = "gp", m = 4), data = dat)
plot(getViz(t4))
set.seed(1)
bt1 <- brms::brm(y ~ s(x1, x2, bs="gp", m = 1), data = dat, chains = 4, cores = 4, iter = 1000)
plot(marginal_smooths(bt1))
set.seed(1)
bt4 <- brms::brm(y ~ s(x1, x2, bs ="gp", m = 4), data = dat, chains = 4, cores = 4, iter = 2000, 
                 control = list(adapt_delta = 0.8,
                                max_treedepth = 13))
plot(marginal_smooths(bt4))

fullmodel_gam <- gam(site ~ s(lon, lat, bs = "gp", m = 3) + dem + temp + rain + 
                       distance_water + frostdays + sunhours + tpi + slope, 
                     family = binomial, 
                     data = evidence)
plot(getViz(fullmodel_gam))

fullmodel_brm <- brm(site ~ s(lon, lat, bs = "gp", m = 3) + dem + temp + rain + 
             distance_water + frostdays + sunhours + tpi + slope, 
           family = bernoulli, 
           data = evidence, 
           chains = 4, 
           cores = 4,
           iter = 2000,
           control = list(adapt_delta = 0.8,
                          max_treedepth = 13))

plot(marginal_smooths(fullmodel_brm))

saveRDS(fullmodel_gam, "fullmodel_gam.RDS")
saveRDS(fullmodel_brm, "fullmodel_brm.RDS")
