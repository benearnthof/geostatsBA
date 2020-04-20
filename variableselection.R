# variableselection with mgcv
evidence <- readRDS("Daten/evidence.csv")
predictors <- stack(c(
  "Daten/dem.grd",
  "Daten/temp_raster.grd",
  "Daten/rain_raster.grd",
  "Daten/water_raster.grd",
  "Daten/frostdays_raster.grd",
  "Daten/sunhours_raster.grd",
  "Daten/tpi_raster.grd",
  "Daten/slope_raster.grd",
  "Daten/aspect_raster.grd"
))
library(mgcv)
fullmodel <- gam(site ~ s(lon, lat, bs = "gp", m = 2) + dem + temp + rain + 
                   distance_water + frostdays + sunhours + tpi + slope + as.factor(aspect),
                 family = binomial, data = evidence, method = "REML")
# all aspect terms are not significant, only aspect = 3 is highly significant. 
# decide to exclude from the model for simplicity
# runtime about 15 seconds on my laptop
summary(fullmodel)
# deviance explained: 31.7%
library(mgcViz)
plot(getViz(fullmodel))

fullmodel <- gam(site ~ s(lon, lat, bs = "gp", m = 2) + dem + temp + 
                   distance_water + frostdays + tpi + slope,
                 family = binomial, data = evidence, method = "REML")

summary(fullmodel)
plot(getViz(fullmodel))





# full smooth terms to see what smooth terms should be included in the model
model2 <- gam(site ~ s(lon, lat, bs = "gp", m = 3) + s(dem) + s(temp) + s(rain) + 
                s(distance_water) + s(frostdays) + s(sunhours) + s(tpi) + s(slope),
              family = binomial, data = evidence)
# runs about 5 minutes on my laptop
summary(model2)
# only slope seems to have one effective degree of freedom
plot(getViz(model2))
# there definitely seems to be a spatial effect
# dem has a slight effect, but very broad confidence bounds past 1000 meters
# the effect of temp is linear past a point where the confidence bounds no longer
# include 0
# effect of rain may also be sufficient as linear term
# distance water seems to have a smooth effect for the most part
# ki of frost days smooth term covers 0 for large chunk
# sunhours should perhaps be included as factor variable
# the effect of tpi is very close to 0 and almost linear
# the effect of slope is linear but significant

# lets try a more parsimonious model
model3 <- gam(site ~ s(lon, lat, bs = "gp", m = 3) + s(dem) + s(temp) + rain + 
                s(distance_water) + frostdays + sunhours + tpi + slope,
              family = binomial, data = evidence)
summary(model3)

# rain and sunhours no longer significant
plot(getViz(model3))
# only distance water seems to have a smooth effect that does not have incredibly
# broad confidence bounds or an almost linear effect. 
# lets try penaliziing the effects of model2

model4 <- gam(site ~ s(lon, lat, bs = "gp", m = 3) + s(dem) + s(temp) + s(rain) + 
                s(distance_water) + s(frostdays) + s(sunhours) + s(tpi) + s(slope),
              family = binomial, data = evidence, select = TRUE)
summary(model4)
plot(getViz(model4))
# rain and distance water should be included as smooth terms
AIC(fullmodel, model2, model3, model4)


# trying out resampling
library(mgcv)
library(caret)

set.seed(0)

dat <- gamSim(1, n = 400, dist = "normal", scale = 2)

b <- train(y ~ x0 + x1 + x2 + x3, 
           data = dat,
           method = "gam",
           trControl = trainControl(method = "cv", number = 1),
           tuneGrid = data.frame(method = "GCV.Cp", select = FALSE)
)

print(b)
summary(b$finalModel)
# s(lon, lat, bs = "gp", m = 2) +
evd <- evidence
evd$site <- as.factor(evd$site)

evd <- evd[sample(nrow(evd), 1000), ]

b <- train(site ~ dem + temp + rain + distance_water + 
             frostdays + sunhours + tpi + slope,
           data = evd,
           method = "gam", 
           trControl = trainControl(method = "cv", number = 5),
           tuneGrid = data.frame(method = "GCV.Cp", select = TRUE)
)
# r + s + r:s
print(b)
summary(b$finalModel)

evd <- evd[sample(nrow(evd), 1000), ]

wot <- gam(site ~ s(lon, lat, bs = "gp", m = 5), data = evidence,
           family = binomial)
plot(getViz(wot))
summary(wot)
defaulttime <- system.time(gam(site ~ s(lon, lat, bs = "gp", m = 3), data = evidence,
                               family = binomial))
remltime <- system.time(gam(site ~ s(lon, lat, bs = "gp", m = 3), data = evidence,
                            family = binomial, method = "REML"))

wot2 <- gam(site ~ s(lon, lat), data = evidence,
                   family = binomial, method = "REML")
plot(getViz(wot2))
wot3 <- gam(site ~ s(lon, lat, bs = "tp"), data = evidence,
            family = binomial)
plot(getViz(wot3))

# Suppose that method = "repeatedcv", number = 10 and repeats = 3,
# then three separate 10-fold cross-validations are used 
# as the resampling scheme.

nrow(evd)

test <- brms::brm(site ~ gp(lon, lat), data = evd, 
                  family = bernoulli,
                  chains = 2, 
                  cores = 2, 
                  iter = 1000, 
                  control = list(adapt_delta = 0.8, 
                                 max_treedepth = 13))

saveRDS(test, file = "2chain1000.RDS")

# starting worker pid=6744 on localhost:11499 at 22:10:05.055
# Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 6690
# Chain 2:  Elapsed Time: 18795.1 seconds (Warm-up)
# Chain 2:                6555.11 seconds (Sampling)
# Chain 2:                25350.2 seconds (Total)

gp_1000_2c <- readRDS("2chain1000.RDS")

plot(gp_1000_2c)

library(rstan)

lon <- seq(from = 9, to = 13.9, length.out = 50)
lat <- seq(from = 47.3, to = 50.6, length.out = 50)

newdata <- expand.grid(lon, lat)
colnames(newdata) <- c("lon", "lat")

preds <- predict(gp_1000_2c, newdata = newdata, nsamples = 10, probs = c(0.5))
head(preds)

newdata$pred <- preds[,1]

library(ggplot2)
ggplot(newdata, aes(x = lon, y = lat, col = pred)) +
  geom_point()

library(raster)
rast <- newdata
colnames(rast) <- c("x", "y", "z")
rast <- rasterFromXYZ(rast, crs="+proj=utm +zone=32 +ellps=WGS84 +units=m +no_defs", digits=5)
plot(rast)
# works better

library(rasterVis)
library(RColorBrewer)
library(viridis)
colr <- colorRampPalette(brewer.pal(11, 'RdYlBu'))
levelplot(rast, 
          margin=FALSE,                       # suppress marginal graphics
          colorkey=list(
            space='bottom',                   # plot legend at bottom
            labels=list(at=seq(0,1,0.1), font=4)      # legend ticks and labels 
          ),    
          par.settings=list(
            axis.line=list(col='transparent') # suppress axes and legend outline
          ),
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=viridis,                   # colour ramp
          at=seq(0, 1, len=101))

# now lets try pure stan predictions:
# step 1: estimate a model in stan
# step 2: reuse stan code with empty model block to just make predictions with 
# set parameters

evidence <- readRDS("Daten/evidence.csv")
set.seed(1)
traindata <- evidence[sample(nrow(evidence), 400),] 
nrow(traindata)

test <- brms::brm(site ~ gp(lon, lat), data = traindata, 
                  family = bernoulli,
                  chains = 2, 
                  cores = 2, 
                  iter = 400, 
                  control = list(adapt_delta = 0.8, 
                                 max_treedepth = 13))
# 10 steps 500 seconds
# 1085.43 seconds (Total)

# generating stancode
stancode <- brms::make_stancode(site ~ gp(lon, lat), data = traindata, 
                                family = bernoulli,
                                chains = 2, 
                                cores = 2, 
                                iter = 400, 
                                control = list(adapt_delta = 0.8, 
                                               max_treedepth = 13))

# parameters of the model 
summary(test)

params <- test$fit@sim$samples
params1 <- params[[1]]
params2 <- params[[2]]

p1 <- as.data.frame.list(params1)
p2 <- as.data.frame.list(params2)

p1 <- p1[, -grep("zgp_", colnames(p1))]
p2 <- p2[, -grep("zgp_", colnames(p2))]

# newdata 
lon <- seq(from = 9, to = 13.9, length.out = 50)
lat <- seq(from = 47.3, to = 50.6, length.out = 50)

newdata <- expand.grid(lon, lat)
colnames(newdata) <- c("lon", "lat")

parameters <- rbind(p1, p2)
parameters <- colMeans(parameters)

standata <- make_standata(site ~ gp(lon, lat), data = traindata, 
                          family = bernoulli,
                          chains = 2, 
                          cores = 2, 
                          iter = 400, 
                          control = list(adapt_delta = 0.8, 
                                         max_treedepth = 13))
standata
nrow(newdata)
names(standata)

head(newdata)
head(standata$Xgp_1)
# modified stancode is in predictions.stan
# we need to append the following variables to the standata
# int<lower=1> N_pred;  // number of observations
# vector[Dgp_1] Xgp_1_pred[N]; // covariates of the GP
str(standata)
nrow(newdata)
standata$N_pred <- as.integer(nrow(newdata))
str(standata)
temp <- standata$Xgp_1
str(temp)
nd <- as.matrix(newdata)
colnames(nd) <- colnames(temp)
standata$Xgp_1_pred <- nd
str(standata)
# should work in theory lets sample that bad boy
# using fit and predict to make predictions
library(rstan)
predict_test <- stan(file = "predictions.stan", 
                     data = standata,
                     iter = 400,
                     chains = 2,
                     cores = 2)

predict_test
wot <- extract(predict_test)
res <- colMeans(wot$y_pred)
newdata$z <- res
# does not work. 
# last resort: simulate gp manually with estimated parameters. ill do that tomorrow.

rast <- newdata
colnames(rast) <- c("x", "y", "z")
rast <- rasterFromXYZ(rast, crs="+proj=utm +zone=32 +ellps=WGS84 +units=m +no_defs", digits=5)
plot(rast)

# well that seems to have gone wrong
# 
gp_1000_2c

nrow(evd)
test <- evd
test$site <- as.factor(test$site)
m_evd <- gam(site ~ s(lon, lat, bs = "gp", m = 2),
             family = binomial,
             data = evd)
summary(m_evd)

plot(m_evd)
p <- predict(m_evd, newdata = newdata, type = "response")

rast <- cbind(newdata, p)
colnames(rast) <- c("x", "y", "z")
rast <- raster::rasterFromXYZ(rast, crs="+proj=utm +zone=32 +ellps=WGS84 +units=m +no_defs", digits=5)
plot(rast)
