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