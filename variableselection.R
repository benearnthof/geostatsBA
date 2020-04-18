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
fullmodel <- gam(site ~ s(lon, lat, bs = "gp", m = 3) + dem + temp + rain + 
                   distance_water + frostdays + sunhours + tpi + slope + as.factor(aspect),
                 family = binomial, data = evidence)
# all aspect terms are not significant, only aspect = 3 is highly significant. 
# decide to exclude from the model for simplicity
# runtime about 15 seconds on my laptop
summary(fullmodel)
# deviance explained: 31.7%
library(mgcViz)
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
