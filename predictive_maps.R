# predictive maps
model <- readRDS("b_expo_2000_iso.RDS")
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

newdata <- as.data.frame(predictors)
coords <- coordinates(predictors)
nd <- cbind(newdata, coords)
nd <- na.omit(nd)
nrow(nd)
# 120839 points not too much 
rast <- cbind(nd$x, nd$y, nd$dem)
colnames(rast) <- c("x", "y", "z")
rast <- rasterFromXYZ(rast, crs="+proj=utm +zone=32 +ellps=WGS84 +units=m +no_defs", digits=5)
plot(rast)
# transformation worked
# splitting newdata into lists with 25 entries each to avoid crashing R session
names <- colnames(nd)
names[10:11] <- c("lon", "lat")
colnames(nd) <- names
# number of lines per split
n <- 100
nr <- nrow(nd)
# this uses 8 GB of ram because lists seem to have enormous memory requirements
splits <- split(nd, rep(1:ceiling(nr/n), each=n, length.out=nr))
preds <- list()
# length(splits) = 1209
for (i in 1:length(splits)) {
  preds[[i]] <- predict(model, newdata = splits[[i]], nsamples = 50, probs = c(0.5))
  print(i)
}
# the whole loop takes about 1 hour per model for 120000 points
object.size(preds)
# just 3.7mb in ram, seems to be feasible to make predictions for whole raster

df <- do.call(rbind.data.frame, preds)
head(df)

rast <- data.frame(x = nd$lon, y = nd$lat, z = df$Estimate)
colnames(rast) <- c("x", "y", "z")
library(raster)
rast <- rasterFromXYZ(rast, crs="+proj=utm +zone=32 +ellps=WGS84 +units=m +no_defs", digits=5)
plot(rast)

saveRDS(rast, file = "b_expo_2000_iso_predictions50.RDS")

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
