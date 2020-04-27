evidence <- readRDS("Daten/evidence.csv")

library(brms)
library(mgcv)

set.seed(1)
evd <- evidence[sample(nrow(evidence), 500),]
nrow(evd)

gp_500 <- brms::brm(site ~ gp(lon, lat), data = evd, 
                     family = bernoulli,
                     chains = 4, 
                     cores = 4, 
                     iter = 1000, 
                     control = list(adapt_delta = 0.8, 
                                    max_treedepth = 13))
saveRDS(gp_500, "gp_500.RDS")
# longest gradient evaluation for 1000 steps: 526.45 seconds
# longest sampling time for single chain
# Chain 2:  Elapsed Time: 1238.15 seconds (Warm-up)
# Chain 2:                572.664 seconds (Sampling)
# Chain 2:                1810.81 seconds (Total)

lon <- seq(from = 9, to = 13.9, length.out = 10)
lat <- seq(from = 47.3, to = 50.6, length.out = 10)

newdata <- expand.grid(lon, lat)
colnames(newdata) <- c("lon", "lat")

preds <- predict(gp_500, newdata = newdata, nsamples = 10, probs = c(0.5))
head(preds)

newdata$pred <- preds[,1]

library(ggplot2)
ggplot(newdata, aes(x = lon, y = lat, col = pred)) +
  geom_point()

library(raster)

rast <- newdata
colnames(rast) <- c("x", "y", "z")
rast <- rasterFromXYZ(rast, crs="+proj=utm +zone=32 +ellps=WGS84 +units=m +no_defs", digits=5)
# works better

# lets try splitting up predictions into multiple steps

lon <- seq(from = 9, to = 13.9, length.out = 100)
lat <- seq(from = 47.3, to = 50.6, length.out = 100)

newdata <- expand.grid(lon, lat)
colnames(newdata) <- c("lon", "lat")

n <- 25
nr <- nrow(newdata)
splits <- split(newdata, rep(1:ceiling(nr/n), each=n, length.out=nr))
preds <- list()

for (i in 1:length(splits)) {
  preds[[i]] <- predict(gp_500, newdata = splits[[i]], nsamples = 50, probs = c(0.5))
  print(i)
}
object.size(preds)
# just 500kb in ram, seems to be feasible to make predictions for whole raster

df <- do.call(rbind.data.frame, preds)
head(newdata)
head(df)
newdata$z <- df$Estimate

rast <- newdata
colnames(rast) <- c("x", "y", "z")
library(raster)
rast <- rasterFromXYZ(rast, crs="+proj=utm +zone=32 +ellps=WGS84 +units=m +no_defs", digits=5)
plot(rast)

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
          col.regions=colr,                   # colour ramp
          at=seq(0, 1, len=101))

# splitting the raster to a list and then predicting individually seems to work fine.

# lets try and parallelize this shit

nrow(newdata)
lon <- seq(from = 9, to = 13.9, length.out = 4)
lat <- seq(from = 47.3, to = 50.6, length.out = 4)

nd <- expand.grid(lon, lat)
colnames(nd) <- c("lon", "lat")

get_brms_predictions <- function(nd, nsamples = 25, splitsize = 25, cores = 4,
                                 model, prob = 0.5) {
  n <- splitsize
  nr <- nrow(nd)
  splits <- split(nd, rep(1:ceiling(nr/n), each = n, length.out = nr))
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  foreach::foreach(i = seq_along(splits)) %dopar% {
    newd <- splits[[1]]
    preds <- predict(model, newdata = newd, nsamples = nsamples, probs = c(prob))
    sapply(gp_500, predict, newdata = newd, nsamples = 10, probs = 0.5)
    predict(model, newdata = newd, nsamples = 10, probs = 0.5)
    purrr::map(gp_500, brms::predict.brmsfit, newdata = newd, nsamples = 10, probs = 0.5)
    newdata$pred <- preds$Estimate
  }
  
  
  parallel::stopCluster(cl)
}

# lets try and set knots for approximate gaussian processes
evidence <- readRDS("Daten/evidence.csv")
library(brms)
library(mgcv)
set.seed(1)
evd <- evidence[sample(nrow(evidence), 2000),]
nrow(evd)

gp_2000_aprox30 <- brms::brm(site ~ gp(lon, lat, k = 30, c = 5/4), data = evd, 
                    family = bernoulli,
                    chains = 4, 
                    cores = 4, 
                    iter = 1000, 
                    control = list(adapt_delta = 0.8, 
                                   max_treedepth = 12))
# Chain 4: 1000 transitions using 10 leapfrog steps per transition would 
# take 241.53 seconds. 
# 2000 points with approximate gp seems to be more than twice as fast as 
# doing 500 points and exact gp. 
# chain 2 fastest at 3017 seconds total
# chain 3 slowest at 4256 seconds total

saveRDS(gp_2000_aprox30, "gp_2000_aprox.RDS")

# there are 900 latent variable terms in the model fit
# 30 by 30 => approximation is a lot better with large data sets

# trying to predict a raster of the original dimensions
lon <- seq(from = 9, to = 13.9, length.out = 396)
lat <- seq(from = 47.3, to = 50.6, length.out = 584)

newdata <- expand.grid(lon, lat)
colnames(newdata) <- c("lon", "lat")

n <- 25
nr <- nrow(newdata)
splits <- split(newdata, rep(1:ceiling(nr/n), each=n, length.out=nr))
preds <- list()

# 9251
for (i in 1:length(splits)) {
  preds[[i]] <- predict(gp_500_aprox30, newdata = splits[[i]], nsamples = 25, probs = c(0.5))
  print(i)
}
# because we limit the amount of latent terms not only is modeling faster, 
# predicting is also quite a lot faster. 

df <- do.call(rbind.data.frame, preds)
head(newdata)
head(df)
newdata$z <- df$Estimate

rast <- newdata
colnames(rast) <- c("x", "y", "z")
library(raster)
rast <- rasterFromXYZ(rast, crs="+proj=utm +zone=32 +ellps=WGS84 +units=m +no_defs", digits=5)
par(mar = c(1, 1, 1, 1))
plot(rast)

saveRDS(rast, file = "hires_pred_raster1.RDS")

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
