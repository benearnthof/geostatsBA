# I'm only interested in settlements from the iron age
# get raw data and predictorstack from preprocessing.R
predictors <- readRDS(file = "Daten/predictors.RDS")
raster::unique(raw_data$Epoche)
presence <- raw_data[raw_data$Epoche %in% c("Hallstattzeit", "Latènezeit"),]
nrow(presence)
sites <- dplyr::select(presence, lng_wgs84, lat_wgs84)
head(sites)
# selecting unique rows to combat biased models
sites <- unique(sites[,c("lng_wgs84", "lat_wgs84")])
# 6206 sites remaining
# creating buffers around points to sample nonsites correctly

temp <- sampleRandom(predictors[[1]], 1000, sp = T)
sp_sites <- sp::SpatialPoints(coords      = sites[,c("lng_wgs84","lat_wgs84")], # order matters
                              proj4string = predictors@crs)
utm_sites <- spTransform(sp_sites, CRS("+proj=utm +zone=32 +ellps=WGS84"))
sf_sites <- st_as_sf(utm_sites)

sites_buff_1000m <- st_buffer(sf_sites, dist = 1000)
temp_utm <- spTransform(temp, CRS("+proj=utm +zone=32 +ellps=WGS84"))

# tmp <- temp_utm[sites_buff_1000m$geometry,]

sp_polygons_buffer <- sf::as_Spatial(sites_buff_1000m$geometry)
crs(sp_polygons_buffer) <- crs(temp_utm)
over(temp_utm, sp_polygons_buffer)
nrow(temp_utm)
ret <- temp_utm[!is.na(over(temp_utm, sp_polygons_buffer)),]
# this returns the points that fall within the buffer zone. 
# lets try to increase the radius
sites_buff_2500m <- st_buffer(sf_sites, dist = 2500)
sp_polygons_2500 <- sf::as_Spatial(sites_buff_2500m$geometry)
crs(sp_polygons_2500) <- crs(temp_utm)
ret2500_within <- temp_utm[!is.na(over(temp_utm, sp_polygons_2500)),]
nrow(ret2500_within)
ret2500_within <- spTransform(ret2500_within, crs(temp))
ret2500_without <- temp_utm[is.na(over(temp_utm, sp_polygons_2500)),]
ret2500_without <- spTransform(ret2500_without, crs(temp))

mapview(sp_sites) +
  mapview(ret2500_without, color = "red", fill = "red")
# works as expected

# function to sample around the points with a given buffer
buffsample <- function(ssize = 1000, distance = 1000, within = FALSE, returnsize = 1000) {
  ret <- new("SpatialPoints",                                                          
             coords = structure(numeric(0), .Dim = c(0L, 2L),                          
                                .Dimnames = list(NULL, c("coords.x1", "coords.x2"))),  
             bbox = structure(c(1, 1, 1, 1), .Dim = c(2L, 2L),                         
                              .Dimnames = list(c("coords.x1", "coords.x2"),
                                               c("min", "max"))),
             proj4string = new("CRS", projargs = "+proj=utm +zone=32 +ellps=WGS84"))
  sp_sites <- sp::SpatialPoints(coords = sites[,c("lng_wgs84","lat_wgs84")], proj4string = predictors@crs)
  utm_sites <- spTransform(sp_sites, CRS("+proj=utm +zone=32 +ellps=WGS84"))
  sf_sites <- st_as_sf(utm_sites)
  sites_buff <- st_buffer(sf_sites, dist = distance)
  
  sp_polygons_buffer <- sf::as_Spatial(sites_buff$geometry)
  while (nrow(ret@coords) < returnsize) {
    sample <- sampleRandom(predictors[[1]], ssize, sp = T)
    sample_utm <- spTransform(sample, CRS("+proj=utm +zone=32 +ellps=WGS84"))
    crs(sp_polygons_buffer) <- crs(sample_utm)
    over(sample_utm, sp_polygons_buffer)
    nrow(sample_utm)
    if (within == TRUE) {
      tmp <- sample_utm[!is.na(over(sample_utm, sp_polygons_buffer)),]
    } else {
      tmp <- sample_utm[is.na(over(sample_utm, sp_polygons_buffer)),]
    }
    ret <- maptools::spRbind(ret, tmp)
  }
  ret <- spTransform(ret, crs(sp_sites))
  ret <- ret[1:returnsize,]
  ret <- coordinates(ret)
  colnames(ret) <- c("lng_wgs84", "lat_wgs84")
  ret <- as.data.frame.matrix(ret)
  return(ret)
}
set.seed(123)
funtest <- buffsample(ssize = 2000, distance = 1500, returnsize = 10000)
funtest5k <- buffsample(ssize = 2000, distance = 1500, returnsize = 5000)
sites_buff_1500m <- st_buffer(sf_sites, dist = 1500)
nrow(funtest)

coordinates(funtest5k) <- c("lng_wgs84", "lat_wgs84")
crs(funtest5k) <- CRS("+proj=utm +zone=32 +ellps=WGS84 +units=m +no_defs")
mapview(sites_buff_1500m) + 
  mapview(funtest5k, alpha.regions = 0.1, color = "red")
# works as intended and performs reasonably well
# finalizing evidence

# function to automate generation of evidence data
generateEvidence <- function(sitesdata, nonsitesdata, predictorstack = predictors) {
  # selecting site points
  sites_temp <- sitesdata
  sites_temp$lon <- as.numeric(as.vector(sites_temp$lng_wgs84))
  sites_temp$lat <- as.numeric(as.vector(sites_temp$lat_wgs84))
  # convert to spatial data in order to extract the predictor values for all points
  coordinates(sites_temp) <- c("lon","lat")
  proj4string(sites_temp) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
  nonsites_temp <- nonsitesdata
  nonsites_temp$lon <- as.numeric(as.vector(nonsites_temp$lng_wgs84))
  nonsites_temp$lat <- as.numeric(as.vector(nonsites_temp$lat_wgs84))
  coordinates(nonsites_temp) <- c("lon", "lat")
  proj4string(nonsites_temp) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
  # extracting predictor values for sites and nonsites
  sSP <- SpatialPoints(sites_temp@coords)
  nsSP <- SpatialPoints(nonsites_temp@coords)
  values_sites <- raster::extract(predictorstack, sSP)
  values_nonsites <- raster::extract(predictorstack, nsSP)
  # converting back to data.frame for modeling
  coords_sites <- sites_temp@coords
  coords_sites <- as.data.frame(coords_sites)
  coords_nonsites <- nonsites_temp@coords
  coords_nonsites <- as.data.frame(coords_nonsites)
  values_sites <- as.data.frame(values_sites)
  values_nonsites <- as.data.frame(values_nonsites)
  values_sites$site <- 1
  values_nonsites$site <- 0
  values_sites$lon <- coords_sites$lon
  values_sites$lat <- coords_sites$lat
  values_nonsites$lon <- coords_nonsites$lon
  values_nonsites$lat <- coords_nonsites$lat
  evidence <- rbind(values_sites, values_nonsites)
  evidence <- na.omit(evidence)
  return(evidence)
}

# function to draw sets of equal size
finalizeEvidence <- function(evd){
  siedl_pts <- filter(evd, site == 1)
  nons_pts <- filter(evd, site == 0)
  sz <- nrow(siedl_pts)
  nons_pts_sub <- sample_n(nons_pts, size = sz)
  temp <- rbind(siedl_pts, nons_pts_sub)
  return(temp)
}

evidence <- generateEvidence(sitesdata = sites, nonsitesdata = funtest, predictorstack = predictors)
evidence <- finalizeEvidence(evidence)

saveRDS(evidence, file = "Daten/evidence.csv")
evidence <- readRDS("Daten/evidence.csv")
# fitting gam and trying out kriging
names(evidence)

require(mgcv)
glmfit <- glm(site ~ dem + temp + rain + distance_water + frostdays + sunhours + 
                tpi + slope + as.factor(aspect), 
              family = binomial(), 
              data = evidence)

# predictive mapping 
df <- as.data.frame(predictors)
df[c("x","y")] <- coordinates(predictors)
pdata <- predict(glmfit, newdata = df, type = "response")
df$pdata <- pdata

x_pred <- predictors
x_pred$pred <- pdata

# predictive plot for glm
plot(x_pred$pred)

evidence$aspect <- as.factor(evidence$aspect)

fit <- gam(site ~ dem + temp + rain + distance_water + frostdays + sunhours + 
                 tpi + slope + aspect, 
               family = binomial, 
               data = evidence)

# how to include factor variable in predictors?
# df_gam <- df[complete.cases(df),]
preds <- predictors
# preds$aspect <- as.factor(preds$aspect)
pdatagam <- predict(predictors, fit, type = "response")
plot(pdatagam)

## gp spline default range
gp <- gam(site ~ s(lon, lat , bs="gp") + dem + temp + rain + 
                distance_water + frostdays + sunhours + tpi + slope + as.factor(aspect), 
              family = binomial, 
              data = evidence)  

fit <- brm(site ~ s(lon, lat) + dem + temp + rain + 
             distance_water + frostdays + sunhours + tpi + slope, 
           family = binomial, 
           data = evidence, chains = 4, cores = 4)

# modellwahl => aic
# k für gp smooth => was ist einfluss von k?

vis.gam(gp, view = c("lon", "lat"))

gp2 <- gam(site ~ s(lon, lat , bs="gp", k=50) + dem + temp + rain + 
               s(distance_water) + frostdays + sunhours + tpi + slope, 
              family = binomial, 
              data = evidence)  
vis.gam(gp2, view = c("lon", "lat"))
plot(gp2)
draw(gp2)
termplot(gp2)

# am 7. januar bei frau höfer abgeben
preds <- predictors
preds$lon <- coordinates(predictors)[,1]
preds$lat <- coordinates(predictors)[,2]
test <- predict(preds, gp2, type = "response")
plot(test)

# Low rank Gaussian process smooths
# 
# Gaussian process/kriging models based on simple covariance functions can be 
# written in a very similar form to thin plate and Duchon spline models (e.g. 
# Handcock, Meier, Nychka, 1994), and low rank versions produced by the eigen 
# approximation method of Wood (2003). Kammann and Wand (2003) suggest a 
# particularly simple form of the Matern covariance function with only a single 
# smoothing parameter to estimate, and this class implements this and other 
# similar models.
# 
# Usually invoked by an s(...,bs="gp") term in a gam formula. Argument m selects
# the covariance function, sets the range parameter and any power parameter. If
# m is not supplied then it defaults to NA and the covariance function 
# suggested by Kammann and Wand (2003) along with their suggested range 
# parameter is used. Otherwise m[1] between 1 and 5 selects the correlation 
# function from respectively, spherical, power exponential, and Matern with 
# kappa = 1.5, 2.5 or 3.5. m[2] if present specifies the range parameter, with
# non-positive or absent indicating that the Kammann and Wand estimate should
# be used. m[3] can be used to specify the power for the power exponential 
# which otherwise defaults to 1.
# 
# 1 spherical
# 2 power exponential
# 3 matern kappa 1.5
# 4 matern kappa 2.5
# 5 matern kappa 3.5
# k = 0.5 entspricht exponential modell
# https://stats.stackexchange.com/questions/322523/what-is-the-rationale-of-the-mat%C3%A9rn-covariance-function

