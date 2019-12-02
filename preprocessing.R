# This script serves for preprocessing the data.
# All results are available in the folder "Daten" only run this script once.
# Analysis based on a Dataset i got from:
# Dr. Stephan Luecke & Dr. Caroline von Nicolai. The dataset can be thought of as
# "Open Data" since it is part of the dissertation of Peer Fender (2017)
# Link to the dissertation:
# https://archiv.ub.uni-marburg.de/ubfind/Record/urn:nbn:de:hebis:04-z2017-0774/Holdings#tabnav
# Link to the Dataset:
# https://pma.gwi.uni-muenchen.de:8888/sql.php?server=9&db=vfpa_eisenzeit&table=fender_2017
# Username: EZD_statistik; Password: statistik@EZD
# I exported the data from the database as a .csv and did minimal preprocessing which i did not save.
# All i saved is the .RDS file.

# Importing the Data and cleaning everything up
raw_data <- readRDS("Daten/Daten_mit_Epoche.RDS")
raw_data <- raw_data[,!colnames(raw_data) %in% c("lng", "lat", "wkb_geom", "wkb_geom_wgs84")]
raw_data$lng_wgs84 <- as.numeric(as.vector(raw_data$lng_wgs84))
raw_data$lat_wgs84 <- as.numeric(as.vector(raw_data$lat_wgs84))
raw_data <- raw_data[order(raw_data$lng_wgs84, raw_data$lat_wgs84),]

# Loading data from Downloaded raster files and through R's interface for loading
# Raster and vector data

ger_file <- raster::getData("GADM", country = "Germany", level = 1)
plot(ger_file)
bay_file <- ger_file[match(toupper("Bayern"),toupper(ger_file$NAME_1)),]
# border of bavaria => useful for masking
plot(bay_file)

# Generating random raster as a basis for masking
bay_raster <- raster(xmn = 8.975, xmx = 13.84167, ymn = 47.26667, ymx = 50.56667, nrow = 396, ncol = 584)
bay_raster[] <- runif(396*584)
plot(bay_raster)
# masking
bay_mask <- mask(bay_raster, bay_file)
plot(bay_mask)

# Importing a DEM for germany through 'raster::getData'
dem <- raster::getData(name = "alt", country = "Germany")
names(dem) <- "dem"
# masking
dem <- crop(dem, bay_mask)
dem <- mask(dem, bay_mask)
plot(dem)
writeRaster(dem, "Daten/dem", overwrite = TRUE)

# Percipitation and average temperature
# these files are 72.3 MB large in total so it can take a while until the download is finished
weatherdata <- raster::getData("worldclim", var = "bio", res = 0.5, lon =  12.101624, lat = 49.013432)
# at the moment I just added average temperature and rain but variables like solar
# radiation might also be of interest later on
weatherdata <- weatherdata[[c(1, 12)]]
names(weatherdata) <- c("temp", "rain")

# masking and converting the temperature units from .1 of a degree to whole degrees celsius
weatherdata <- crop(weatherdata, bay_mask)
weatherdata <- mask(weatherdata, bay_mask)
weatherdata <- stack(weatherdata[[1]] / 10, weatherdata[[2]])
plot(weatherdata)

writeRaster(weatherdata[[1]], "Daten/temp_raster", overwrite = TRUE)
writeRaster(weatherdata[[2]], "Daten/rain_raster", overwrite = TRUE)

# calculating terrain => slope aspect and tpi
# dem <- predictors$dem
env_data <- terrain(dem, opt = c("slope", "aspect", "tpi"), unit = "degrees")
env_data$aspect <- ceiling((env_data$aspect + 360/8/2)/(360/8))
env_data$aspect[env_data$aspect > 8] <- 1
plot(env_data)
# Turn into Cathegorial Variable
# since for aspect the values for 360 and 0 degrees should be the same i convert to 8
# categories => North North-East East South-East South...

env_data <- crop(env_data, bay_mask)
env_data <- mask(env_data, bay_mask)
plot(env_data)
writeRaster(env_data[[1]], "Daten/tpi_raster", overwrite = TRUE)
writeRaster(env_data[[2]], "Daten/slope_raster", overwrite = TRUE)
writeRaster(env_data[[3]], "Daten/aspect_raster", overwrite = TRUE)

# calculating distance to the nearest water.
# source for shapefiles https://biogeo.ucdavis.edu/data/diva/wat/DEU_wat.zip
river_shape <- shapefile("Daten/DEU_water_lines_dcw.shx")
lake_shape <- shapefile("Daten/DEU_water_areas_dcw.shx")
# masking here takes ~ 1minute
river_raster <- mask(bay_raster, river_shape)
lake_raster <- mask(bay_raster, lake_shape)
# todo: add visualizations
# calculating distances
# since these functions calculate distance values for every single raster pixel
# they take a very long time (~1 hour on my computer)
# only run these lines once, the result is available in Daten/distance_water.RData
# distance_river <- distance(river_raster)
# distance_lake <- distance(lake_raster)
# save(distance_river, distance_lake, file = "Daten/distance_water.RData")
load("Daten/distance_water.RData")

# i am interested in the minimum distance to bodies of water
distance_water <- min(distance_lake, distance_river)
names(distance_water) <- "distance_water"
distance_water <- mask(distance_water, bay_mask)
plot(distance_water)
writeRaster(distance_water, "Daten/water_raster", overwrite = TRUE)
# should be easier to peer review

# frostdays per year
# (Source: https://maps.dwd.de/geoserver/web/wicket/bookmarkable/org.geoserver.web.demo.MapPreviewPage?0)
# https://maps.dwd.de/geoserver/web/wicket/bookmarkable/org.geoserver.web.demo.MapPreviewPage?0
# download as .tif and place in Daten directory

frostdays <- raster("Daten/dwd-Frostdays_annual_map_normals_1971_30.tif")
names(frostdays) <- "frostdays"
# reprojecting
frostdays <- projectRaster(frostdays, crs = crs(bay_mask))
frostdays <- raster::resample(frostdays, bay_mask)
frostdays <- mask(frostdays, bay_mask)
writeRaster(frostdays, "Daten/frostdays_raster", overwrite = TRUE)

# average sunhours per day per year
# (Source: https://maps.dwd.de/geoserver/web/wicket/bookmarkable/org.geoserver.web.demo.MapPreviewPage?0)
# https://maps.dwd.de/geoserver/web/wicket/bookmarkable/org.geoserver.web.demo.MapPreviewPage?0

sunhours <- raster("Daten/dwd-SDMS_17_1971_30.tif")
names(sunhours) <- "sunhours"
# reprojecting
sunhours <- projectRaster(sunhours, crs = crs(bay_mask))
sunhours <- raster::resample(sunhours, bay_mask)
sunhours <- mask(sunhours, bay_mask)

writeRaster(sunhours, "Daten/sunhours_raster", overwrite = TRUE)

# stacking all the rasters in a single predictorstack

predictors <- stack(c(
  "Daten/dem.grd",
  "Daten/temp_raster.grd",
  "Daten/rain_raster.grd",
  "Daten/water_raster.grd",
  "Daten/frostdays_raster.grd",
  "Daten/sunhours_raster.grd",
  "Daten/tpi_raster.grd",
  "Daten/slope_raster.grd",
  "Daten/aspect_raster.grd"))
