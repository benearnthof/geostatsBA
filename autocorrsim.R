# Accounting for spatial autocorrelation ====
# Step 1: Data generation
require(raster)
require(gstat)
require(lattice)
require(nlme)

set.seed(100)

# function to draw from a multivariate normal distribution
rmvn <- function(n, mu = 0, V = matrix(1)) {
  p <- length(mu)
  if (any(is.na(match(dim(V), p)))) {
    stop("Dimension problem!")
  }
  D <- chol(V)
  t(matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p)))
}

# setting up a square region
simgrid <- expand.grid(1:60, 1:60)
n <- nrow(simgrid)

# distance matrix
distance <- as.matrix(dist(simgrid))

# random variable to fill grid
delta <- 0.05
X <- rmvn(1, rep(0, n), exp(-delta * distance))

# converting to raster
Xraster <- rasterFromXYZ(cbind(simgrid[, 1:2] - 0.5, X))

# theoretic distance correlation:
plot(1:100, exp(-delta * 1:100), type = "l", xlab = "Dist", ylab = "Corr")
plot(Xraster)

elev <- raster::getData("alt", country = "CHE", mask = TRUE)
plot(elev)
extent(Xraster)
extent(elev)

# cropping raster
e <- as(extent(7.45, 7.95, 46.6, 47.1), "SpatialPolygons")
crs(e) <- crs(elev)
r <- crop(elev, e)
plot(r)
r@data@values <- r@data@values / 1000
plot(r)

shift <- shift(r, x = -7.45, y = -46.6)
plot(shift)
shift <- scale(shift)
plot(shift)
extent(shift) <- c(0, 60, 0, 60)
plot(shift)
length(shift@data@values)
elev <- shift
elevpx <- as(elev, "SpatialPixelsDataFrame")

# simulating the abundance distribution of a species ====
# quadratic effect of elevation
# Also a second set of expected values by including an effect of another
# spatially autocorrelated covariate

# defining the coefficients
beta0 <- 1
beta1 <- 3
beta2 <- -2

# abundance as a quadratic function of elevation
lambda1 <- exp(beta0 + beta1 * values(elev) + beta2 * values(elev)^2)
# abundance as a quadratic function of elevation + other spatially
# autocorrelated covariate
lambda2 <- exp(beta0 + beta1 * values(elev) + beta2 * values(elev)^2 + values(Xraster))

# Plot the results
par(mfrow = c(1, 2))
plot(values(elev), lambda1, cex = 0.5, main = expression(lambda == f(elevation)))
plot(values(elev), lambda2, cex = 0.5, main = expression(lambda == f(elevation, X)))

# points of the second response curve are much more dispersed around the quadratic
# effect of elevation (as expected)

# compute actual counts using the second set of expected values.
# => This is what we would observe if we would visit every pixel and count all
# individuals from that species with detection probability = 1

# Determine the actual counts
counts <- rpois(n, lambda2)
counts <- rasterFromXYZ(cbind(coordinates(elev), counts))

par(mfrow = c(1, 1))
plot(counts, main = "Abundance distribution")

# lets randomly choose 200 pixels and extract the related abundance values
# => presence data only available at these 200 pixels

coords <- coordinates(counts)
fulldata <- data.frame(coords,
  elevation = extract(elev, coords),
  Xvar = extract(Xraster, coords), counts = extract(counts, coords)
)
sites <- sample(1:n, 200)
sampledata <- fulldata[sites, ]
head(sampledata)
# visualizing the results
plot(counts, main = "Abundance distribution with sampling sites")
points(sampledata[, 1:2], pch = 16)

# Analysis using GLS ====
# Best way to model the residual spatial autocorrelation would be to use generalized
# least squares.
# The GLS directly models the spatial covariance structure in the variance covariance matrix
# using parametric functions.
# Let's use this model as a baseline

m2 <- gls(log1p(counts) ~ elevation + I(elevation^2), data = sampledata)
vario2 <- Variogram(m2, form = ~ x + y, resType = "pearson")
plot(vario2, smooth = TRUE, ylim = c(0.4, 1.3))

# semi-variance increases with distance => spatial autocorrelation is present in the residuals

# fitting the model with spatial correlation structure => gls function using correlation argument
# nugget argument is for the intercept

m3 <- gls(log1p(counts) ~ elevation + I(elevation^2), 
          correlation = corExp(form = ~ x + y, nugget = TRUE), data = sampledata)
m4 <- gls(log1p(counts) ~ elevation + I(elevation^2),
          correlation = corGaus(form = ~ x + y, nugget = TRUE), data = sampledata)
m5 <- gls(log1p(counts) ~ elevation + I(elevation^2), 
          correlation = corSpher(form = ~ x + y, nugget = TRUE), data = sampledata)
m6 <- gls(log1p(counts) ~ elevation + I(elevation^2), 
          correlation = corLin(form = ~ x + y, nugget = TRUE), data = sampledata)
m7 <- gls(log1p(counts) ~ elevation + I(elevation^2), 
          correlation = corRatio(form = ~x + y, nugget = TRUE), data = sampledata)

AIC(m2, m3, m4, m5, m6, m7)

# corRatio and corExp seem to be the best models
# should be expected since we generated missing covariate with the exponential function

vario3 <- Variogram(m3, from = ~x + y, resType = "pearson")
plot(vario3, smooth = FALSE)

# looks pretty solid 
# plotting a sample variogram of the normalized residuals
# These residuals are standardized residuals pre-multiplied by the inverse square-root
# factor fo the estimated error correlation matrix. 
# Thus, if the residual spatial autocorrelation was properly accounted for by m3, 
# we should not see any trend in this new variogram. 
vario4 <- Variogram(m3, form = ~x + y, resType = "normalized", maxDist = 60)
plot(vario4, smooth = FALSE)

# more or less horizontal band of points => the model seems to fit the data reasonably well
# lets start to make inferences based on these results

# Accounting for spatial autocorrelation using the gls function allows us to get 
# unbiased parameter and uncertainty estimates, the modelled correlation structure
# however is not directly use when making predictions. If we want to use the extra 
# information provided by the spatial autocorrelation to make better predictions, 
# we need more complex models. 

resp <- expm1(predict(m3, newdata = data.frame(elevation = values(elev), x = coordinates(elev)[,1], 
                                               y = coordinates(elev)[2])))
resp <- rasterFromXYZ(cbind(coordinates(elev), resp))
plot(resp, main = "Predicted abundances")

# looks promising so far

# Analysis using ordinary kriging ====
# (Gaussian process regression or best linear unbiased prediction) 
# Basic Idea:
# The value of interest at an unknown point is computed as a weighted average 
# of the sampled neighbours. The weights are defined by the variogram model. 
# Ordinary kriging assumes that the data comes from a multivariate normal 
# distribution => log transform the generated example data

datsp <- sampledata
coordinates(datsp) <- c("x", "y")

# Sample variogram of the data without covariate
vario <- variogram(log1p(counts) ~ 1, datsp)
plot(vario)

# Ordinary kriging only requires the relationship between similarity and distance
# no covariates are introduced as of yet. This relationship is computed from the 
# sample variogram. To get predictions we need to fit a theoretical variogram 
# model to our sample variogram. 
# Fit a non-linear regression to the variogram using fit.variogram function of 
# gstat package. The vgm function allows us to specify the theoretical model 
# => exponential, gaussian, spherical, etc. 
# and also specify initial values for the variogram parameters.

# Fit theoretical variogram model to our sample variogram
counts.ivgm <- vgm(model = "Exp", nugget = 0, range = 15, psill = 1)
counts.vgm <- fit.variogram(vario, model = counts.ivgm)
counts.vgm

plot(vario, counts.vgm, pch = "+", main = "log1p(counts)", cex = 2)

# lets make predictions using ordinary kriging
# ~1 indicates ordinary kriging, second argument is SpatialPointsDataFrame containing
# the locations of the sampled sites with the corresponding values of the response variable
# third argument is Spatial object => spatial pixels data frame containing the 
# coordinates of the sites for which we want to have predictions. 
# last argument is fitted variogram model. 

# krige function returns SpatialPixelsDataFrame with 2 columns 
# Predictions and prediction variances. 
# => need to back transform the predictions after the log transform

# Ordinary kriging 
# crs(datsp) <- crs(elevpx)
# does not work properly lets ignore crs for now
crs(datsp) <- NA
crs(elevpx) <- NA
counts.ok <- krige(log1p(counts) ~ 1, locations = datsp, newdata = elevpx, model = counts.vgm)

preds <- stack(counts.ok)
plot(preds, main = c("Ordinary kriging predictions", "Prediction variance"))

# well that looks like shit lmao
# back transform the predictions
plot(exp(preds[[1]]) - 1, main = "Back-transformed Predictions")

# at least we tried. 
# no covariate => predictions are strongly smoothed and most of the time not very accurate
# => This method does not work well if the distribution of a species is mostly determined 
# by covariates with sharp boundaries. 

# Analysis using universal kriging (regression kriging) ====
# (Universal Kriging; Kriging with external drift.)
# Combines information from the covariates and the residual spatial autocorrelation
# Also assumes the data is drawn from a multivariate normal distribution. 
# Model fitting is identical to the GLS with covariates, but prediction is different. 

# Fit a simple linear model
m <- lm(log1p(counts) ~ elevation + I(elevation^2), data = datsp)
# Compute the sample variogram of the residuals: residual spatial
# autocorrelation
vario.rs <- variogram(residuals(m) ~ 1, datsp)

# Fit a theoretical variogram model to the residuals, and give initial
# values
counts.rivgm <- vgm(model = "Exp", nugget = 0, range = 10, psill = 0.5)
counts.rvgm <- fit.variogram(vario.rs, model = counts.rivgm)
counts.rvgm

plot(vario.rs, counts.rvgm, pc = "+", main = "Residuals", cex = 2)

# we got the theoretical model describing the residual spatial autocorrelation
# use krige function with covariates in the formula
# When using universal kriging the third argument must not only contain the coordinates 
# of the predicted sites, but also the value of the covariates for these sites. 

# regression-kriging
names(elevpx@data) <- "elevation"
counts.uk <- krige(log1p(counts) ~ elevation + I(elevation^2), locations = datsp,
                   newdata = elevpx, model = counts.rvgm)
preds <- stack(counts.uk)
plot(preds, main = c("Universal kriging predictions", "Prediction variance"))

# looks a lot better already!
plot(exp(preds[[1]]) - 1, main = "Predictions (back-transformed")
# wow it seems we are able to account for most of the structure in the data!