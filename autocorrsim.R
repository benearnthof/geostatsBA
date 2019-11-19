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