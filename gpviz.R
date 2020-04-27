# Demo of Gaussian process regression with R
# James Keirstead
# 5 April 2012

# Chapter 2 of Rasmussen and Williams's book `Gaussian Processes
# for Machine Learning' provides a detailed explanation of the
# math for Gaussian process regression.  It doesn't provide
# much in the way of code though.  This Gist is a brief demo
# of the basic elements of Gaussian process regression, as
# described on pages 13 to 16.


# Load in the required libraries for data manipulation
# and multivariate normal distribution
require(MASS)
require(plyr)
require(reshape2)
require(ggplot2)

# Set a seed for repeatable plots
set.seed(12345)

# Calculates the covariance matrix sigma using a
# simplified version of the squared exponential function.
#
# Although the nested loops are ugly, I've checked and it's about
# 30% faster than a solution using expand.grid() and apply()
#
# Parameters:
# 	X1, X2 = vectors
# 	l = the scale length parameter
# Returns:
# 	a covariance matrix
calcSigma <- function(X1, X2, l = 1, type = c("exp", "matern")) {
  Sigma <- matrix(rep(0, length(X1) * length(X2)), nrow = length(X1))

  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      if (type == "exp") {
        # exponential
        Sigma[i, j] <- exp(-0.5 * (abs(X1[i] - X2[j]) / l)^2)
      } else if (type == "matern1") {
        # matern
        Sigma[i, j] <- exp(-1 * (abs(X1[i] - X2[j]) / l))
      } else if (type == "matern3") {
        Sigma[i, j] <- (1 + ((sqrt(3) * abs(X1[i] - X2[j]))/l)) * (exp(-1 * ((sqrt(3) * abs(X1[i] - X2[j]))/l)))
      } else if (type == "matern5") {
        Sigma[i, j] <- (1 + ((sqrt(5) * abs(X1[i] - X2[j]))/l) + (5 * (abs(X1[i] - X2[j]))^2)/3*l^2) * (exp(-1 * ((sqrt(5) * abs(X1[i] - X2[j]))/l)))
      }
      
    }
  }
  return(Sigma)
}

# 1. Plot some sample functions from the Gaussian process
# as shown in Figure 2.2(a)

# Define the points at which we want to define the functions
x.star <- seq(-5, 5, len = 300)

# Calculate the covariance matrix
sigma <- calcSigma(x.star, x.star, type = "matern5")

# Generate a number of functions from the process
n.samples <- 3
set.seed(1)
values <- matrix(rep(0, length(x.star) * n.samples), ncol = n.samples)
for (i in 1:n.samples) {
  # Each column represents a sample from a multivariate normal distribution
  # with zero mean and covariance sigma
  values[, i] <- mvrnorm(1, rep(0, length(x.star)), sigma)
}
values <- cbind(x = x.star, as.data.frame(values))
values <- melt(values, id = "x")

# Plot the result
fig2a <- ggplot(values, aes(x = x, y = value)) +
  #geom_rect(xmin = -Inf, xmax = Inf, ymin = -3, ymax = 3, fill = "grey90") +
  geom_line(aes(group = variable, color = variable), size = 0.7) +
  theme_bw() +
  scale_y_continuous(lim = c(-3, 3), name = "f(x)") +
  xlab("x") +
  theme(legend.position = "none")

fig2a



# 2. Now let's assume that we have some known data points;
# this is the case of Figure 2.2(b). In the book, the notation 'f'
# is used for f$y below.  I've done this to make the ggplot code
# easier later on.
f <- data.frame(
  x = c(-4, -3, -1, 0, 2, 4.5),
  y = c(-2, 0, 1, 2, -1, -2.5)
)

# Calculate the covariance matrices
# using the same x.star values as above
x <- f$x
k.xx <- calcSigma(x, x, type = "exp")
k.xxs <- calcSigma(x, x.star, type = "exp")
k.xsx <- calcSigma(x.star, x, type = "exp")
k.xsxs <- calcSigma(x.star, x.star, type = "exp")

# These matrix calculations correspond to equation (2.19)
# in the book.
f.star.bar <- k.xsx %*% solve(k.xx) %*% f$y
cov.f.star <- k.xsxs - k.xsx %*% solve(k.xx) %*% k.xxs

# This time we'll plot more samples.  We could of course
# simply plot a +/- 2 standard deviation confidence interval
# as in the book but I want to show the samples explicitly here.
n.samples <- 50
values <- matrix(rep(0, length(x.star) * n.samples), ncol = n.samples)
for (i in 1:n.samples) {
  values[, i] <- mvrnorm(1, f.star.bar, cov.f.star)
}
values <- cbind(x = x.star, as.data.frame(values))
values <- melt(values, id = "x")

# Plot the results including the mean function
# and constraining data points
hardxd <- data.frame(x.star, f.star.bar)
fig2b <- ggplot(values, aes(x = x, y = value)) +
  geom_line(aes(group = variable), colour = "black", alpha = 0.4, size = 0.7) +
  geom_line(data = hardxd, aes(x = x.star, y = f.star.bar), colour = "red") +
  geom_point(data = f, aes(x = x, y = y)) +
  theme_bw() +
  scale_y_continuous(lim = c(-4, 4), name = "f(x)") +
  scale_x_continuous(lim = c(-5, 5)) +
  xlab("x")

fig2b

# 3. Now assume that each of the observed data points have some
# normally-distributed noise.

# The standard deviation of the noise
sigma.n <- 0.1

# Recalculate the mean and covariance functions
f.bar.star <- k.xsx %*% solve(k.xx + sigma.n^2 * diag(1, ncol(k.xx))) %*% f$y
cov.f.star <- k.xsxs - k.xsx %*% solve(k.xx + sigma.n^2 * diag(1, ncol(k.xx))) %*% k.xxs

# Recalulate the sample functions
values <- matrix(rep(0, length(x.star) * n.samples), ncol = n.samples)
for (i in 1:n.samples) {
  values[, i] <- mvrnorm(1, f.bar.star, cov.f.star)
}
values <- cbind(x = x.star, as.data.frame(values))
values <- melt(values, id = "x")

xd <- cbind(x.star, f.bar.star)
xd <- as.data.frame.matrix(xd)
names(xd) <- c("x", "y")
# Plot the result, including error bars on the observed points
ggplot(values, aes(x = x, y = value)) +
  geom_line(aes(group = variable), colour = "black", alpha = 0.35) +
  geom_line(data = xd, aes(x = x, y = y), colour = "red") +
  geom_errorbar(data = f, aes(x = x, y = NULL, ymin = y - 2 * sigma.n, ymax = y + 2 * sigma.n), width = 0.2, color = "blue") +
  geom_point(data = f, aes(x = x, y = y)) +
  theme_bw() +
  scale_y_continuous(lim = c(-3, 3), name = "output, f(x)") +
  xlab("input, x")


# sampling from a 2d gaussian process
library("ggplot2")
library("plgp")

# kernel function
rbf_D <- function(X, l = 1, eps = sqrt(.Machine$double.eps)) {
  D <- plgp::distance(X)
  Sigma <- exp(-D / l)^2 + diag(eps, nrow(X))
}
# number of samples
nx <- 40
x <- seq(-5, 5, length = nx)
# grid of pairwise values
X <- expand.grid(x, x)
# compute squared exponential kernel on pairwise values
Sigma <- rbf_D(X, l = 4)

# sample from multivariate normal with mean zero, sigma = sigma
set.seed(5)
Y <- MASS::mvrnorm(1, rep(0, dim(Sigma)[1]), Sigma)

# plot results
pp <- data.frame(y = Y, x1 = X[, 1], x2 = X[, 2])
# ggplot(pp,aes(x=x1,y=x2)) +
#   geom_raster(aes(fill=y), interpolate = TRUE) +
#   geom_contour(aes(z=y), bins = 12, color = "gray30",
#                size = 0.5, alpha = 0.5) +
#   coord_equal() +
#   scale_x_continuous(expand=c(0,0)) +
#   scale_y_continuous(expand=c(0,0)) +
#   scale_fill_viridis_c(option = "viridis")

library(plotly)
library(plot3D)

require(akima)
require(rgl)
s <- interp(pp$x1, pp$x2, pp$y)
# surface3d(s$x, s$y, s$z)

# persp(s$x, s$y, s$z)


persp3D(s$x, s$y, s$z,
  scale = FALSE, theta = 60, ticktype = "detailed",
  facets = TRUE, colkey = FALSE, phi = 30
)
