# generalized additive models code
require(gamair)
data(engine)
attach(engine)
plot(size, wear, xlab = "Engine capacity", ylab = "Wear index")

# write a function defining tent function basis
tf <- function(x, xj, j) {
  ## generate jth tent function from set defined by knots xj
  dj <- xj * 0
  dj[j] <- 1
  approx(xj, dj, x)$y
}

# and use it to write an R function that will take a sequence of knots and
# an array of x values to produce a model matrix for the piecewise linear
# function.
tf.X <- function(x, xj) {
  ## tent function basis matrix given data x
  ## and knot sequence xj
  nk <- length(xj)
  n <- length(x)
  X <- matrix(NA, n, nk)
  for (j in 1:nk) X[, j] <- tf(x, xj, j)
  X
}

# fit a rank k = 6 basis with evenly spread knots over range of data
sj <- seq(min(size), max(size), length = 6) # generate knots
X <- tf.X(size, sj) # generate model matrix
X
# fit a lm in every basis interval
b <- lm(wear ~ X - 1)
# generate prediction data
s <- seq(min(size), max(size), length = 200)
# get prediction matrix
xp <- tf.X(s, sj)
plot(size, wear)
# plot predicted
lines(s, xp %*% coef(b))

# write a function for fitting a penalized piecewise smoother
prs.fit <- function(y, x, xj, sp) {
  X <- tf.X(x, xj) # get model matrix
  D <- diff(diag(length(xj)), differences = 2) # sqrt penalty
  X <- rbind(X, sqrt(sp) * D) # augmented model matrix
  Y <- c(y, rep(0, nrow(D))) # augmented data
  lm(Y ~ X - 1)
}

# k = 20, evenly spread knots
sj <- seq(min(size), max(size), length = 20) ## knots
b <- prs.fit(wear, size, sj, 10) ## penalized fit
plot(size, wear) ## plot data
Xp <- tf.X(s, sj) ## prediction matrix
lines(s, Xp %*% coef(b)) ## plot the smooth

# simple direct search for the GCV optimal smoothing parameter
rho <- seq(-9, 11, length = 90)
n <- length(wear)
V <- rep(NA, 90)

# loop through all smoothing parameters
for (i in seq_len(90)) {
  b <- prs.fit(wear, size, sj, exp(rho[i]))
  trF <- sum(influence(b)$hat[1:n]) # extract EDF
  rss <- sum((wear - fitted(b)[1:n])^2) # residual SS
  V[i] <- n * rss / (n - trF)^2
}

plot(V)
which.min(V)
lambda <- exp(rho[54])
plot(size, wear)
b <- prs.fit(wear, size, sj, lambda)
lines(s, Xp %*% coef(b))

## Additive models ====

# produce constrained model matrices for X and D
tf.XD <- function(x, xk, cmx = NULL, m = 2) {
  ## get X and D subject to constraint
  nk <- length(xk)
  X <- tf.X(x, xk)[, -nk] ## basis matrix
  D <- diff(diag(nk), differences = m)[, -nk] ## root penalty
  if (is.null(cmx)) cmx <- colMeans(X)
  X <- sweep(X, 2, cmx) ## subtract cmx from columns
  list(X = X, D = D, cmx = cmx)
}

attach(trees)
require(mgcv)
model <- gam(Volume ~ s(Girth) + s(Height))
predicted <- predict(model, newdata = list(Girth = Girth, Height = Height))

plot(predicted ~ Volume)
summary(model)

model2 <- brm(Volume ~ t2(Girth, Height), data = trees)
plot(model2)
require(brms)
ms <- marginal_smooths(model2)
plot(ms, stype = "raster")

model3 <- gam(Volume ~ s(Height) + s(Girth), family = Gamma(link = log))
plot(model3)
gam.check(model3)

plot(density(Volume))
hist(Volume)
hist(log(Volume))

library(mgcViz)
v <- getViz(model3)
plot(v)

b <- gam(Volume ~ t2(Height, Girth))
b <- getViz(b)
plot(b) + l_fitRaster() + l_fitContour() + l_points()
ck1 <- check1D(b, "Girth")
ck1 + l_densCheck()

ck1 <- check2D(b, x1 = "Height", x2 = "Girth")
ck1 + l_gridCheck2D(gridFun = mean)

# cubic regression splines
ct2 <- gam(Volume ~ s(Height,bs="cr") + s(Girth,bs="cr"), 
           family = Gamma(link = log), data = trees)
ct2
plot(ct2)


bspline <- function(x,k,i,m=2) {
  # evaluate ith B-spline basis function of order m at the
  # values in x, given knot locations in k
 if (m==-1) { # base of recursion
  res <- as.numeric(x<k[i+1]&x>=k[i])
  } else { # construct from call to lower order basis
  z0 <- (x-k[i])/(k[i+m+1]-k[i])
  z1 <- (k[i+m+2]-x)/(k[i+m+2]-k[i+1])
  res <- z0*bspline(x,k,i,m-1)+ z1*bspline(x,k,i+1,m-1)
}
  res
}
