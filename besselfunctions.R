# bessel functions
nus <- c(0:5, 10, 20)

x <- seq(0, 4, length.out = 501)
plot(x, x, ylim = c(0, 6), ylab = "", type = "n",
     main = "Bessel Functions  I_nu(x)")
for(nu in nus) lines(x, besselI(x, nu = nu), col = nu + 2)
legend(0, 6, legend = paste("nu=", nus), col = nus + 2, lwd = 1)

x <- seq(0, 40, length.out = 801); yl <- c(-.5, 1)
plot(x, x, ylim = yl, ylab = "", type = "n",
     main = "Bessel Functions  J_nu(x)")
abline(h=0, v=0, lty=3)
for(nu in nus) lines(x, besselJ(x, nu = nu), col = nu + 2)
legend("topright", legend = paste("nu=", nus), col = nus + 2, lwd = 1, bty="n")

## Negative nu's --------------------------------------------------
xx <- 2:7
nu <- seq(-10, 9, length.out = 2001)
## --- I() --- --- --- ---
matplot(nu, t(outer(xx, nu, besselI)), type = "l", ylim = c(-50, 200),
        main = expression(paste("Bessel ", I[nu](x), " for fixed ", x,
                                ",  as ", f(nu))),
        xlab = expression(nu))
abline(v = 0, col = "light gray", lty = 3)
legend(5, 200, legend = paste("x=", xx), col=seq(xx), lty=1:5)

## --- J() --- --- --- ---
bJ <- t(outer(xx, nu, besselJ))
matplot(nu, bJ, type = "l", ylim = c(-500, 200),
        xlab = quote(nu), ylab = quote(J[nu](x)),
        main = expression(paste("Bessel ", J[nu](x), " for fixed ", x)))
abline(v = 0, col = "light gray", lty = 3)
legend("topright", legend = paste("x=", xx), col=seq(xx), lty=1:5)

## ZOOM into right part:
matplot(nu[nu > -2], bJ[nu > -2,], type = "l",
        xlab = quote(nu), ylab = quote(J[nu](x)),
        main = expression(paste("Bessel ", J[nu](x), " for fixed ", x)))
abline(h=0, v = 0, col = "gray60", lty = 3)
legend("topright", legend = paste("x=", xx), col=seq(xx), lty=1:5)


##---------------  x --> 0  -----------------------------
x0 <- 2^seq(-16, 5, length.out=256)
plot(range(x0), c(1e-40, 1), log = "xy", xlab = "x", ylab = "", type = "n",
     main = "Bessel Functions  J_nu(x)  near 0\n log - log  scale") ; axis(2, at=1)
for(nu in sort(c(nus, nus+0.5)))
  lines(x0, besselJ(x0, nu = nu), col = nu + 2, lty= 1+ (nu%%1 > 0))
legend("right", legend = paste("nu=", paste(nus, nus+0.5, sep=", ")),
       col = nus + 2, lwd = 1, bty="n")

x0 <- 2^seq(-10, 8, length.out=256)
plot(range(x0), 10^c(-100, 80), log = "xy", xlab = "x", ylab = "", type = "n",
     main = "Bessel Functions  K_nu(x)  near 0\n log - log  scale") ; axis(2, at=1)
for(nu in sort(c(nus, nus+0.5)))
  lines(x0, besselK(x0, nu = nu), col = nu + 2, lty= 1+ (nu%%1 > 0))
legend("topright", legend = paste("nu=", paste(nus, nus + 0.5, sep = ", ")),
       col = nus + 2, lwd = 1, bty="n")

x <- x[x > 0]
plot(x, x, ylim = c(1e-18, 1e11), log = "y", ylab = "", type = "n",
     main = "Bessel Functions  K_nu(x)"); axis(2, at=1)
for(nu in nus) lines(x, besselK(x, nu = nu), col = nu + 2)
legend(0, 1e-5, legend=paste("nu=", nus), col = nus + 2, lwd = 1)

yl <- c(-1.6, .6)
plot(x, x, ylim = yl, ylab = "", type = "n",
     main = "Bessel Functions  Y_nu(x)")
for(nu in nus){
  xx <- x[x > .6*nu]
  lines(xx, besselY(xx, nu=nu), col = nu+2)
}
legend(25, -.5, legend = paste("nu=", nus), col = nus+2, lwd = 1)

## negative nu in bessel_Y -- was bogus for a long time
curve(besselY(x, -0.1), 0, 10, ylim = c(-3,1), ylab = "")
for(nu in c(seq(-0.2, -2, by = -0.1)))
  curve(besselY(x, nu), add = TRUE)
title(expression(besselY(x, nu) * "   " *
                   {nu == list(-0.1, -0.2, ..., -2)}))

# animation done in python

## Matern spline default range
gp <- gam(site ~ s(lon, lat , bs="gp", k=50) + dem + temp + rain + 
            distance_water + frostdays + sunhours + tpi + slope + as.factor(aspect), 
          family = binomial, 
          data = evidence)  

vis.gam(gp, view = c("lon", "lat"))

gp2 <- gam(site ~ s(lon, lat , bs="gp", k=50) + dem + temp + rain + 
             s(distance_water) + frostdays + sunhours + tpi + slope, 
           family = binomial, 
           data = evidence)  
vis.gam(gp2, view = c("lon", "lat"))
plot(gp2)
draw(gp2)
termplot(gp2)

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

