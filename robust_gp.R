# robust gaussian processes in STAN
library(rstan)
# utilityfile from 
# https://github.com/betanalpha/knitr_case_studies/blob/master/gaussian_processes/gp_part1/gp_utility.R
c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

c_light_trans <- c("#DCBCBC80")
c_light_highlight_trans <- c("#C7999980")
c_mid_trans <- c("#B97C7C80")
c_mid_highlight_trans <- c("#A2505080")
c_dark_trans <- c("#8F272780")
c_dark_highlight_trans <- c("#7C000080")

c_light_teal <- c("#6B8E8E")
c_mid_teal <- c("#487575")
c_dark_teal <- c("#1D4F4F")

# Plot Gaussian process realizations
plot_gp_realizations <- function(fit, data, true, title) {
  params <- extract(fit)
  I <- length(params$f_predict[,1])
  
  c_superfine <- c("#8F272705")
  plot(1, type="n", xlab="x", ylab="y", main=title,
       xlim=c(-10, 10), ylim=c(-10, 10))
  for (i in 1:I)
    lines(data$x_predict, params$f_predict[i,], col=c_superfine)
  
  points(data$x_predict, data$y_predict, col="white", pch=16, cex=0.6)
  points(data$x_predict, data$y_predict, col=c_mid_teal, pch=16, cex=0.4)
  lines(true$x_total, true$f_total, lwd=2, xlab="x", ylab="y")
  points(data$x, data$y, col="white", pch=16, cex=1.2)
  points(data$x, data$y, col="black", pch=16, cex=0.8)
}

# Plot Gaussian process predictive realizations
plot_gp_pred_realizations <- function(fit, data, true, title) {
  params <- extract(fit)
  I <- length(params$y_predict[,1])
  
  plot(1, type="n", xlab="x", ylab="y", main=title,
       xlim=c(-10, 10), ylim=c(-10, 10))
  for (i in 1:I)
    lines(data$x_predict, params$y_predict[i,], col=c_superfine)
  
  points(data$x_predict, data$y_predict, col="white", pch=16, cex=0.6)
  points(data$x_predict, data$y_predict, col=c_mid_teal, pch=16, cex=0.4)
  lines(true$x_total, true$f_total, lwd=2, xlab="x", ylab="y")
  points(data$x, data$y, col="white", pch=16, cex=1.2)
  points(data$x, data$y, col="black", pch=16, cex=0.8)
}

# Plot Gaussian process quantiles
plot_gp_quantiles <- function(fit, data, true, title) {
  params <- extract(fit)
  
  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:length(data$x_predict),
                 function(n) quantile(params$f_predict[,n], probs=probs))
  
  plot(1, type="n", main=title,
       xlab="x", ylab="y", xlim=c(-10, 10), ylim=c(-10, 10))
  polygon(c(data$x_predict, rev(data$x_predict)), c(cred[1,], rev(cred[9,])),
          col = c_light, border = NA)
  polygon(c(data$x_predict, rev(data$x_predict)), c(cred[2,], rev(cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(data$x_predict, rev(data$x_predict)), c(cred[3,], rev(cred[7,])),
          col = c_mid, border = NA)
  polygon(c(data$x_predict, rev(data$x_predict)), c(cred[4,], rev(cred[6,])),
          col = c_mid_highlight, border = NA)
  lines(data$x_predict, cred[5,], col=c_dark, lwd=2)
  
  points(data$x_predict, data$y_predict, col="white", pch=16, cex=0.6)
  points(data$x_predict, data$y_predict, col=c_mid_teal, pch=16, cex=0.4)
  lines(true$x_total, true$f_total, lwd=2, xlab="x", ylab="y", col="black")
  points(data$x, data$y, col="white", pch=16, cex=1.2)
  points(data$x, data$y, col="black", pch=16, cex=0.8)
}

# Plot Gaussian process predictive quantiles
plot_gp_pred_quantiles <- function(fit, data, true, title) {
  params <- extract(fit)
  
  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:length(data$x_predict),
                 function(n) quantile(params$y_predict[,n], probs=probs))
  
  plot(1, type="n", main=title,
       xlab="x", ylab="y", xlim=c(-10, 10), ylim=c(-10, 10))
  polygon(c(data$x_predict, rev(data$x_predict)), c(cred[1,], rev(cred[9,])),
          col = c_light, border = NA)
  polygon(c(data$x_predict, rev(data$x_predict)), c(cred[2,], rev(cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(data$x_predict, rev(data$x_predict)), c(cred[3,], rev(cred[7,])),
          col = c_mid, border = NA)
  polygon(c(data$x_predict, rev(data$x_predict)), c(cred[4,], rev(cred[6,])),
          col = c_mid_highlight, border = NA)
  lines(data$x_predict, cred[5,], col=c_dark, lwd=2)
  
  points(data$x_predict, data$y_predict, col="white", pch=16, cex=0.6)
  points(data$x_predict, data$y_predict, col=c_mid_teal, pch=16, cex=0.4)
  lines(true$x_total, true$f_total, lwd=2, xlab="x", ylab="y", col="black")
  points(data$x, data$y, col="white", pch=16, cex=1.2)
  points(data$x, data$y, col="black", pch=16, cex=0.8)
}

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

alpha_true <- 3
rho_true <- 5.5
sigma_true <- 2

N_total = 501
x_total <- 20 * (0:(N_total - 1)) / (N_total - 1) - 10

simu_data <- list(alpha=alpha_true, rho=rho_true, sigma=sigma_true,
                  N=N_total, x=x_total)

# sampling from the gaussian process
simu_fit <- stan(file='simu_gauss.stan', data=simu_data, iter=1,
                 chains=1, seed=494838, algorithm="Fixed_param")

# extracting from the simulation fit
f_total <- extract(simu_fit)$f[1,]
y_total <- extract(simu_fit)$y[1,]

true_realization <- data.frame(f_total, x_total)
names(true_realization) <- c("f_total", "x_total")

observed_idx <- c(50*(0:10)+1)
N = length(observed_idx)
x <- x_total[observed_idx]
y <- y_total[observed_idx]

plot(x_total, f_total, type="l", lwd=2, xlab="x", ylab="y",
     xlim=c(-10, 10), ylim=c(-10, 10))
points(x_total, y_total, col="white", pch=16, cex=0.6)
points(x_total, y_total, col=c_mid_teal, pch=16, cex=0.4)
points(x, y, col="white", pch=16, cex=1.2)
points(x, y, col="black", pch=16, cex=0.8)

N_predict <- N_total
x_predict <- x_total
y_predict <- y_total

# save the predicted and true data 
stan_rdump(c("N", "x", "y",
             "N_predict", "x_predict", "y_predict"),
           file="gp.data.R")

data <- read_rdump("gp.data.R")

stan_rdump(c("f_total", "x_total", "sigma_true"), file="gp.truth.R")

# construct the prior data generating process conditioned on the true realization
# of the gaussian process

f_data <- list(sigma=sigma_true, N=N_total, f=f_total)
dgp_fit <- stan(file='simu_gauss_dgp.stan', data=f_data, iter=1000, warmup=0,
                chains=1, seed=5838298, refresh=1000, algorithm="Fixed_param")

plot_gp_pred_quantiles(dgp_fit, data, true_realization,
                       "True Data Generating Process Quantiles")
