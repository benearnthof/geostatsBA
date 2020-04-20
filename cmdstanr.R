# installing cmdstanr to use covariance functions
devtools::install_github("stan-dev/cmdstanr")
library(cmdstanr)
# requires working version of cmdstan
# install_cmdstan()
# verify that version can be found
cmdstan_version()

# executing example model
file <- file.path(cmdstan_path(), "examples", "bernoulli", "bernoulli.stan")
mod <- cmdstan_model(file)

# R object 
mod$print()
mod$exe_file()

data_list <- list(N = 10, y = c(0,1,0,0,0,0,0,0,0,1))
fit <- mod$sample(
  data = data_list, 
  seed = 123, 
  num_chains = 2, 
  num_cores = 2
)

fit$summary()
fit$cmdstan_diagnose()
fit$cmdstan_summary()

# creating stanfit object
stanfit <- rstan::read_stan_csv(fit$output_files())
print(stanfit)

library(mgcv)
library(brms)
set.seed(1)
dat <- gamSim(1, n = 100)
# power exponential 
m_exponential <- gam(y ~ s(x1, x2, bs = "gp", m = 2), data = dat)
plot(m_exponential)

m_matern_1_5 <- gam(y ~ s(x1, x2, bs = "gp", m = 3), data = dat)
plot(m_matern_1_5)

m_matern_3_5 <- gam(y ~ s(x1, x2, bs = "gp", m = 5), data = dat)
plot(m_matern_3_5)

b_exponential <- brm(y ~ gp(x1, x2), data = dat, 
                        chains = 2, cores = 2, iter = 400) # uses power exponential by default
# 150 seconds for 1000 transitions using 10 steps
# between 10 and 20 seconds on the server
summary(b_exponential)
plot(b_exponential)

# seems good now lets try customizing our stancode

stancode <- brms::make_stancode(y ~ gp(x1, x2), data = dat, 
                                chains = 2, cores = 2, iter = 400)
standata <- brms::make_standata(y ~ gp(x1, x2), data = dat, 
                                chains = 2, cores = 2, iter = 400)

# editing covariance function in stancode and then using cmdstanr
stancode

file <- paste0(getwd(), "/cmdstan_matern32.stan")
mod <- cmdstan_model(file)

data_mod <- list()

for (i in 1:length(standata)) {
  data_mod[[i]] <- standata[[i]]
}
names(data_mod) <- names(standata)

fit <- mod$sample(
  data = data_mod, 
  seed = 123, 
  num_chains = 2, 
  num_cores = 2,
  num_samples = 400,
  num_warmup = 200
)

fit$summary()
fit$cmdstan_summary()
test <- fit$draws()

# something is fucky

file <- paste0(getwd(), "/cmdstan_exponential.stan")
mod <- cmdstan_model(file)

fit_expo <- mod$sample(
  data = data_mod, 
  seed = 123, 
  num_chains = 2, 
  num_cores = 2,
  num_samples = 400,
  num_warmup = 200
)

b_exponential
# sdgp(gpx1x2)       3.70    0.68  
# lscale(gpx1x2)     0.13    0.02 
# Intercept          7.76    1.33
# sigma              2.19    0.19
fit_expo$summary()

# sdgp_1[1]         3.85    0.800 here median was closer
# lscale_1[1,1]     0.134   0.0248
# Intercept         7.49    1.34
# sigma             2.21    0.211

# estimates are reasonably close for the relatively small sample size

draws_expo <- fit_expo$draws()

expo_stanfit <- rstan::read_stan_csv(fit_expo$output_files())
print(expo_stanfit)
twochainz <- expo_stanfit@sim$samples
str(twochainz[[1]])
is.data.frame(twochainz[[1]])
is.data.frame(twochainz[[2]])

head2d = function(x, n = 6L, ...) {
  indvecs = lapply(seq_along(dim(x)), function(i) {
    if (length(n) >= i) {
      ni = n[i]
    } else {
      ni =  dim(x)[i]
    }
    if (ni < 0L)
      ni = max(nrow(x) + ni, 0L)
    else
      ni = min(ni, dim(x)[i])
    seq_len(ni)
  })
  lstargs = c(list(x),indvecs, drop = FALSE)
  do.call("[", lstargs)
}

head2d(twochainz[[1]], n = c(5,5))
chain1 <- twochainz[[1]]
chain2 <- twochainz[[2]]
df1 <- chain1[, -grep("zgp_", colnames(chain1))]
df1 <- df1[, -grep("mu.", colnames(df1))]
df1 <- df1[, 1:4]

getrels <- function(chain) {
  ret <- chain[, -grep("zgp_", colnames(chain))]
  ret <- ret[, -grep("mu.", colnames(ret))]
  ret <- ret[, 1:4]
  return(ret)
}

c1 <- getrels(twochainz[[1]])
c2 <- getrels(twochainz[[2]])

chainz <- rbind(c1, c2)
colnames(chainz) <- c("Intercept", "Sigma_GP", "Lscale", "Sigma")
mlt <- reshape2::melt(chainz)

library(ggplot2)
ggplot(mlt, aes(x = value)) +
  geom_density() +
  facet_wrap(~variable, scales = "free", ncol = 1)

plot(b_exponential)

# lets compare the matern covariance functions
# lets do more data and more samples 

set.seed(1)
dat <- gamSim(1, n = 300)

m_matern_1_5 <- gam(y ~ s(x1, x2, bs = "gp", m = 3), data = dat)
plot(m_matern_1_5) # Intercept estimate 7.7020 std error 0.1263


standata <- brms::make_standata(y ~ gp(x1, x2), data = dat, 
                                chains = 2, cores = 2, iter = 400)

file <- paste0(getwd(), "/cmdstan_matern32.stan")
mod <- cmdstan_model(file)

data_mod <- list()

for (i in 1:length(standata)) {
  data_mod[[i]] <- standata[[i]]
}
names(data_mod) <- names(standata)

fit <- mod$sample(
  data = data_mod, 
  seed = 123, 
  num_chains = 2, 
  num_cores = 2,
  num_samples = 400,
  num_warmup = 200
)

set.seed(1)
evd <- evidence[sample(nrow(evidence), size = 500),]
standata <- brms::make_standata(site ~ gp(lon, lat), data = evd,
                                family = bernoulli)

file <- paste0(getwd(), "/cmdstan_matern32.stan")
mod <- cmdstan_model(file)

data_mod <- list()

for (i in 1:length(standata)) {
  data_mod[[i]] <- standata[[i]]
}
names(data_mod) <- names(standata)

fit <- mod$sample(
  data = data_mod, 
  seed = 123, 
  num_chains = 4, 
  num_cores = 4,
  num_samples = 400,
  num_warmup = 200
)
