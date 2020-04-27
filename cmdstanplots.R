# installing cmdstanr to use covariance functions
#devtools::install_github("stan-dev/cmdstanr")
#devtools::install_github("cran/StanHeaders")
#install.packages("rstan")
library(cmdstanr)
library(rstan)
library(brms)
#install.packages("brms")
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
set.seed(1)
dat <- gamSim(1, n = 500)
m_matern32<- gam(y ~ s(x1, x2, bs = "gp", m = 3), data = dat)
summary(m_matern32)

stancode <- brms::make_stancode(y ~ gp(x1, x2, k = 30, c = 5/4), data = dat)
standata <- brms::make_standata(y ~ gp(x1, x2, k = 30, c = 5/4), data = dat)

modeldata <- list()
for (i in 1:length(standata)) {
  modeldata[[i]] <- standata[[i]]
}
names(modeldata) <- names(standata)
setwd("/home/rstudio/geostatsBA/modeldump")
file <- paste0(getwd(), "/test.stan")
setwd("/home/rstudio/geostatsBA/modeldump")
mod <- cmdstan_model(file)

b_matern32 <- mod$sample(
  data = modeldata, 
  seed = 123, 
  num_chains = 4, 
  num_cores = 4,
  num_samples = 400,
  num_warmup = 400,
  refresh = 1
)

# 300 data points 4 chains
# 400 samples 400 warmup finished with mean chain execution time of 
# 9434.5 seconds = 2.6 hours
# # 30 of 1600 transitions ended with a divergence on default 0.8 adapt delta
saveRDS(b_matern32, file = "toydata_matern32.RDS")
b_matern32 <- readRDS("toydata_matern32.RDS")

draws_matern32 <- m_matern32$draws()
draws <- b_matern32$draws()
dim(draws)[2]

chains <- list()
for (i in 1:dim(draws)[2]) {
  chains[[i]] <- draws[, i, ]
}

library(purrr)

chains <- purrr::map(chains, as.data.frame)

for (i in 1:length(chains)) {
  chains[[i]] <- chains[[i]][, -grep(pattern="zgp_1",colnames(chains[[i]]))]
  chains[[i]] <- chains[[i]][, -grep(pattern="mu_generated",colnames(chains[[i]]))]
}

chains <- purrr::map(chains, subset, select = 2:5)
purrr::map(chains, head, n = 5)

dta <- do.call(cbind, chains)
tmp <- list()
for (i in 1:(length(colnames(dta))/length(chains))) {
  tmp[[i]] <- dta[,seq(i, ncol(dta), 4)]
}

mlt <- purrr::map(tmp, reshape2::melt)
# function to generate lineplot from such data
df <- mlt[[1]]
plt_chains <- function(df, nwalks = 4, variable = "Intercept") {
  df$rown <- rep(seq(from = 1, to = (nrow(df)/nwalks)), times = 4)
  get_labels <- function(n) {
    labels <- as.character(1:n)
    labels
  }
  plt <- ggplot(df, aes(x = rown, y = value, group = variable, col = variable)) +
    geom_line() +
    scale_color_manual(values = c("#313695", "#74add1", "#a50026", "#f46d43"),
                       labels = get_labels(nwalks)) +
    theme_bw() +
    xlab("Iteration") +
    ylab(variable) +
    guides(color=guide_legend(title="Chain"))
  plt
}

# Intercept, Sigma GP, l-scale, sigma
plt_chains(mlt[[1]], variable = "Intercept")
plt_chains(mlt[[2]], variable = "Sigma GP")
plt_chains(mlt[[3]], variable = "l-scale")
plt_chains(mlt[[4]], variable = "Sigma")

test <- tmp[[1]]
colnames(test) <- c("Intercept", "Intercept", "Intercept", "Intercept")

# getting all relevant values into a matrix
mat <- matrix(nrow = 1600, ncol = 4)
colnames(mat) <- c("Intercept", "SigmaGP", "l-scale", "Sigma")
for (i in 1:length(tmp)) {
  mat[, i] <- reshape2::melt(tmp[[i]])$value
}

# plotting this matrix
bayesplot::mcmc_areas(mat, pars = c("Intercept", "SigmaGP"), prob = 0.9)
bayesplot::mcmc_areas(mat, pars = c("l-scale", "Sigma"), prob = 0.9)

plotchains <- chains
for (i in 1:length(plotchains)) {colnames(plotchains[[i]]) <- c("Intercept", "SigmaGP", "l-scale", "Sigma")}
plotchains <- purrr::map(plotchains, as.matrix.data.frame)
overlay_plt <- bayesplot::mcmc_dens_overlay(plotchains, pars = c("Sigma"))
overlay_plt +
  scale_color_manual(values = c("#313695", "#74add1", "#a50026", "#f46d43"))

# df.new = dta[,seq(1, ncol(dta), 4)]

get_mlt <- function(draws) {
  c1 <- draws[, 1, ]
  c2 <- draws[, 2, ]
  
  c1 <- as.data.frame(c1)
  c2 <- as.data.frame(c2)
  
  c1 <- c1[,-grep(pattern="1.zgp_1",colnames(c1))]
  c2 <- c2[,-grep(pattern="2.zgp_1",colnames(c2))]
  
  c1 <- c1[, 2:5]
  c2 <- c2[, 2:5]
  
  colnames(c1) <- c("Intercept", "SGP", "LScale", "Sigma")
  colnames(c2) <- c("Intercept", "SGP", "LScale", "Sigma")
  
  twochainz <- rbind(c1, c2)
  mlt <- reshape2::melt(twochainz)
  return(mlt)
}

mlt <- get_mlt(draws_matern32)

ggplot(mlt, aes(x = value)) +
  geom_density() +
  facet_wrap(~variable, scales = "free", ncol = 1)

b_matern32$cmdstan_diagnose()

setwd("/home/rstudio/geostatsBA")
file <- paste0(getwd(), "/cmdstanr_matern_52.stan")
mod <- cmdstan_model(file)
############## matern 3 2 for predictive maps
set.seed(1)
dat <- gamSim(1, n = 100)
m_matern32<- gam(y ~ s(x1, x2, bs = "gp", m = 3), data = dat)
summary(m_matern32)
plot(m_matern32)

# stancode <- brms::make_stancode(y ~ gp(x1, x2), data = dat)
standata <- brms::make_standata(y ~ gp(x1, x2), data = dat)

modeldata <- list()
for (i in 1:length(standata)) {
  modeldata[[i]] <- standata[[i]]
}
names(modeldata) <- names(standata)

b_matern_pred <- mod$sample(
  data = modeldata, 
  seed = 123, 
  num_chains = 2, 
  num_cores = 2,
  num_samples = 200,
  num_warmup = 200,
  refresh = 1
)

draws_matern52 <- b_matern52$draws()

mlt_matern_52 <- get_mlt(draws_matern52)

ggplot(mlt_matern_52, aes(x = value)) +
  geom_density() +
  facet_wrap(~variable, scales = "free", ncol = 1)

# modeling the gaussian process smooth here
rm(list = ls())

evidence <- readRDS("Daten/evidence.csv")

library(brms)
# taking a sample of 1000 points
set.seed(1)
evd <- evidence[sample(nrow(evidence), 1000),]

b_smooth <- brm(site ~ s(lon, lat),
                data = evd, 
                family = bernoulli, 
                chains = 4, 
                cores = 4, 
                iter = 1000, 
                control = list(adapt_delta = 0.8, 
                               max_treedepth = 13))
# 4 seconds for 1000 steps

summary(b_smooth)
plot(b_smooth)
plt_cs_500 <- conditional_smooths(b_smooth)
plt_cs_500

b_smooth_fullevidence <- brm(site ~ s(lon, lat),
                             data = evidence, 
                             family = bernoulli, 
                             chains = 4, 
                             cores = 4, 
                             iter = 1000, 
                             control = list(adapt_delta = 0.8, 
                                            max_treedepth = 13))
# about 270 seconds for 1000 steps
# chain 2 clocking in at 5 minutes for 100 iterations
# 50% 12 minutes in
# chain 2: 20.5 minutes total
# chain 3: 1358 seconds
# chain 4: 1377 seconds
# chain 1: 1473 seconds

plt_cs_12222 <- conditional_smooths(b_smooth_fullevidence)
plt_cs_12222

evd_500 <- evd[1:500, ]

b_gp_500 <- brm(site ~ gp(lon, lat),
                data = evd_500, 
                family = bernoulli, 
                chains = 4, 
                cores = 4, 
                iter = 1000, 
                control = list(adapt_delta = 0.8, 
                               max_treedepth = 13))
saveRDS(b_gp_500, "b_gp_500.RDS")

plt_gp_500 <- conditional_effects(b_gp_500)

# starting worker pid=16506 on localhost:11911 at 20:07:26.750
# 1000 transitions using 10 leapfrog steps per transition would take 2871.49 seconds.
# 
# Chain 1:  Elapsed Time: 6477.24 seconds (Warm-up)
# Chain 1:                2909.06 seconds (Sampling)
# Chain 1:                9386.3 seconds (Total)
# Chain 1: 
#   Chain 4: Iteration: 1000 / 1000 [100%]  (Sampling)
# Chain 4: 
#   Chain 4:  Elapsed Time: 6499.09 seconds (Warm-up)
# Chain 4:                2940.27 seconds (Sampling)
# Chain 4:                9439.36 seconds (Total)
# Chain 4: 
#   Chain 3: Iteration: 1000 / 1000 [100%]  (Sampling)
# Chain 3: 
#   Chain 3:  Elapsed Time: 7034.45 seconds (Warm-up)
# Chain 3:                2708.78 seconds (Sampling)
# Chain 3:                9743.24 seconds (Total)
# 

evidence <- readRDS("Daten/evidence.csv")

set.seed(1)
evd <- evidence[sample(nrow(evidence), 1000), ]
evd$site <- as.factor(evd$site)
library(caret)
library(mgcv)
b <- train(site ~ lon + lat + dem + temp + rain + distance_water + 
             frostdays + sunhours + tpi + slope,
           family = binomial,
           data = evd,
           method = "gam", 
           trControl = trainControl(method = "cv", number = 100),
           tuneGrid = data.frame(method = "GCV.Cp", select = TRUE)
)
# r + s + r:s
print(b)
summary(b$finalModel)
