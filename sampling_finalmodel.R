# final modeling in brms and cmdstanr
library(mgcv)
library(brms)
library(cmdstanr)
evidence <- readRDS("Daten/evidence.csv")
expo <- gam(site ~ s(lon, lat, bs = "gp", m = 2) + s(dem) + temp + s(rain) + 
              distance_water + sunhours + s(tpi) + slope,
            family = binomial, data = evidence, select = TRUE, method = "REML")

# change link to bernoulli for more efficient sampling
# doing approximate latent gaussian process to speed up sampling by a lot
stancode <- brms::make_stancode(site ~ gp(lon, lat, k = 30, c = 5/4) + s(dem) + temp + s(rain) + 
                                  distance_water + sunhours + s(tpi) + slope,
                                family = bernoulli, data = evidence)

standata <- brms::make_standata(site ~ gp(lon, lat, k = 30, c = 5/4) + s(dem) + temp + s(rain) + 
                                  distance_water + sunhours + s(tpi) + slope,
                                family = bernoulli, data = evidence)

set.seed(1)
evd_200 <- evidence[sample(nrow(evidence), size = 200),]
# can we get away with max treedepth of 10?
# no we cannot
b_expo_200 <- brms::brm(site ~ gp(lon, lat, k = 30, c = 5/4) + s(dem) + temp + s(rain) + 
                          distance_water + sunhours + s(tpi) + slope,
                        family = bernoulli,
                        data = evd_200,
                        chains = 4, 
                        cores = 4, 
                        iter = 200, 
                        control = list(adapt_delta = 0.8, 
                                       max_treedepth = 10)
                        )
saveRDS(b_expo_200, "final_expo_200.RDS")
# times: 
# 20 sec for 1000 iterations
# longest chain took 158 seconds

set.seed(1)
evd_400 <- evidence[sample(nrow(evidence), size = 400),]

b_expo_400 <- brms::brm(site ~ gp(lon, lat, k = 30, c = 5/4) + s(dem) + temp + s(rain) + 
                          distance_water + sunhours + s(tpi) + slope,
                        family = bernoulli,
                        data = evd_400,
                        chains = 4, 
                        cores = 4, 
                        iter = 400, 
                        control = list(adapt_delta = 0.8, 
                                       max_treedepth = 12)
)
saveRDS(b_expo_400, "final_expo_400.RDS")
# times: 17, 14, 52, 52 seconds for 1k steps each, huh server seems to be relatively busy
# elapsed time ~ 22 minutes
# max_treedepth 12 seems to be the way to go
# lets see what duration we get for 1k points

set.seed(1)
evd_1000 <- evidence[sample(nrow(evidence), size = 1000),]

b_expo_1000 <- brms::brm(site ~ gp(lon, lat, k = 30, c = 5/4) + s(dem) + temp + s(rain) + 
                          distance_water + sunhours + s(tpi) + slope,
                        family = bernoulli,
                        data = evd_1000,
                        chains = 4, 
                        cores = 4, 
                        iter = 400, 
                        control = list(adapt_delta = 0.8, 
                                       max_treedepth = 12)
)
saveRDS(b_expo_1000, "final_expo_1000.RDS")
# about 120 seconds for 1000 transitions, not too bad considering we more than 
# doubled the amount of data
# started at 22:40:45
# 40 warmup iterations clocking in 10 minutes later at 22:51
# 50% at 23:54

set.seed(1)
b_expo_12222 <- brms::brm(site ~ gp(lon, lat, k = 30, c = 5/4) + s(dem) + temp + s(rain) + 
                           distance_water + sunhours + s(tpi) + slope,
                         family = bernoulli,
                         data = evidence,
                         chains = 4, 
                         cores = 4, 
                         iter = 1000, 
                         control = list(adapt_delta = 0.8, 
                                        max_treedepth = 12)
)
saveRDS(b_expo_12222, "final_expo_12222.RDS")
# model seems to be severely misspecified, slow mixing chains etc.
# 1400 secs for 1000 transitions
# started at 23:20:12.340
# only chain 3 seems to do sampling properly





set.seed(1)
evd_1000 <- evidence[sample(nrow(evidence), size = 1000),]

standata <- brms::make_standata(site ~ gp(lon, lat),
                      family = bernoulli,
                      data = evd_1000,
                      chains = 4, 
                      cores = 4, 
                      iter = 1000, 
                      control = list(adapt_delta = 0.8, 
                                     max_treedepth = 12))
modeldata <- list()

for (i in 1:length(standata)) {
  modeldata[[i]] <- standata[[i]]
}

names(modeldata) <- names(standata)

library(cmdstanr)
file <- paste0(getwd(), "/finalmodel_matern32.stan")
mod <- cmdstan_model(file)

b_matern32_12222 <- mod$sample(
  data = modeldata, 
  seed = 123, 
  num_chains = 4, 
  num_cores = 4,
  num_samples = 500,
  num_warmup = 500,
  adapt_delta = 0.8,
  max_depth = 12,
  refresh = 10
)
saveRDS(b_matern32_12222, file = "b_matern32_12222.RDS")
# this one did not converge aswell

# matern52
file <- paste0(getwd(), "/finalmodel_matern52.stan")
mod_2 <- cmdstan_model(file)

b_matern52_12222 <- mod_2$sample(
  data = modeldata, 
  seed = 123, 
  num_chains = 4, 
  num_cores = 4,
  num_samples = 500,
  num_warmup = 500,
  adapt_delta = 0.8,
  max_depth = 12,
  refresh = 10
)
saveRDS(b_matern52_12222, file = "b_matern52_12222.RDS")
# 15% after 10 hours

# matern 52 model does not run, seems to be misspecified
set.seed(1)
evd_3000 <- evidence[sample(nrow(evidence), size = 3000),]
standata_matern52 <- brms::make_standata(site ~ gp(lon, lat),
                                family = bernoulli,
                                data = evd_3000,
                                chains = 4, 
                                cores = 4, 
                                iter = 1000, 
                                control = list(adapt_delta = 0.8, 
                                               max_treedepth = 12))
stancode_matern52 <- brms::make_stancode(site ~ gp(lon, lat),
                                         family = bernoulli,
                                         data = evd_3000,
                                         chains = 4, 
                                         cores = 4, 
                                         iter = 1000, 
                                         control = list(adapt_delta = 0.8, 
                                                        max_treedepth = 12))


# matern 32 without smooth terms
set.seed(1)
evd <- evidence[sample(nrow(evidence), size = 400),]
stancode_matern32 <- brms::make_stancode(site ~ gp(lon, lat),
                                         family = bernoulli,
                                         data = evd,
                                         chains = 4, 
                                         cores = 4, 
                                         iter = 1000, 
                                         control = list(adapt_delta = 0.8, 
                                                        max_treedepth = 12))
standata <- brms::make_standata(site ~ gp(lon, lat),
                      family = bernoulli,
                      data = evd,
                      chains = 4, 
                      cores = 4, 
                      iter = 1000, 
                      control = list(adapt_delta = 0.8, 
                                     max_treedepth = 12))
modeldata <- list()

for (i in 1:length(standata)) {
  modeldata[[i]] <- standata[[i]]
}

names(modeldata) <- names(standata)

library(cmdstanr)
file <- paste0(getwd(), "/cmdstan_matern32_1000.stan")
mod <- cmdstan_model(file)

b_matern32_1000_novariables <- mod$sample(
  data = modeldata, 
  seed = 123, 
  num_chains = 4, 
  num_cores = 4,
  num_samples = 200,
  num_warmup = 200,
  adapt_delta = 0.8,
  max_depth = 12,
  refresh = 10
)
saveRDS(b_matern32_1000_novariables, file = "b_matern32_1000_novariables.RDS")
b_matern32_1000_novariables$output_files()
b_matern32_1000_novariables$cmdstan_diagnose()
b_matern32_1000_novariables$summary()



# doing another run without smooths and a fifth of the data
set.seed(1)
evd_3000 <- evidence[sample(nrow(evidence), size = 3000),]
b_expo_3000 <- brms::brm(site ~ gp(lon, lat, k = 30, c = 5/4) + distance_water + dem,
                          family = bernoulli,
                          data = evd_3000,
                          chains = 4, 
                          cores = 4, 
                          iter = 1000, 
                          control = list(adapt_delta = 0.8, 
                                         max_treedepth = 12))
saveRDS(b_expo_3000, file = "b_expo_3000_nosmooths.RDS")
# saveRDS(b_expo_5000, file = "b_expo_2000_nosmooths.RDS")
# started sampling 23/04/20 at 21:05

plot(b_expo_5000)
tmp <- readRDS(file = "b_expo_2000_nosmooths.RDS")

######## running matern 52 without covariables
modeldata <- list()

for (i in 1:length(standata_matern52)) {
  modeldata[[i]] <- standata_matern52[[i]]
}

names(modeldata) <- names(standata_matern52)

library(cmdstanr)
file <- paste0(getwd(), "/finalmodel_matern52_novariables.stan")
mod <- cmdstan_model(file)

b_matern52_3000_novariables <- mod$sample(
  data = modeldata, 
  seed = 123, 
  num_chains = 4, 
  num_cores = 4,
  num_samples = 500,
  num_warmup = 500,
  adapt_delta = 0.8,
  max_depth = 12,
  refresh = 10
)
saveRDS(b_matern52_3000_novariables, file = "b_matern52_3000_novariables.RDS")


set.seed(1)
evd_3000 <- evidence[sample(nrow(evidence), size = 3000),]
b_expo_3000_2 <- brms::brm(site ~ gp(lon, lat, k = 30, c = 5/4) + distance_water + dem +
                           temp + rain,
                         family = bernoulli,
                         data = evd_3000,
                         chains = 4, 
                         cores = 4, 
                         iter = 1000, 
                         control = list(adapt_delta = 0.8, 
                                        max_treedepth = 12))
saveRDS(b_expo_3000_2, file = "b_expo_3000_nosmooths_2.RDS")


set.seed(1)
evd_2000 <- evidence[sample(nrow(evidence), size = 2000),]
b_expo_2000_iso <- brms::brm(site ~ gp(lon, lat, k = 30, c = 5/4),
                           family = bernoulli,
                           data = evd_2000,
                           chains = 4, 
                           cores = 4, 
                           iter = 1000, 
                           control = list(adapt_delta = 0.8, 
                                          max_treedepth = 12))
saveRDS(b_expo_2000_iso, file = "b_expo_2000_iso.RDS")



saveRDS(b_matern32_1000_novariables, file = "bm321knovar.RDS") 
b_matern32_1000_novariables$summary()
