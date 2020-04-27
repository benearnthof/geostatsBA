evidence <- readRDS("Daten/evidence.csv")

library(brms)
library(mgcv)

set.seed(1)
evd <- evidence[sample(nrow(evidence), 1000),]
nrow(evd)

gp_1000 <- brms::brm(site ~ gp(lon, lat), data = evd, 
                  family = bernoulli,
                  chains = 4, 
                  cores = 4, 
                  iter = 1000, 
                  control = list(adapt_delta = 0.8, 
                                 max_treedepth = 13))

# longest gradient evaluation for 1000 steps: 2999.72 seconds
# longest sampling time for single chain
# Chain 3:  Elapsed Time: 8753.02 seconds (Warm-up)
# Chain 3:                3365.59 seconds (Sampling)
# Chain 3:                12118.6 seconds (Total)
# dont forget to save this one
saveRDS(gp_1000, "gp_1000.RDS")

# lets try approximate gp with full data set
set.seed(1)

evd <- evidence
nrow(evd)
gp_1000_approx_30 <- brms::brm(site ~ gp(lon, lat, k = 30, c = 5/4), data = evd, 
                     family = bernoulli,
                     chains = 4, 
                     cores = 4, 
                     iter = 1000, 
                     control = list(adapt_delta = 0.8, 
                                    max_treedepth = 13))
saveRDS(gp_1000_approx_30, file = "gp_12222_approx_30.RDS")
# starting worker pid=25332 on localhost:11090 at 16:23:24.116
# Chain 3: 1000 transitions using 10 leapfrog steps per 
# transition would take 1252.35 seconds.
# first 100 iterations in chain 2 at 18:16 ~2 hours for 100 iterations not bad
# chain 3 and 4 clocking in 100 iters at 18:26
# chain 1 comin in last at 18:40
# fastest chain 60k seconds

set.seed(1)
# running model without smooth terms to see if that alleviates the problem 
# of very slow convergence
b_expo_12222_nosmooths <- brms::brm(site ~ gp(lon, lat, k = 30, c = 5/4) + dem + temp + rain + 
                            distance_water + sunhours + tpi + slope,
                          family = bernoulli,
                          data = evidence,
                          chains = 4, 
                          cores = 4, 
                          iter = 1000, 
                          control = list(adapt_delta = 0.8, 
                                         max_treedepth = 12))
saveRDS(b_expo_12222_nosmooths, "b_expo_12222_nosmooths.RDS")
fit <- b_expo_12222_nosmooths$fit
str(fit@sim$samples)
c2 <- fit@sim$samples[[2]]
length(c2)
