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

# dont forget to save this one
saveRDS(gp_1000, "gp_1000.RDS")
