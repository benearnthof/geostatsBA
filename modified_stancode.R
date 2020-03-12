testevidence2 <- sample_n(evidence, size = 250, replace = FALSE)

stancode <- brms::make_stancode(site ~ s(lon, lat, bs = "gp", m = 3) + dem + temp + rain + 
              gp(distance_water) + frostdays + sunhours + tpi + slope,
            family = bernoulli, data = testevidence2, 
            chains = 2, cores = 2, iter = 1000, 
            control = list(adapt_delta = 0.8, 
                           max_treedepth = 13))

standata <- brms::make_stancode(site ~ s(lon, lat, bs = "gp", m = 3) + dem + temp + rain + 
                                  gp(distance_water) + frostdays + sunhours + tpi + slope,
                                family = bernoulli, data = testevidence2, 
                                chains = 2, cores = 2, iter = 1000, 
                                control = list(adapt_delta = 0.8, 
                                               max_treedepth = 13))

# paste matern 32 covariance function in there
require(rstan)
require(StanHeaders)
brm_custom_1d <- rstan::stan_model("gp_matern_32.stan", 
                                   allow_undefined = TRUE,
                                   includes = paste0('\n#include "', file.path(getwd(),'gp_matern_32_cov.hpp'), '"\n'))
