# installing cmdstanr to use covariance functions
devtools::install_github("stan-dev/cmdstanr")
library(cmdstanr)
# requires working version of cmdstan
install_cmdstan()
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

