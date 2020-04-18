# installing cmdstanr to use covariance functions
devtools::install_github("stan-dev/cmdstanr")
library(cmdstanr)
# requires working version of cmdstan
install_cmdstan()
# verify that version can be found
cmdstan_version()

