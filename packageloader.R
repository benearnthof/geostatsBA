# packages to be loaded
packages <- c(
  "broom", "sf", "geosphere", "raster", "dplyr", "spData", "remotes", "mlr",
  "vegan", "mapview", "ggplot2", "shiny", "zoo", "rbenchmark", "tmaptools",
  "shinyjs", "parallelMap", "rasterVis", "arm", "latticeExtra", "grid", "pROC",
  "maptools", "doParallel", "shinythemes", "devtools", "mgcv", "brms", "gstat",
  "lattice", "nlme", "spdep", "schoenberg"
)

if (!require(schoenberg)) devtools::install_github('gavinsimpson/schoenberg')

# checking which packages have been installed
packs <- lapply(packages, FUN = function(packages) {
  do.call("require", list(packages))
})
packs <- unlist(packs, use.names = F)
# generating list of packages yet to be installed
instpacks <- packages[!packs]

# installing all packages that have not yet been installed
lapply(instpacks, FUN = function(instpacks) {
  do.call("install.packages", list(instpacks))
})

# should return a vector of TRUE entries - one entry for every successfully loaded package
check <- unlist(lapply(packages, FUN = function(packages) {
  do.call("require", list(packages))
}))

failed <- which(check == FALSE)
failed <- packages[failed]

if (identical(character(0), failed)) {
  print("All packages loaded successfully.")
} else {
  cat("Packages", "\n", failed, "\n", "could not be loaded.")
}
