# advanced packageloading utility for the bioconductor ecosystem 
# as provided @ https://www.huber.embl.de/msmb/install_packages.R
options(install.packages.check.source = "no")
options(install.packages.compile.from.source = "never")
Sys.setenv(R_REMOTES_UPGRADE = "never")

## Function to install packages one at a time with indication of time left
## Overall probably slower than install.packages if everything works
## but doesn't require downloading all packages first before trying to install any
installer_with_progress <- function(pkgs) {
  
  if(length(pkgs) == 0) { invisible(return(NULL)) }
  
  toInstall <- pkgs
  bp <- progress::progress_bar$new(total = length(toInstall),
                                   format = "Installed :current of :total (:percent ) - current package: :package",
                                   show_after = 0,
                                   clear = FALSE)
  
  length_prev <- length(toInstall)
  fail <- NULL
  while(length(toInstall)) {
    pkg <- toInstall[1]
    bp$tick(length_prev - length(toInstall),  tokens = list(package = pkg))
    length_prev <- length(toInstall)
    tryCatch(
      suppressMessages( BiocManager::install(pkg, quiet = TRUE, update = FALSE, ask = FALSE, type = "binary") ),
      error = function(e) { fail <<- c(fail, pkg) },
      warning = function(w) { fail <<- c(fail, pkg) },
      ## remove current package, otherwise we loop in event of failure
      ## update the list to reflect any dependencies that are now installed
      finally = { toInstall <- setdiff(toInstall, installed.packages()[, "Package"]) }
    )
  }
  bp$tick(length_prev - length(toInstall),  tokens = list(package = "DONE!"))
  
  return(fail)
}

## these packages are needed prior to the installation
if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(c('BiocManager'), repos = "https://cloud.r-project.org",
                   quiet = TRUE, update = FALSE, ask = FALSE, type = "both")
}
## update any existing packages
BiocManager::install(update = TRUE, ask = FALSE)

if(!requireNamespace("remotes", quietly = TRUE)) {
  install.packages(c('remotes'), quiet = TRUE, update = FALSE, ask = FALSE, type = "both")
}
if(!requireNamespace("magrittr", quietly = TRUE)) {
  BiocManager::install('magrittr', quiet = TRUE, update = FALSE, ask = FALSE, type = "both")
}
if(!requireNamespace("progress", quietly = TRUE)) {
  BiocManager::install('progress', quiet = TRUE, update = FALSE, ask = FALSE, type = "both")
}

## structSSI is currently deprecated and has been removed from CRAN for now (24-06-2020)
## This will install a CRAN version by default if it reappears, otherwise use an archive version
## Update 17-05-2021: This isn't coming back to CRAN any time soon, so lets use the GitHub version
if(!requireNamespace("structSSI", quietly = TRUE)) {
  BiocManager::install('krisrs1128/structSSI', upgrade = FALSE, quiet = TRUE, ask = FALSE, type = "both")
}

## list of packages required for each chapters
chapter_pkgs <- readRDS(url("https://www.huber.embl.de/msmb/chapter_pkgs.rds"))

## subset a selection of chapters if specified
if(exists('chapter_index') && is.numeric(chapter_index)) {
  chapter_pkgs <- chapter_pkgs[ chapter_index ]
}

for(i in seq_along(chapter_pkgs)) {
  message("### CHAPTER: ", i, " ###")
  pkgsAvailable  = installed.packages()[, "Package"]
  pkgsToInstall = setdiff(chapter_pkgs[[i]], pkgsAvailable)
  BiocManager::install(pkgsToInstall, update = FALSE, upgrade = FALSE, ask = FALSE, type = "both")
}

## report packages no installed
## find only those not currently installed
pkgsAvailable  = installed.packages()[, "Package"]
pkgsNeeded = unique(unlist(chapter_pkgs))
pkgsToInstall = setdiff(pkgsNeeded, pkgsAvailable)
if(length(pkgsToInstall)) {
  message("The following packages failed to install: \n",
          paste(pkgsToInstall, collapse = ", "))
  message("You can try re-running this installation script.\n",
  "It will only try to install the missing packages.\n",
  "This may make it easier to see the information R gives about why the installation failed.\n",
  "Please contact mike.smith@embl.de if you need additional help.")
}

Sys.unsetenv("R_REMOTES_UPGRADE")
