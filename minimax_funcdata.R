# find rows containing maxima and minima in functional data
# dta are the y values of a functional data frame
test <- crack_simulation(crack_proportion = c(80, 20, 0), crack_fleet = purrr::rdunif(100, 7500, 1), crack_nruns = 200)
dta <- test$ret
dta <- dta[, 2:(ncol(dta) - 1)]
nms <- colnames(dta)[apply(dta, 1, which.max)]
maxi <- unique(nms)
nms <- colnames(dta)[apply(dta, 1, which.min)]
mini <- unique(nms)
extrema <- unique(c(maxi, mini))

# lets see how the number of extrema increases with sample size

results <- numeric(length = 12)
sizes <- c(10, 25, 50, 100, 150, 200, 250, 300, 500, 1000, 2000, 5000)
for (i in 1:12) {
  tmp <- crack_simulation(crack_proportion = c(80, 20, 0), crack_fleet = purrr::rdunif(100, 7500, 1), crack_nruns = sizes[i])
  dta <- tmp$ret[, 2:(ncol(tmp$ret) - 1)]
  nms <- colnames(dta)[apply(dta, 1, which.max)]
  maxi <- unique(nms)
  nms <- colnames(dta)[apply(dta, 1, which.max)]
  mini <- unique(nms)
  extrema <- unique(c(maxi, mini))
  results[i] <- length(extrema)
  print(i)
}
results
# the number of

# from this, the question arises: Is there a bound for the number of rows containing
# extreme values of a random walk? We need some sensible rule to judge the number
# of columns that are returned in the final data frame. We want to convey as little
# information as possible without leaving out important data.
# Let's start by coding up a more efficient simulation

fnc <- function(n) sample(c(1L, -1L), n, replace = TRUE)
cumsum(fnc(1000L))

results <- numeric(length = 1000)
sizes <- c(10, 25, 50, 100, 150, 200, 250, 300, 500, 1000, 2000, 5000, 10000, 25000, 50000)
n <- 1000L
reps <- 25
totalresults <- matrix(ncol = 1000, nrow = reps)

for (k in 1:reps) {
  for (i in 1:n) {
    tmp <- matrix(0, nrow = n, ncol = i)
    colnames(tmp) <- 1:i
    for (j in 1:ncol(tmp)) {
      tmp[, j] <- cumsum(fnc(n))
    }
    maxi <- unique(colnames(tmp)[apply(tmp, 1, which.max)])
    mini <- unique(colnames(tmp)[apply(tmp, 1, which.min)])
    extrema <- unique(c(maxi, mini))
    results[i] <- length(extrema)
    print(i)
  }
  totalresults[k, ] <- results
  print(k)
}

results1k <- results
results100 <- results

colnames(totalresults) <- 1:1000
rownames(totalresults) <- 1:25

dta <- t(totalresults)

dta <- as.data.frame.matrix(dta)
dta$y <- 1:1000

plot.new()
index <- 1:1000
for (i in 1:nrow(totalresults)) {
  plot(totalresults[i,] ~ index, add = TRUE)
}

wot <- lm(y~., data = dta)
summary(wot)
clmeans <- colMeans(totalresults)
dta <- data.frame(y = 1:1000, x = clmeans)

wot <- gam(x~s(y), data = dta)
plot(wot)
summary(wot)

predict.gam(wot, newdata = list(y = 500))
plot(clmeans~index)
points(x = 500, y = 48.7, col = "red")
pdata <- predict.gam(wot, newdata = list(y = 1:1000))
points(x = 1:1000, y = pdata, col = "red")
points(x = 1:1000, y = 2*sqrt(x), col = "blue")
points(x = 1:1000, y = 8*log(x), col = "green")
points(x = 1:1000, y = exp(1) * 3 * log(x), col = "magenta")
points(x = 1:1000, y = log(1000) * log(x))
# it seems we can use this model to predict the amount of random walks that hit
# a maximum or minimum
# this may not be the ideal or most elegant way to go about solvint this problem
# but i think a simple decision rule like the following would be suited to tackle
# the issue at hand
# if nruns <= 10 return all runs
# if nruns > 10 <= 500 return all runs that yielded extreme values
# if nruns > 500 return the first 50 runs that yielded extreme values
# in any case we should also return the average of all runs.

# another model that seems to fit quite well is 8*log(x) or exp(1) * 3 * log(x)
# but they overestimate the amount of extreme paths for the lower amount of paths.

# what i forgot to mention, is that the number of paths that reach extreme values
# is, of course, dependent on the length of the random walks aswell.
# Given infinite length one might expect every single path to reach extreme values
# infinitely many times.