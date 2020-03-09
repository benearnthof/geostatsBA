# lets try implementing covariance functions for brms
set.seed(1)
dat <- mgcv::gamSim(1, n = 250, scale = 2)

fit_spherical <- gam(y ~ s(x1, x2, bs = "gp", m = 1), data = dat)
fit_exponential <- gam(y ~ s(x1, x2, bs = "gp", m = 2), data = dat)
fit_matern_1_5 <- gam(y ~ s(x1, x2, bs = "gp", m = 3), data = dat)
fit_matern_2_5 <- gam(y ~ s(x1, x2, bs = "gp", m = 4), data = dat)
fit_matern_3_5 <- gam(y ~ s(x1, x2, bs = "gp", m = 5), data = dat)

library(mgcViz)
plot(getViz(fit_spherical))
plot(getViz(fit_exponential))
plot(getViz(fit_matern_1_5))
plot(getViz(fit_matern_2_5))
plot(getViz(fit_matern_3_5))

viz <- getViz(fit_matern_1_5)

# fit and plot a two-dimensional smooth term
fit2 <- brm(y ~ gp(x1, x2), data = dat, chains = 2, cores = 2, iter = 1000)
ms <- marginal_smooths(fit2)
plot(ms, stype = "contour")
plot(ms, stype = "raster")

ms2 <- marginal_effects(fit2, nsamples = 200, spaghetti = TRUE)
plot(ms2, ask = FALSE, points = TRUE)

smoothobject <- s(dat$x1, dat$x2, bs = "gp", m = 3)
# smooth variables from debugging function
gpsmooth <- smooth.construct.gp.smooth.spec(smoothobject, data = smoothdata, knots = smoothknots)
gpsmooth$gp.defn

# this is the way it is implemented in mgcv
gpT <- function(x) {
  ## T matrix for Kamman and Wand Matern Spline...
  cbind(x[,1]*0+1,x)
} ## gpT

gpE <- function(x,xk,defn = NA) {
  ## Get the E matrix for a Kammann and Wand Matern spline.
  ## rho is the range parameter... set to K&W default if not supplied
  ind <- expand.grid(x=1:nrow(x),xk=1:nrow(xk))
  ## get d[i,j] the Euclidian distance from x[i] to xk[j]... 
  E <- matrix(sqrt(rowSums((x[ind$x,,drop=FALSE]-xk[ind$xk,,drop=FALSE])^2)),nrow(x),nrow(xk))
  rho <- -1; k <- 1
  if ((length(defn)==1&&is.na(defn))||length(defn)<1) { type <- 3 } else
    if (length(defn)>0) type <- round(defn[1])
  if (length(defn)>1) rho <- defn[2]
  if (length(defn)>2) k <- defn[3]
  
  if (rho <= 0) rho <- max(E) ## approximately the K & W choise
  E <- E/rho
  if (!type%in%1:5||k>2||k<=0) stop("incorrect arguments to GP smoother")
  if (type>2) eE <- exp(-E)
  E <- switch(type,
              (1 - 1.5*E + 0.5 *E^3)*(E <= 1), ## 1 spherical 
              exp(-E^k), ## 2 power exponential
              (1 + E) * eE, ## 3 Matern k = 1.5
              eE + (E*eE)*(1+E/3), ## 4 Matern k = 2.5
              eE + (E*eE)*(1+.4*E+E^2/15) ## 5 Matern k = 3.5
  )
  attr(E,"defn") <- c(type,rho,k)
  E
} ## gpE

smooth.construct.gp.smooth.spec <- function(object,data,knots)
  ## The constructor for a Kamman and Wand (2003) Matern Spline, and other GP smoothers.
  ## See also Handcock, Meier and Nychka (1994), and Handcock and Stein (1993).
{ ## deal with possible extra arguments of "gp" type smooth
  xtra <- list()
  
  if (is.null(object$xt$max.knots)) xtra$max.knots <- 2000 
  else xtra$max.knots <- object$xt$max.knots 
  if (is.null(object$xt$seed)) xtra$seed <- 1 
  else xtra$seed <- object$xt$seed 
  
  ## now collect predictors
  x <- array(0,0)
  
  for (i in 1:object$dim) {
    xx <- data[[object$term[i]]]
    if (i==1) n <- length(xx) else 
      if (n!=length(xx)) stop("arguments of smooth not same dimension")
    x<-c(x,xx)
  }
  
  if (is.null(knots)) { knt <- 0; nk <- 0}
  else { 
    knt <- array(0,0)
    for (i in 1:object$dim) { 
      dum <- knots[[object$term[i]]]
      if (is.null(dum)) { knt <- 0; nk <- 0; break} # no valid knots for this term
      knt <- c(knt,dum)
      nk0 <- length(dum)
      if (i > 1 && nk != nk0) 
        stop("components of knots relating to a single smooth must be of same length")
      nk <- nk0
    }
  }
  if (nk>n) { ## more knots than data - silly.
    nk <- 0
    warning("more knots than data in an ms term: knots ignored.")
  }
  
  xu <- uniquecombs(matrix(x,n,object$dim),TRUE) ## find the unique `locations'
  if (nrow(xu) < object$bs.dim) stop(
    "A term has fewer unique covariate combinations than specified maximum degrees of freedom")
  ## deal with possibility of large data set
  if (nk==0) { ## need to create knots
    nu <- nrow(xu)  ## number of unique locations
    if (n > xtra$max.knots) { ## then there *may* be too many data      
      if (nu > xtra$max.knots) { ## then there is really a problem 
        seed <- try(get(".Random.seed",envir=.GlobalEnv),silent=TRUE) ## store RNG seed
        if (inherits(seed,"try-error")) {
          runif(1)
          seed <- get(".Random.seed",envir=.GlobalEnv)
        }
        kind <- RNGkind(NULL)
        RNGkind("default","default")
        set.seed(xtra$seed) ## ensure repeatability
        nk <- xtra$max.knots ## going to create nk knots
        ind <- sample(1:nu,nk,replace=FALSE)  ## by sampling these rows from xu
        knt <- as.numeric(xu[ind,])  ## ... like this
        RNGkind(kind[1],kind[2])
        assign(".Random.seed",seed,envir=.GlobalEnv) ## RNG behaves as if it had not been used
      } else { 
        knt <- xu; nk <- nu
      } ## end of large data set handling
    } else { knt <- xu;nk <- nu } ## just set knots to data
  } 
  
  x <- matrix(x,n,object$dim)
  knt <- matrix(knt,nk,object$dim)
  
  ## centre the covariates...
  
  object$shift <- colMeans(x)
  x <- sweep(x,2,object$shift)
  knt <- sweep(knt,2,object$shift)
  
  ## Get the E matrix...
  E <- gpE(knt,knt,object$p.order)
  object$gp.defn <- attr(E,"defn")
  
  def.k <- c(10,30,100)
  dd <- ncol(knt)
  if (object$bs.dim[1] < 0) object$bs.dim <- ncol(knt) + 1 + def.k[dd] ## default basis dimension 
  if (object$bs.dim < ncol(knt)+2) {
    object$bs.dim <- ncol(knt)+2
    warning("basis dimension reset to minimum possible")
  }
  object$null.space.dim <- ncol(knt) + 1
  
  k <- object$bs.dim - object$null.space.dim   
  
  if (k < nk) {
    er <- slanczos(E,k,-1) ## truncated eigen decomposition of E
    D <- diag(c(er$values,rep(0,object$null.space.dim))) ## penalty matrix
  } else { ## no point using eigen-decomp
    D <- matrix(0,object$bs.dim,object$bs.dim)
    D[1:k,1:k] <- E  ## penalty
    er <- list(vectors=diag(k)) ## U is identity here
  }
  rm(E)
  
  object$S <- list(S=D)
  
  object$UZ <- er$vectors ## UZ - (original params) = UZ %*% (working params)
  
  object$knt = knt ## save the knots
  object$df <- object$bs.dim
  object$rank <- k
  class(object)<-"gp.smooth"
  
  object$X <- Predict.matrix.gp.smooth(object,data)
  
  object
} ## end of smooth.construct.gp.smooth.spec

# lets find out what the cov_exp_quad function in brms outputs
library(brms)
debug(brm)
fit_brm_expo <- brm(y ~ gp(x1, x2), data = dat, chains = 2, cores = 2, iter = 1000)
plot(fit_brm_expo)

stancode <- make_stancode(y ~ gp(x1, x2), data = dat, chains = 2, cores = 2, iter = 1000)

# we need to expand the cov_exp_quad function
# we need to expand the spd_cov_exp_quad function
# and add the call to the gaussian process computation 

# exponential-quadratic covariance matrix
# not vectorized over parameter values
cov_exp_quad <- function(x, x_new = NULL, sdgp = 1, lscale = 1) {
  Dls <- length(lscale)
  if (Dls == 1L) {
    # one dimensional or isotropic GP
    diff_quad <- brms:::diff_quad(x = x, x_new = x_new)
    out <- sdgp^2 * exp(-diff_quad / (2 * lscale^2))
  } else {
    # multi-dimensional non-isotropic GP
    diff_quad <- brms:::diff_quad(x = x[, 1], x_new = x_new[, 1])
    out <- sdgp^2 * exp(-diff_quad / (2 * lscale[1]^2))
    for (d in seq_len(Dls)[-1]) {
      diff_quad <- brms:::diff_quad(x = x[, d], x_new = x_new[, d])
      out <- out * exp(-diff_quad / (2 * lscale[d]^2))
    }
  }
  out
}

# spectral density function
# vectorized over parameter values
spd_cov_exp_quad <- function(x, sdgp = 1, lscale = 1) {
  NB <- NROW(x)
  D <- NCOL(x)
  Dls <- NCOL(lscale)
  out <- matrix(nrow = length(sdgp), ncol = NB)
  if (Dls == 1L) {
    # one dimensional or isotropic GP
    constant <- sdgp^2 * (sqrt(2 * pi) * lscale)^D
    neg_half_lscale2 <- -0.5 * lscale^2
    for (m in seq_len(NB)) {
      out[, m] <- constant * exp(neg_half_lscale2 * sum(x[m, ]^2))
    }
  } else {
    # multi-dimensional non-isotropic GP
    constant <- sdgp^2 * sqrt(2 * pi)^D * matrixStats::rowProds(lscale)
    neg_half_lscale2 = -0.5 * lscale^2
    for (m in seq_len(NB)) {
      x2 <- as_draws_matrix(x[m, ]^2, dim = dim(lscale))
      out[, m] <- constant * exp(rowSums(neg_half_lscale2 * x2))
    }
  }
  out
}

# compute the mth eigen value of an approximate GP
eigen_val_cov_exp_quad <- function(m, L) {
  ((m * pi) / (2 * L))^2
}

# compute the mth eigen function of an approximate GP
eigen_fun_cov_exp_quad <- function(x, m, L) {
  x <- as.matrix(x)
  D <- ncol(x)
  stopifnot(length(m) == D, length(L) == D)
  out <- vector("list", D)
  for (i in seq_cols(x)) {
    out[[i]] <- 1 / sqrt(L[i]) * 
      sin((m[i] * pi) / (2 * L[i]) * (x[, i] + L[i]))
  }
  Reduce("*", out)
}

fit_exponential_1d <- gam(y ~ s(x2, bs = "gp", m = 2), data = dat)
brm_exponential_1d <- brm(y ~ gp(x2), data = dat, chains = 2, cores = 2, iter = 1000)

plot(fit_exponential_1d)
plot(brm_exponential_1d)
me <- marginal_effects(brm_exponential_1d, nsamples = 200, spaghetti = TRUE, nug = 0.1)
plot(me)

fit_matern_1_5_1d <- gam(y ~ s(x2, bs = "gp", m = 3), data = dat)

plot(fit_matern_1_5_1d)

# spectral density of matern covariance process
# n = dimension 
# pi = pi
# nu = half integer
# rho = lscale = l = covariance parameter
# f = s = value of x

# spectralGP for spectral density function needed for matern covariance
install.packages("spectralGP")
library(spectralGP)

debug(smooth.construct.gp.smooth.spec)
fit_matern_1_5_1d <- gam(y ~ s(x2, bs = "gp", m = 3), data = dat)
# stored E and knt through debug
E
# E is the matern covariance matrix obtained through gam
knt
?matern.specdens()
# omega is two column matrix or a vector. 
# we used a vector to calculate covariances in the 1D example
# here it is just dat$x2
# param => rho and k, range and differentiability parameter. 
# here given by 
attributes(E)$defn[2:3]
# d = dimension of the domain
matern.specdens(dat$x2, param = list(0.991199, 1), d = 1)
debug(matern.specdens)
undebug(matern.specdens)
matern.specdens(dat$x2, param = c(0.991199, 1), d = 1)

# fields for matern to calculate matern covariance quickly
install.packages("fields")
library(fields)

?Matern
# returns a vector of covariances.

exp_brm <- cov_exp_quad(dat$x2)
exp_fields <- Exp.simple.cov(dat$x2)

exp_cov_simple <- function(x, x_new = NULL, sdgp = 1, lscale = 1, C = NA, marginal = FALSE) {
  if (!is.null(x_new)) {
    x_new <- x
  }
  if (is.na(C[1]) & !marginal) {
    return(sdgp^2 * exp(-brms:::diff_quad(x, x_new)/(2 * lscale)))
  }
  if (!is.na(C[1])) {
    return(sdgp^2 * exp(-brms:::diff_quad(x, x_new)/(2 * lscale)) %*% C)
  }
  if (marginal) {
    return(rep(1, nrow(x)))
  }
}

exp_mix <- exp_cov_simple(dat$x2)
all.equal(exp_brm, exp_mix)

# lets see if we can modify the stancode with the "new" covariance function
# models
fit_exponential_1d <- gam(y ~ s(x2, bs = "gp", m = 2), data = dat)
set.seed(1)
brm_exponential_1d <- brm(y ~ gp(x2), data = dat, chains = 2, cores = 2, iter = 1000)
brm_exponential_1d_spline <- brm(y ~ s(x2), data = dat, chains = 2, cores = 2, iter = 1000)


# now lets modify stancode 
stancode <- make_stancode(y ~ gp(x2), data = dat, chains = 2, cores = 2, iter = 1000)
stancode2 <- make_stancode(y ~ s(x2), data = dat, chains = 2, cores = 2, iter = 1000)
standata <- make_standata(y ~ gp(x2), data = dat, chains = 2, cores = 2, iter = 1000)
standata2 <- make_standata(y ~ s(x2), data = dat, chains = 2, cores = 2, iter = 1000)

require(rstan)
require(StanHeaders)
brm_custom_1d <- stan_model("onedimGP.stan", allow_undefined = TRUE,
                            includes = paste0('\n#include "', file.path(getwd(),'gp_exponential_cov.hpp'), '"\n'))
brm_custom_1d <- stan_model("onedimGP.stan")
set.seed(1)
fit <- sampling(brm_custom_1d, data = standata, iter = 1000, chains = 2, cores = 2)

set.seed(1)
brm_spline <- stan_model("onedimspline.stan")
brm_spline_fit <- sampling(brm_spline, data = standata2)

# we need to construct a full brms object to be able to plot marginal smooths
# step into the brm function to see how the object is being constructed.

plot(fit_exponential_1d)
plot(brm_exponential_1d)
me <- marginal_effects(brm_exponential_1d, nsamples = 200, spaghetti = TRUE, nug = 0.1)
plot(me)

# how to proceeed from here: 
# generate stancode for models from formula
# enter functions from the math library to calculate different covariance matrices
# compile that to a stanmodel
# convert that to a brms object
# use marginal_smooths for smooth plots.
# y ~ s(x2), data = dat, chains = 2, cores = 2, iter = 1000
formula <- "y ~ s(x2)"
data = dat 
family = gaussian() 
prior = NULL
autocor = NULL 
data2 = NULL 
cov_ranef = NULL
sample_prior = "no"
sparse = NULL
knots = NULL
stanvars = NULL
stan_funs = NULL
fit = NA
save_ranef = TRUE
save_mevars = FALSE
save_all_pars = FALSE
inits = "random"
chains = 2
iter = 1000 
warmup = floor(iter / 2)
thin = 1
cores = 2
control = NULL
algorithm = "sampling"
future = FALSE
silent = TRUE
seed = 1
save_model = NULL
stan_model_args = list()
save_dso = TRUE
file = NULL 

  
  # validate arguments later passed to Stan
  dots <- list()
  testmode <- isTRUE(dots$testmode)
  dots$testmode <- NULL
  algorithm <- "sampling"
  silent <- brms:::as_one_logical(silent)
  iter <- brms:::as_one_numeric(iter)
  warmup <- brms:::as_one_numeric(warmup)
  thin <- brms:::as_one_numeric(thin)
  chains <- brms:::as_one_numeric(chains)
  cores <- brms:::as_one_numeric(cores)
  future <- brms:::as_one_logical(future) && chains > 0L
  seed <- brms:::as_one_numeric(seed, allow_na = TRUE)
  if (is.character(inits) && !inits %in% c("random", "0")) {
    inits <- get(inits, mode = "function", envir = parent.frame())
  }
  
  if (is.brmsfit(fit)) {
    # re-use existing model
    x <- fit
    icnames <- c("loo", "waic", "kfold", "R2", "marglik")
    x[icnames] <- list(NULL)
    sdata <- standata(x)
    x$fit <- rstan::get_stanmodel(x$fit)
  } else {  
    # build new model
    formula <- brms:::validate_formula(
      formula, data = data, family = family, 
      autocor = autocor, sparse = sparse
    )
    family <- brms:::get_element(formula, "family")
    bterms <- parse_bf(formula)
    data.name <- brms:::substitute_name(data)
    #data <- validate_data(data, bterms = bterms)
    #data2 <- validate_data2(
      #data2, bterms = bterms, 
      #brms:::get_data2_autocor(formula)
    #)
    prior <- brms:::check_prior(
      prior, formula = formula, data = data,
      sample_prior = sample_prior, warn = FALSE
    )
    # initialize S3 object
    x <- brms:::brmsfit(
      formula = formula, family = family, data = data, 
      data.name = data.name, prior = prior, 
      cov_ranef = cov_ranef, stanvars = stanvars, 
      stan_funs = stan_funs, algorithm = algorithm
    )
    x$ranef <- brms:::tidy_ranef(bterms, data = x$data)  
    x$exclude <- brms:::exclude_pars(
      x, save_ranef = save_ranef, 
      save_mevars = save_mevars,
      save_all_pars = save_all_pars
    )
    # x$model <- make_stancode(
    #   formula, data = data, prior = prior, 
    #   cov_ranef = cov_ranef, sample_prior = sample_prior, 
    #   knots = knots, stanvars = stanvars, stan_funs = stan_funs, 
    #   save_model = save_model
    # )
    x$model <- stancode2
    # generate Stan data before compiling the model to avoid
    # unnecessary compilations in case of invalid data
    # sdata <- make_standata(
    #   formula, data = data, prior = prior, data2 = data2,
    #   cov_ranef = cov_ranef, sample_prior = sample_prior,
    #   knots = knots, stanvars = stanvars
    # )
    sdata <- standata2
    stopifnot(is.list(stan_model_args))
    silence_stan_model <- !length(stan_model_args)
    stan_model_args$model_code <- x$model
    # if (!isTRUE(save_dso)) {
    #   warning2("'save_dso' is deprecated. Please use 'stan_model_args'.")
    #   stan_model_args$save_dso <- save_dso
    # }
    message("Compiling the C++ model")
    x$code <- brms:::eval_silent(
      do_call(rstan::stan_model, stan_model_args),
      silent = silence_stan_model, type = "message"
    )
  }
  
  args <- brms:::nlist(
    object = x$fit, data = sdata, pars = x$exclude, 
    include = FALSE, algorithm, iter, seed
  )
  args[names(dots)] <- dots
  message("Start sampling")
  #args$algorithm <- "sampling"
  if (args$algorithm == "sampling") {
    args$algorithm <- NULL
    c(args) <- brms:::nlist(
      init = inits, warmup, thin, control, 
      show_messages = !silent
    )
    # if (future) {
    #   if (cores > 1L) {
    #     brms:::warning2("Argument 'cores' is ignored when using 'future'.")
    #   }
    #   args$chains <- 1L
    #   futures <- fits <- vector("list", chains)
    #   for (i in seq_len(chains)) {
    #     args$chain_id <- i
    #     if (is.list(inits)) {
    #       args$init <- inits[i]
    #     }
    #     futures[[i]] <- future::future(
    #       brms::do_call(rstan::sampling, args), 
    #       packages = "rstan"
    #     )
    #   }
    #   for (i in seq_len(chains)) {
    #     fits[[i]] <- future::value(futures[[i]]) 
    #   }
    #   x$fit <- rstan::sflist2stanfit(fits)
    #   rm(futures, fits)
    # } else {
    args <- brms:::nlist(chains, cores)
    args <- brms:::nlist(
      object = x$fit, data = sdata, pars = x$exclude, 
      include = FALSE, algorithm, iter, seed
    )
    # need to recompile to overwrite x$fit 
    args$algorithm <- "NUTS"
      x$fit <- rstan::sampling(object = x$code, data = sdata, iter = 1000, seed = 1, chains = 2)
      
      me <- marginal_smooths(x, nsamples = 200, spaghetti = TRUE, nug = 0.1)
      plot(me)
    

  if (!testmode) {
    x <- rename_pars(x)
  }
  if (!is.null(file)) {
    write_brmsfit(x, file)
  }
  x
                }
