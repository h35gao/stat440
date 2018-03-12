#--- diagnostic functions -------------------------------------------------

# trace plot
mcmc.trace <- function(x, thin = 1, theta.name, add = FALSE, ...) {
  n <- length(x)
  iter <- seq(1, n, by = thin)
  if(missing(theta.name)) theta.name <- ""
  if(!add) {
    plot(x = iter, y = x[iter], type = "l",
         xlab = "Iteration", ylab = theta.name, main = theta.name, ...)
  } else {
    lines(x = iter, y = x[iter], ...)
  }
}

# acf plot
mcmc.acf <- function(x, lag.max = NULL, theta.name, add = FALSE, ...) {
  if(missing(theta.name)) theta.name <- ""
  aa <- acf(x = x, lag.max = lag.max, plot = FALSE)
  lag <- drop(aa$lag)
  aa <- drop(aa$acf)
  if(!add) {
    plot(x = lag, y = aa,
         main = theta.name, ylab = "Autocorrelation", xlab = "Lag", ...)
  } else {
    points(x = lag, y = aa, ...)
  }
}

# density plot
mcmc.density <- function(x, type.dens = c("density", "histogram"),
                         theta.name, add = FALSE, ...) {
  if(missing(theta.name)) theta.name <- ""
  type.dens <- match.arg(type.dens)
  if(type.dens == "density") {
    X <- density(x)
    if(!add) {
      plot(X, xlab = theta.name, main = theta.name, ...)
    } else lines(X, ...)
  } else if(type.dens == "histogram") {
    hist(x, xlab = theta.name, main = theta.name, ...)
  }
}

# effective sample size
effect.size <- function(x) {
  aa <- acf(x, plot = FALSE)
  length(x)/(2*sum(aa$acf)-1)
}

#--- generic Metropolized IID sampler --------------------------------------

#' @param X.prop matrix of proposals
#' @param ld.prop log-density of proposals in proposal distribution
#' @param ld.targ log-density of proposals in target distribution.  Invalid draws should be flagged with \code{-Inf}.
#' @param burn Number of samples from output to discard.  In other words, output as dimension \code{(nrow(X.prop) - burn) x ncol(X.prop)}
#' @param acc.out Logical, whether or not to return acceptance rate.
#' @return MCMC output with possibly also the acceptance rate.
miid.sampler <- function(X.prop, ld.prop, ld.targ,
                         burn, acc.out = FALSE, ind.out = FALSE) {
  nsamples <- nrow(X.prop)
  if(missing(burn)) burn <- min(floor(nsamples)/10, 1e3)
  X <- matrix(NA, nsamples, ncol(X.prop))
  colnames(X) <- colnames(X.prop)
  II <- rep(NA, nsamples) # indices of output rows of X.prop
  paccept <- 0 # acceptance probability
  # initialize
  X.old <- X.prop[1,]
  X[1,] <- X.old
  lr <- ld.targ - ld.prop
  lr.old <- lr[1]
  II.old <- 1
  II[1] <- 1
  for(ii in 2:nsamples) {
    X.new <- X.prop[ii,]
    if(ld.targ[ii] > -Inf) {
      lr.new <- lr[ii]
      lacc <- lr.new-lr.old # log acceptance rate
      if(lacc > 0 || runif(1) < exp(lacc)) {
        X.old <- X.new
        lr.old <- lr.new
        II.old <- ii
        paccept <- paccept+1
      }
    }
    # storage
    X[ii,] <- X.old
    II[ii] <- II.old
  }
  # delete burn-in
  ans <- X[(burn+1):nsamples,]
  paccept <- paccept/nsamples # acceptance rate
  message("MH Acceptance Rate: ", round(paccept*100), "%")
  if(acc.out || ind.out) {
    ans <- list(X = ans)
    if(acc.out) ans <- c(ans, list(accept = paccept))
    if(ind.out) ans <- c(ans, list(ind = II[(burn+1):nsamples]))
  }
  ans
}
