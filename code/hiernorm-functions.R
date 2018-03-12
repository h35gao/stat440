#--- functions for the 8 schools analysis ----------------------------------

# model is:
#
# xi ~ind N(mui, sig2i),   sig2i = sigmai^2
# mui ~iid N(lambda, tau^2)
#
# model parameters are mu and theta = (lambda, tau)
# mu are the so-called "random effects".
# theta are called "hyperparameters".

#--- frequentist inference: profile likelihood ---------------------------------

#' Loglikelihood of hyperparameters \code{theta = (lambda, tau)}.
#'
#' @param lambda Scalar; random-effects mean.
#' @param tau Scalar; random-effects standard deviation.
#' @param x Vector of school sample means.
#' @param sig2 Vector of (true) variances for each sample mean.
#' @return Loglikelihood evaluated at \code{theta}.
loglik.theta <- function(lambda, tau, x, sig2) {
  sum(dnorm(x, mean = lambda, sd = sqrt(sig2 + tau^2), log = TRUE))
}

#' Conditional MLE of \code{lambda} given \code{tau}.
#'
#' @param tau Scalar; random-effects standard deviation.
#' @param x Vector of school sample means.
#' @param sig2 Vector of (true) variances for each sample mean.
#' @return Scalar; conditional MLE of \code{lambda} given \code{tau}.
lambda.prof <- function(tau, x, sig2) {
  wtau <- 1/(sig2 + tau^2)
  sum(wtau*x)/sum(wtau)
}

#' Profile loglikelihood of \code{tau}.
#'
#' @param tau Scalar; random-effects standard deviation.
#' @param x Vector of school sample means.
#' @param sig2 Vector of (true) variances for each sample mean.
#' @return Scalar; value of the profile loglikelihood.
lprof.tau <- function(tau, x, sig2) {
  V <- sig2 + tau^2
  wtau <- 1/V
  ltau <- sum(wtau*x)/sum(wtau)
  sum(dnorm(x, mean = ltau, sd = sqrt(V), log = TRUE))
}

#' Calculate MLE of \code{theta = (lambda, tau)}.
#'
#' @param x Vector of school sample means.
#' @param sig2 Vector of (true) variances for each sample mean.
#' @param interval Interval over which to maximize the profile likelihood for \code{tau}.
#' @return The MLE of \code{theta}.
theta.mle <- function(x, sig2, interval = c(0, 50)) {
  # maximize profile likelihood
  tau.mle <- optimize(f = lprof.tau, x = x, sig2 = sig2,
                      interval = interval, maximum = TRUE)$max
  # get profiled lambda
  theta.mle <- c(lambda.prof(tau = tau.mle, x = x, sig2 = sig2), tau.mle)
  names(theta.mle) <- c("lambda", "tau")
  theta.mle
}

#--- bayesian inference --------------------------------------------------------

#' Conditional distribution of random-effects mean.
#'
#' @param lambda Scalar; random-effects mean.
#' @param tau Scalar; random-effects standard deviation.
#' @param x Vector of school sample means.
#' @param sig2 Vector of (true) variances for each sample mean.
#' @return Log of the conditional distribution \code{p(lambda | tau, x, sig2)}.
#' @details The prior distribution is of the form \code{pi(lambda, tau) = g(tau)}, i.e., does not depend on \code{lambda}.
lpcond.lambda <- function(lambda, tau, x, sig2) {
  wtau <- 1/(sig2 + tau^2)
  s2tau <- 1/sum(wtau)
  ltau <- sum(wtau*x)*s2tau
  dnorm(lambda, ltau, sqrt(s2tau), log = TRUE)
}


#' Marginal distribution of random-effects standard deviation.
#'
#' @param tau Scalar; random-effects standard deviation.
#' @param x Vector of school sample means.
#' @param sig2 Vector of (true) variances for each sample mean.
#' @return Log of the marginal distribution of \code{tau}, *excluding* the prior (See Details).
#' @details The prior distribution is of the form \code{pi(lambda, tau) = g(tau)}, i.e., does not depend on \code{lambda}.  Thus, \code{lpmarg.tau} returns the log of \code{p(tau | x, sig2) / g(tau)}.
lpmarg.tau <- function(tau, x, sig2) {
  V <- sig2 + tau^2
  wtau <- 1/V
  s2tau <- 1/sum(wtau)
  ltau <- sum(wtau*x)*s2tau
  .5 * log(s2tau) + sum(dnorm(x, ltau, sqrt(V), log = TRUE))
}

#' Sample from the conditional distribution \code{p(lambda | tau, x, sig2)}.
#'
#' @param n Number of random draws to return.
#' @param tau Random-effects standard deviation.  Can be a scalar, or a vector of length \code{n}.
#' @param x Vector of school sample means.
#' @param sig2 Vector of (true) variances for each sample mean.
#' @return Vector of random draws.
#' @details If \code{tau} is scalar, draws are iid.  If \code{tau} is a vector, then \code{lambda}'s are independently drawn from the corresponding distribution, i.e., \code{lambda_i ~ind p(lambda | tau_i, x, sig2)}.
rcond.lambda <- function(n, tau, x, sig2) {
  # since we'll have to perform calculations for each sample n,
  # more efficient to do colSums than rowSums,
  # since elements of former are stored contiguously in memory.
  wtau <- matrix(tau, nrow = length(x), ncol = n, byrow = TRUE)
  wtau <- 1/(sig2 + wtau^2)
  s2tau <- 1/colSums(wtau) # sigma^2(tau)
  ltau <- colSums(wtau*x)*s2tau # lambda_hat(tau) = sum(w * x)
  rnorm(n, ltau, sqrt(s2tau))
}

#' Sample from the conditional distribution \code{p(mu | lambda, tau, x, sig2)}.
#'
#' @param n Number of random draws to return.
#' @param lambda Random-effects mean.  Can be a scalar, or a vector of length \code{n}.
#' @param tau Random-effects standard deviation.  Can be a scalar, or a vector of length \code{n}.
#' @param x Vector of school sample means.
#' @param sig2 Vector of (true) variances for each sample mean.
#' @return Matrix of in which each column is a random draw of \code{mu}.
rcond.mu <- function(n, lambda, tau, x, sig2) {
  K <- length(x)
  Tau <- matrix(tau, nrow = K,
                ncol = n, byrow = TRUE)
  Lambda <- matrix(lambda, nrow = K, ncol = n)
  Btau <- sig2/(sig2 + Tau^2)
  mu <- rnorm(n*K, mean = Btau*Lambda + (1-Btau)*x, sqrt((1-Btau)*sig2))
  matrix(mu, K, n)
}

# plot multiple density functions on the same plot
# everything else gets added later
multi.dens <- function(X, col, xlim, ylim, ...) {
  nvar <- ncol(X)
  if(missing(col)) col <- rainbow(nvar+1)[1:nvar]
  suppressWarnings(dlist <- apply(X, 2, density, ...))
  # empty plot to get axes
  if(missing(xlim)) xlim <- range(sapply(dlist, function(dd) dd$x))
  if(missing(ylim)) ylim <- range(sapply(dlist, function(dd) dd$y))
  plot(x = 0, type = "n", xlim = xlim, ylim = ylim,
       xlab = "", ylab = "")
  # now add lines
  for(ii in 1:nvar) {
    lines(dlist[[ii]], col = col[ii], ...)
  }
  invisible(NULL) # return invisibly
}
