#--- functions for the Weibull distribution ------------------------------------

require(mvtnorm)
require(numDeriv)
source("mcmc-functions.R")

# Model is:

# Y ~ Weibull(gamma, lambda)
#   <=>  Y = lambda * X^gamma,   where   X ~ Expo(1).

#-------------------------------------------------------------------------------

#' Loglikelihood function for Weibull distribution.
#'
#' @param gamma Shape parameter (scalar).
#' @param lambda Scale parameter (scalar).
#' @param logy Vector of log-observations (see Details).
#' @return Loglikelihood a parameter values.
#' @details If \eqn{X \sim Expo(1)}, then \eqn{Y = \lambda X^{1/\gamma} \sim Weibull(\gamma, \lambda)}.
#'
#' Since the likelihood function requires calculation of \code{y^gamma = exp(gamma * log(y))} for the vector of observations \code{y}, and \code{exp/log} are fairly expensive operations, for multiple evaluations computation time is roughly halved by inputting the data on the log-scale.
weibull.loglik <- function(gamma, lambda, logy) {
  n <- length(logy)
  leta <- gamma * log(lambda)
  gamma * sum(logy) - sum(exp(gamma*logy - leta)) + n * (log(gamma) - leta)
}

#-------------------------------------------------------------------------------

#' Random Walk Metropolis sampling for the Weibull parameter posterior distribution.
#'
#' @param nsamples Number of posterior iterations (not including burn-in).
#' @param y Vector of survival times (non-negative).
#' @param theta0 Vector of length 2 giving initial values for \code{theta = (gamma, lambda)}.
#' @param rwsd Vector of length 2 giving the standard deviation of the Metropolis proposal for each component of \code{theta = (gamma, lambda)}.
#' @param burn Integer number of initial iterations to discard as burn-in.  The default value is \code{min(nsamples/10, 1000)}.
#' @param acc.out Logical; whether or not to return the proportion of accepted proposals.
#' @return A \code{nsamples x 2} matrix of posterior draws, or a list with elements \code{Theta} and \code{accept}, the first of which are the posterior draws, the second being the proportion of accepted proposals.
#' @details The prior distribution is flat on \code{theta = (gamma, lambda)}, i.e., \code{pi(theta) ~ 1}.
weibull.rwm <- function(nsamples, y, theta0, rwsd,
                        burn, acc.out = FALSE) {
  logy <- log(y) # preprocess data
  # default burn-in
  if(missing(burn)) burn <- min(floor(nsamples/10), 1e3)
  Theta <- matrix(NA, nsamples, 2) # storage
  colnames(Theta) <- c("gamma", "lambda")
  paccept <- 0 # keep MH acceptance probability
  # initialize
  theta.curr <- theta0
  # only evaluate the posterior once per step
  # for flat prior, logpost = loglik
  lp.curr <- weibull.loglik(theta.curr[1], theta.curr[2], logy)
  for(ii in (-burn+1):nsamples) {
    theta.prop <- rnorm(2, mean = theta.curr, sd = rwsd)
    if(all(theta.prop > 0)) {
      # only bother updating if draw is "valid", otherwise p(theta | y) = 0
      lp.prop <- weibull.loglik(theta.prop[1], theta.prop[2], logy)
      lacc <- lp.prop-lp.curr # log acceptance rate
      if(lacc > 0 || runif(1) < exp(lacc)) {
        # ||: only evaluate 2nd if 1st is true
        # automatic accept if acc = exp(lacc) > 1
        theta.curr <- theta.prop
        lp.curr <- lp.prop
        paccept <- paccept+1 # count number of accepted proposals
      }
    }
    # storage
    if(ii > 0) Theta[ii,] <- theta.curr
  }
  paccept <- paccept/(nsamples+burn) # acceptance rate
  message("MH Acceptance Rate: ", round(paccept*100), "%")
  if(acc.out) {
    ans <- list(Theta = Theta, accept = paccept)
  } else ans <- Theta
  ans
}

#-------------------------------------------------------------------------------

#' Calculate the mode and quadrature of the Weibull loglikelihood.
#'
#' @param y Vector of survival times (non-negative).
#' @param theta0 Vector of length 2 giving initial values for \code{theta = (gamma, lambda)}.
#' @param ... Further arguments to pass to \code{optim} (see Details).
#' @return A list with elements:
#' \describe{
#'   \item{\code{mode}}{The mode of the Weibull loglikelihood: \code{mode = argmax_theta loglik(theta | y)}.}
#'   \item{\code{quad}}{The *negative* quadrature at the mode: \code{quad = -d^2/dtheta^2 loglik(theta = mode | y)}.}
#' }
#' @details Since \code{optim} performs minimization by default, its objective function here is the *negative* loglikelihood, so there's no need to set \code{fnscale = -1}.
weibull.mq <- function(y, theta0, ...) {
  logy <- log(y)
  # so can pass arbitrary controls to optim
  negloglik <- function(theta) {
    if(any(theta < 0)) return(Inf)
    -weibull.loglik(theta[1], theta[2], logy)
  }
  prop.mode <- optim(par = theta0, fn = negloglik, ...)
  if(prop.mode$convergence != 0) warning("optim did not converge.")
  prop.mode <- prop.mode$par
  names(prop.mode) <- c("gamma", "lambda")
  prop.quad <- hessian(func = negloglik, x = prop.mode)
  colnames(prop.quad) <- names(prop.mode)
  rownames(prop.quad) <- names(prop.mode)
  list(mode = prop.mode, quad = prop.quad)
}

#' Metropolized IID sampling for the Weibull posterior parameter distribution.
#'
#' @param nsamples Number of posterior iterations (not including burn-in).
#' @param y Vector of survival times (non-negative).
#' @param theta0 Vector of length 2 giving initial values optim to find the posterior mode of \code{theta = (gamma, lambda)}.
#' @param burn Integer number of initial iterations to discard as burn-in.  The default value is \code{min(nsamples/10, 1000)}.
#' @param acc.out Logical; whether or not to return the proportion of accepted proposals.
#' @return A \code{nsamples x 2} matrix of posterior draws, or a list with elements \code{Theta} and \code{accept}, the first of which are the posterior draws, the second being the proportion of accepted proposals.
#' @details The prior distribution is flat on \code{theta = (gamma, lambda)}, i.e., \code{pi(theta) ~ 1}.
weibull.miid <- function(nsamples, y, theta0,
                         burn, acc.out = FALSE) {
  prop.mq <- weibull.mq(y = y, theta0 = theta0)
  prop.mu <- prop.mq$mode
  prop.V <- solve(prop.mq$quad)
  if(missing(burn)) burn <- min(floor(nsamples/10), 1e3) # default burn-in
  # pre-draw all proposals
  Theta.prop <- rmvnorm(nsamples + burn,
                        mean = prop.mu, sigma = prop.V)
  # pre-compute all proposal densities
  lp.prop <- dmvnorm(Theta.prop, log = TRUE,
                     mean = prop.mu, sigma = prop.V)
  # pre-compute all target densities
  # can vectorize this slightly better, but not by much because of y^lambda
  logy <- log(y)
  logpost <- function(theta) {
    if(any(theta < 0)) return(-Inf)
    weibull.loglik(theta[1], theta[2], logy)
  }
  lp.targ <- apply(Theta.prop, 1, logpost)
  # initialize
  theta.curr <- theta0
  # log-acceptance rate denominator
  lacc.curr <- logpost(theta.curr) - dmvnorm(theta.curr, mean = prop.mu,
                                             sigma = prop.V, log = TRUE)
  Theta <- matrix(NA, nsamples, 2)
  colnames(Theta) <- c("gamma", "lambda")
  paccept <- 0 # keep MH acceptance probability
  # minimal for-loop
  for(ii in (-burn+1):nsamples) {
    II <- burn+ii
    if(lp.targ[II] != -Inf) {
      # only bother if proposal is valid
      theta.prop <- Theta.prop[II,]
      lacc.prop <- lp.targ[II] - lp.prop[II] # log-acceptance rate numerator
      lacc <- lacc.prop - lacc.curr # log acceptance rate
      if(lacc > 0 || runif(1) < exp(lacc)) {
        # automatic accept if acc = exp(lacc) > 1
        theta.curr <- theta.prop
        lacc.curr <- lacc.prop
        paccept <- paccept+1
      }
    }
    # storage
    if(ii > 0) Theta[ii,] <- theta.curr
  }
  paccept <- paccept/(nsamples+burn) # acceptance rate
  message("MH Acceptance Rate: ", round(paccept*100), "%")
  if(acc.out) {
    ans <- list(Theta = Theta, accept = paccept)
  } else ans <- Theta
  ans
}

#-------------------------------------------------------------------------------

#' Metropolis-within-Gibbs sampling for the Weibull posterior parameter distribution.
#'
#' @param nsamples Number of posterior iterations (not including burn-in).
#' @param y Vector of survival times (non-negative).
#' @param theta0 Vector of length 2 giving initial values for \code{theta = (gamma, lambda)}.
#' @param rwsd Vector of length 2 giving the standard deviation of the Metropolis proposal for each component of \code{theta = (gamma, lambda)}.
#' @param burn Integer number of initial iterations to discard as burn-in.  The default value is \code{min(nsamples/10, 1000)}.
#' @param acc.out Logical; whether or not to return the proportion of accepted proposals.
#' @return A \code{nsamples x 2} matrix of posterior draws, or a list with elements \code{Theta} and \code{accept}, the first of which are the posterior draws, the second being the proportion of accepted proposals.
#' @details The prior distribution is flat on \code{theta = (gamma, lambda)}, i.e., \code{pi(theta) ~ 1}.
weibull.mwg <- function(nsamples, y, theta0, rwsd,
                        burn, acc.out = FALSE) {
  logy <- log(y)
  if(missing(burn)) burn <- min(floor(nsamples/10), 1e3) # burnin
  Theta <- matrix(NA, nsamples, 2)
  colnames(Theta) <- c("gamma", "lambda")
  ntheta <- ncol(Theta) # number of parameters
  # keep MH acceptance probability, one for each parameter
  paccept <- rep(0, ntheta)
  names(paccept) <- colnames(Theta)
  # initialize
  theta.curr <- theta0
  # only evaluate the posterior once per acceptance rate calculation
  lp.curr <- weibull.loglik(theta.curr[1], theta.curr[2], logy)
  for(ii in (-burn+1):nsamples) {
    for(jj in 1:ntheta) {
      # similar to rw update, except once for each parameter
      # note that there are ntheta posterior evaluations per cycle
      # also, nested for-loop is very slow in R.
      theta.prop <- theta.curr
      theta.prop[jj] <- theta.prop[jj] + rwsd[jj] * rnorm(1)
      if(all(theta.prop > 0)) {
        # only bother updating if draw is "valid"
        lp.prop <- weibull.loglik(theta.prop[1], theta.prop[2], logy)
        lacc <- lp.prop-lp.curr # log acceptance rate
        if(lacc > 0 || runif(1) < exp(lacc)) {
          # automatic accept if acc = exp(lacc) > 1
          theta.curr <- theta.prop
          lp.curr <- lp.prop
          paccept[jj] <- paccept[jj]+1
        }
      }
    }
    # storage
    if(ii > 0) Theta[ii,] <- theta.curr
  }
  paccept <- paccept/(nsamples+burn) # acceptance rate
  message("MH Acceptance Rate:")
  message("gamma: ", round(paccept[1]*100), "%")
  message("lambda: ", round(paccept[2]*100), "%")
  if(acc.out) {
    ans <- list(Theta = Theta, accept = paccept)
  } else ans <- Theta
  ans
}

#-------------------------------------------------------------------------------

#' Conditional posterior log-density of transformed Weibull parameter \code{eta = 1/lambda^gamma}.
#'
#' @param eta Vector of quantiles at which to evaluate the log-density.
#' @param gamma Value of the Weibull shape parameter on which to condition.
#' @param logy Vector of observations on the log scale.
#' @param alpha0,beta0 Shape and rate parameter of the conjugate prior gamma distribution for \code{eta} (see Details).
#' @param Tgamma Optional precomputed value of \code{sum(y^gamma)}.
#' @return Vector of log-density evaluations.
#' @details The conditional conjugate prior for \code{eta} is a Gamma distribution.  Note the change-of-variables formula between the prior on the \code{(eta,lambda)} and \code{(gamma,lambda)} scales:
#' \preformatted{
#' pi(gamma,lambda) = pi(eta,lambda) * ???.
#' }
#'
#' Also, note that vectorized evaluation is possible (i.e., vectors for each of \code{gamma, alpha0, beta0}) as long as \code{Tgamma} is supplied for each value of \code{gamma}.
weibull.lpcond <- function(eta, gamma, logy, alpha0, beta0, Tgamma) {
  if(missing(Tgamma)) {
    Tgamma <- sum(exp(gamma[1] * logy))
  }
  if(length(Tgamma) != length(gamma)) {
    stop("Vectorization in gamma only possible when Tgamma of same length is supplied.")
  }
  dgamma(x = eta, shape = alpha0 + length(logy),
         rate = beta0 + Tgamma, log = TRUE)
}

#' Conditional posterior draws of the transformed Weibull parameter \code{eta = 1/lambda^gamma}.
#'
#' @param n Number of posterior draws.
#' @param gamma Value of the Weibull shape parameter on which to condition.
#' @param logy Vector of observations on the log scale.
#' @param alpha0,beta0 Shape and rate parameter of the conjugate prior gamma distribution for \code{eta} (see Details).
#' @param Tgamma Optional precomputed value of \code{sum(y^gamma)}.
#' @return Vector of log-density evaluations.
#' @details See \code{\link{weibull.lpcond}}.
weibull.rcond <- function(n, gamma, logy, alpha0, beta0, Tgamma) {
  if(missing(Tgamma)) {
    Tgamma <- sum(exp(gamma[1] * logy))
  }
  if(length(Tgamma) != length(gamma)) {
    stop("Vectorization in gamma only possible when Tgamma of same length is supplied.")
  }
  rgamma(n, shape = alpha0 + length(logy), rate = beta0 + Tgamma)
}


#' Marginal posterior log-density of the Weibull parameter \code{gamma}.
#'
#' @param gamma Weibull shape parameter at which to evaluate the marginal log-posterior.
#' @param logy Vector of observations on the log scale.
#' @param alpha0,beta0 Shape and rate parameter of the conjugate prior gamma distribution for \code{eta} (see Details).
#' @param Tgamma Optional precomputed value of \code{sum(y^gamma)}.
#' @return Log-density evaluation.
#' @details Calculation can also be vectorized as long as \code{Tgamma} is precomputed: see \code{\link{weibull.lpcond}}.
weibull.lpmarg <- function(gamma, logy, alpha0, beta0, Tgamma) {
  if(missing(Tgamma)) {
    Tgamma <- sum(exp(gamma[1] * logy))
  }
  if(length(Tgamma) != length(gamma)) {
    stop("Vectorization in gamma only possible when Tgamma of same length is supplied.")
  }
  n <- length(logy)
  ahat <- alpha0 + n
  bhat <- beta0 + Tgamma
  lgamma(ahat) - ahat * log(bhat) + gamma * sum(logy) + n * log(gamma)
}

#' Collapsed Metropolized IID sampling for the Weibull posterior parameter distribution.
weibull.cmiid <- function(nsamples, y, lambda.rng, burn, acc.out = FALSE) {
  if(missing(burn)) burn <- min(floor(nsamples/10), 1e3) # default burn-in
  # marginal log-posterior
  lpmarg <- function(gamma, Tgamma) {
    weibull.lpmarg(gamma, logy, alpha0 = 1, beta0 = 0, Tgamma = Tgamma)
  }
  # mode-quadrature parameters
  prop.mean <- optimize(lpmarg,
                        lower = lambda.rng[1], upper = lambda.rng[2],
                        maximum = TRUE)$max
  prop.sd <- hessian(func = lpmarg, x = prop.mean)
  prop.sd <- drop(sqrt(-1/prop.sd))
  # proposals
  gamma.prop <- rnorm(nsamples+burn, mean = prop.mean, sd = prop.sd)
  # log-density of proposal distribution
  ld.prop <- dnorm(gamma.prop, mean = prop.mean, sd = prop.sd, log = TRUE)
  # log-density of target distribution (i.e., marginal posterior)
  # step 1: precompute Tgamma. will need these for later
  T.prop <- sapply(gamma.prop, function(gg) sum(exp(gg * logy)))
  ld.targ <- rep(NA, nsamples+burn)
  bad.draws <- gamma.prop < 0 # bad draws
  ld.targ[bad.draws] <- -Inf
  # don't use ifelse for this, as ifelse still evaluates bad.draws,
  # just doesn't return them.
  ld.targ[!bad.draws] <- lpmarg(gamma.prop[!bad.draws], T.prop[!bad.draws])
  # generic MIID on p(gamma | y)
  gpost <- miid.sampler(X.prop = as.matrix(gamma.prop),
                        ld.prop = ld.prop, ld.targ = ld.targ,
                        burn = burn, acc.out = TRUE, ind.out = TRUE)
  gamma.post <- drop(gpost$X) # the posterior draws themselves
  paccept <- gpost$accept # MIID acceptance rate
  # analytic conditional draws of p(eta | gamma, y)
  T.post <- T.prop[gpost$ind] # precomputations for all accepted proposals
  eta.post <- weibull.rcond(nsamples, gamma.post, logy,
                            alpha0 = 1, beta0 = 0, Tgamma = T.post)
  lambda.post <- eta.post^(-1/gamma.post) # transform back to original scale
  Theta <- cbind(gamma = gamma.post, lambda = lambda.post)
  ## message("MH Acceptance Rate: ", round(paccept*100), "%")
  if(acc.out) {
    ans <- list(Theta = Theta, accept = paccept)
  } else ans <- Theta
  ans
}


#-------------------------------------------------------------------------------

# statistics to plot analytic Weibull posterior
# input: sequence of gamma and lambda values, along with the 2D grid of log-posterior density evaluations.  And alpha, the confidence level
# output: same as input + marginal, mean, and confidence interval for each parameter.
weibull.plot.stat <- function(gseq, lseq, lpmat, alpha = .95) {
  # calculate margins
  gldens <- exp(lpmat-max(lpmat)) # unnormalized p(gamma, lambda | y)
  gdens <- rowSums(gldens) # p(gamma | y)
  dg <- gseq[2]-gseq[1]
  gdens <- gdens/sum(gdens)/dg # normalize
  gmean <- sum(gseq*gdens)*dg # posterior mean
  # credible interval
  clim <- c((1-alpha)/2, 1 - (1-alpha)/2)
  gL <- gseq[which.min(abs(cumsum(gdens)*dg - clim[1]))]
  gU <- gseq[which.min(abs(cumsum(gdens)*dg - clim[2]))]
  ldens <- colSums(gldens) # p(lambda | y)
  dl <- lseq[2]-lseq[1]
  ldens <- ldens/sum(ldens)/dl
  lmean <- sum(lseq*ldens)*dl
  lL <- lseq[which.min(abs(cumsum(ldens)*dl - clim[1]))]
  lU <- lseq[which.min(abs(cumsum(ldens)*dl - clim[2]))]
  list(gseq = gseq, gdens = gdens, gmean = gmean, gCI = c(gL, gU),
       lseq = lseq, ldens = ldens, lmean = lmean, lCI = c(lL, lU),
       gldens = gldens)
}
