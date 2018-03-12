#--- stochastic volatility modeling ----------------------------------------

# model:
# dXt = (alpha - (tau*Vt)/2)dt + (tau*Vt)^(1/2) dBt
# dVt = -gamma(Vt - mu)dt + sigma Vt^lambda dBt

# Yt = (Xt, Vt)
# theta = (alpha, gamma, mu, sigma, lambda, tau)

require(numDeriv) # a package for numerical derivatives
source("proflm-functions.R") # for profile likelihood

#' Simulate data from SV model.
#'
#' @param n Length of SV time series.
#' @param Y0 Initial value of SV model, i.e., prior to first observation time.
#' @param dt Time between observations.
#' @param theta Vector of SV model parameters, namely \code{theta = c(alpha, gamma, mu, sigma, lambda, tau)}.
#' @return A \code{n x 2} matrix consisting of discretely observed series \code{(Xt, Vt)}.
sv.sim <- function(n, Y0, dt, theta) {
  # extract parameters
  alpha <- theta[1]
  gamma <- theta[2]
  mu <- theta[3]
  sigma <- theta[4]
  lambda <- theta[5]
  tau <- theta[6]
  # allocate space
  Yt <- matrix(NA, n, 2)
  colnames(Yt) <- c("X", "V")
  Yt[1,] <- Y0
  for(ii in 2:n) {
    Xt <- Yt[ii-1,1] # today's values
    Vt <- Yt[ii-1,2]
    # tomorrow's means
    muX <- Xt + (alpha - tau*Vt/2)*dt
    muV <- Vt + -gamma*(Vt - mu)*dt
    # tomorrow's standard deviations
    sigX <- sqrt(tau*Vt * dt)
    sigV <- sigma * Vt^lambda * sqrt(dt)
    # generate tomorrow with normals
    Yt[ii,1] <- muX + sigX * rnorm(1)
    Yt[ii,2] <- muV + sigV * rnorm(1)
  }
  Yt
}

#' Loglikelihood function for SV model.
#'
#' @param theta Vector of parameter values (see \code{\link{sv.sim}}).
#' @param Yt A \code{n x 2} matrix of discrete time series observations \code{Xt, Vt}.
#' @param dt Time between observations.
#' @return Scalar; loglikelihood evaluated at \code{theta}.
sv.loglik <- function(theta, Yt, dt) {
  # extract parameters
  alpha <- theta[1]
  gamma <- theta[2]
  mu <- theta[3]
  sigma <- theta[4]
  lambda <- theta[5]
  tau <- theta[6]
  n <- nrow(Yt) # number of observations
  Xt <- Yt[-n,1]
  Vt <- Yt[-n,2]
  dXt <- diff(Yt[,1])
  dVt <- diff(Yt[,2])
  # means
  muX <- (alpha - tau*Vt/2)*dt
  muV <- -gamma*(Vt - mu)*dt
  # standard deviations
  sigX <- sqrt(tau*Vt * dt)
  sigV <- sigma * Vt^lambda * sqrt(dt)
  # loglikelihood
  ll <- sum(dnorm(x = dXt, mean = muX, sd = sigX, log = TRUE)) # X
  ll + sum(dnorm(x = dVt, mean = muV, sd = sigV, log = TRUE)) # V
}

#--- profile likelihood --------------------------------------------------------

#' Conditional MLE for SV model.
#'
#' Analytic calculation of \code{phi = (alpha, gamma, mu, sigma)} for given value of \code{eta = (lambda, tau)}.
#'
#' @param eta Value of \code{lambda, tau} on which to condition.
#' @param Yt Matrix of observed data (see \code{\link{sv.loglik}}).
#' @param dt Time between observations.
#' @param eta.out Logical; if \code{FALSE} return conditional MLE \code{phi.hat} only, otherwise return the vector \code{(eta, phi.hat)}.
#' @return The conditional MLE for given \code{eta}.
sv.profmle <- function(eta, Yt, dt, eta.out = TRUE) {
  lambda <- eta[1]
  tau <- eta[2]
  dXt <- diff(Yt[,1])
  dVt <- diff(Yt[,2])
  n <- nrow(Yt)-1
  Vt <- Yt[1:n,2]
  # alpha
  tVt <- tau*Vt
  y <- dXt + tVt/2 * dt
  X <- matrix(dt, nrow = n, ncol = 1)
  diagV <- tVt * dt # diagonal of variance matrix
  suff <- lm.suff(y = y, X = X, V = diagV)
  alpha <- suff$beta[1]
  # gamma, mu, sigma
  y <- dVt
  X <- cbind(Vt*dt, dt)
  diagV <- Vt^(2*lambda) * dt
  suff <- lm.suff(y = y, X = X, V = diagV) # diagonal of variance matrix
  gamma <- -suff$beta[1]
  mu <- suff$beta[2]/gamma
  sigma <- suff$sigma
  mle <- c(alpha, gamma, mu, sigma)
  names(mle) <- c("alpha", "gamma", "mu", "sigma")
  if(eta.out) mle <- c(mle,
                       lambda = as.numeric(lambda), tau = as.numeric(tau))
  mle
}

#' Profile loglikelihood function for SV model.
#'
#' The value of \code{eta = (lambda, tau)} maximizing this function is the MLE \code{eta.hat} of the full likelihood.
#'
#' @param eta Value of \code{lambda, tau}.
#' @param Yt Matrix of observed data (see \code{\link{sv.loglik}}).
#' @param dt Time between observations.
#' @return Scalar; value of the profile likelihood.
#' @details Equivalent to (but more efficient than)
#'
#' \code{sv.loglik(theta = sv.profmle(eta, Yt, dt, eta.out = TRUE), Yt, dt)}.
sv.profll <- function(eta, Yt, dt) {
  lambda <- eta[1]
  tau <- eta[2]
  dXt <- diff(Yt[,1])
  dVt <- diff(Yt[,2])
  n <- nrow(Yt)-1
  Vt <- Yt[1:n,2]
  # X component
  tVt <- tau*Vt
  y <- dXt + tVt/2 * dt
  X <- matrix(dt, nrow = n, ncol = 1)
  diagV <- tVt * dt # diagonal of variance matrix
  suff <- lm.suff(y = y, X = X, V = diagV)
  llX <- lm.profll(suff = suff, sigma = 1)
  # V component
  y <- dVt
  X <- cbind(Vt*dt, dt)
  diagV <- Vt^(2*lambda) * dt # diagonal of variance matrix
  suff <- lm.suff(y = y, X = X, V = diagV)
  llV <- lm.profll(suff = suff)
  llX + llV
}

#' Residuals of SV model.
#'
#' @param Yt Matrix of observed data (see \code{\link{sv.loglik}}).
#' @param dt Time between observations.
#' @param theta Vector of parameter values (see \code{\link{sv.sim}}).
#' @return A matrix of size \code{(nrow(Yt)-1) x 2} of model residuals.
#' @details In the context of the Euler approximation, \code{dXt} and \code{dVt} are independent normals given the past history of the time series.  So if each is subtracted its conditional mean and divided by its standard deviation, the reult is a 2-column matrix of iid \code{N(0,1)}.  This requires \code{theta} to be known.  Substituting the true value of \code{theta} by the MLE leads to approximately normal residuals, which can be used for goodness-of-fit testing.
sv.resid <- function(Yt, dt, theta) {
  # extract parameters
  alpha <- theta[1]
  gamma <- theta[2]
  mu <- theta[3]
  sigma <- theta[4]
  lambda <- theta[5]
  tau <- theta[6]
  n <- nrow(Yt) # number of observations
  Xt <- Yt[-n,1]
  Vt <- Yt[-n,2]
  dXt <- diff(Yt[,1])
  dVt <- diff(Yt[,2])
  # means
  muX <- (alpha - tau*Vt/2)*dt
  muV <- -gamma*(Vt - mu)*dt
  # standard deviations
  sigX <- sqrt(tau*Vt * dt)
  sigV <- sigma * Vt^lambda * sqrt(dt)
  cbind(X = (dXt - muX)/sigX, V = (dVt - muV)/sigV) # residuals
}

#-------------------------------------------------------------------------------

#' Maximum likelihood estimation for the SV model.
#'
#' @param Yt Matrix of observed data (see \code{\link{sv.loglik}}).
#' @param dt Time between observations.
#' @param eta0 Optional initial value for \code{eta = (lambda, tau)}
#' @param var.calc Logical; whether or not to estimate variance of the MLE.
#' @param ... Optional control arguments to \code{optim}.
#' @return The MLE of \code{theta} and optionally an estimate of its variance matrix.
#' @details Uses profile likelihood to reduce the 6D optimization problem to 2D.  Optional initial estimates of \code{eta = (lambda, tau)} are provided to minimize user input to \code{optim}.
sv.fit <- function(Yt, dt, eta0, var.calc = TRUE, ...) {
  # param names
  theta.names <- c("alpha", "gamma", "mu", "sigma", "lambda", "tau")
  if(missing(eta0)) {
    # lambda0 = 1/2 (historical significance)
    # tau0: fit gbm with sigma^2 = tau mean(Vt)
    eta0 <- c(.5, var(diff(Yt[,1]))/dt/mean(Yt[,2]))
  }
  # -ve profile likelihood
  nlp <- function(eta) {
    if(any(eta < 0)) return(Inf)
    -sv.profll(eta = eta, Yt = Yt, dt = dt)
  }
  # optimization
  eta.fit <- optim(par = eta0,
                   fn = nlp, control = list(...))
  if(eta.fit$convergence != 0) stop("optim did not converge.")
  # full mle
  theta.mle <- sv.profmle(eta = eta.fit$par, Yt = Yt, dt = dT)
  names(theta.mle) <- theta.names
  if(var.calc) {
    # variance calculation
    Ve <- hessian(func = function(theta) {
      # parameter restrictions: gamma, mu, sigma, lambda, tau > 0
      if(any(theta[-1] <= 0)) return(Inf)
      -sv.loglik(theta = theta, Yt = Yt, dt = dT)
    }, x = theta.mle)
    Ve <- solve(Ve)
    colnames(Ve) <- theta.names
    rownames(Ve) <- theta.names
    out <- list(mle = theta.mle, var = Ve)
  } else out <- theta.mle
  out
}
