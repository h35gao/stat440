#--- some basic SDE functions --------------------------------------------------

#' Simulate Brownian motion.
#'
#' @param tseq Sequence of timepoints at which to simulate the process.  Must be positive and sorted in increasing order.
#' @param B0 Value of the process at time \code{t = 0}.
#' @return A single random realization of the process at the desired timepoints.
rbm <- function(tseq, B0 = 0) {
  if(!identical(tseq, sort(tseq))) {
    stop("tseq must be sorted in increasing order.")
  }
  delt <- diff(c(0, tseq))
  N <- length(delt)
  Z <- rnorm(N, sd = sqrt(delt))
  cumsum(Z)
}

#' Simulate a Brownian bridge.
#'
#' Generates Brownian motion conditioned on its value at previous and subsequent timepoints.
#'
#' @param tseq Vector of time points at which to generate the process.
#' @param tL Vector of the same length as \code{tseq} giving the time points prior to \code{tseq} at which the proces is observed (see Details).
#' @param Vector of time points after \code{tseq} at which the process is observed (see Details).
#' @param BL Value of the process at times \code{tL}.
#' @param BU Value of the process at times \code{tU}.
#' @return Vector \code{Bt} of random draws from the distribution \code{p(Bt | BL, BU)}.
rbb <- function(tseq, tL, tU, BL, BU) {
  if(any(tseq > tU) || any(tseq < tL)) {
    stop("tseq must be between tL and tU.")
  }
  t1 <- tseq - tL
  t2 <- tU - tseq
  sig2 <- 1/(1/t1 + 1/t2)
  mu <- (BL/t1 + BU/t2) * sig2
  rnorm(length(tseq), mean = mu, sd = sqrt(sig2))
}

#-------------------------------------------------------------------------------

#' Exact transition density of geometric Brownian motion.
#'
#' Calculates the conditional distribution \code{p(St | S0)}, where St is gBM.
#'
#' @param s Value(s) at which to evaluate \code{St}.
#' @param alpha Rate parameter of gBM.
#' @param sigma Volatility parameter of gBM.
#' @param t Time between \code{S0} and \code{St}.
#' @param s0 Initial value of process on which to condition.
#' @param log Logical; whether or not conditional PDF is returned on the log-scale.
#' @return Conditional PDF evaluation(s).
#' @details Geometric Brownian motion is a stochastic process \eqn{S_t} satisfying the stochastic differential equation (SDE)
#' \deqn{
#' d S_t = \alpha S_t dt + \sigma S_t dB_t.
#' }
dgbm <- function(s, alpha, sigma, s0, t, log = FALSE) {
  dlnorm(s, log = log,
         meanlog = (alpha - sigma^2/2)*t + log(S0),
         sdlog = sigma * sqrt(t))
}

#' Euler simulation of geometric Brownian motion.
#'
#' Given an initial value \code{S0} at time \code{t = 0}, generates iid random vectors \code{S = (S1, ..., SN)}, each of which is approximately distributed as the time series of a geometric Brownian motion at times \code{t = (dT, 2*dT, ..., N*dT)}.
#'
#' @param n Number of random vectors (i.e., time series) to generate.
#' @param alpha Rate parameter of gBm.
#' @param sigma Volatility parameter of gBm.
#' @param dT Time between observations.
#' @param N Length of each time series.
#' @return An \code{N x n} matrix of iid time series, each of which is a column.
gbm.sim <- function(n, alpha, sigma, dT, s0, N) {
  SS <- matrix(NA, n, N+1) # always faster to fill matrix by column than by row
  SS[,1] <- s0 # initial value
  for(ii in 2:(N+1)) {
    dB <- rnorm(n, sd = sqrt(dT))
    SS[,ii] <- SS[,ii-1] + alpha*SS[,ii-1] * dT + sigma*SS[,ii-1] * dB
  }
  t(SS[,-1]) # customarily each column is a time series
}
