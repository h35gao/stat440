#--- r/d/p/q methods for the Generalized Gamma Distribution ---------------

# Model is:
# Y ~ Gen-Gamma(alpha, beta, lambda)
#  <=> Y = X^(1/lambda)/beta,   X ~ Gamma(alpha/lambda, 1).

#' Density evaluation.
#'
#' @param x Vector of quantiles at which to evaluate the density.
#' @param alpha Vector of shape parameters.
#' @param beta Vector of rate parameters.
#' @param lambda Vector of power parameters.
#' @param log Logical; whether or not to return the density on the log-scale.
#' @return Vector of density evaluations.
dggamma <- function(x, alpha, beta, lambda, log = FALSE) {
  lbe <- log(x*beta)
  ld <- dgamma(exp(lambda * lbe), shape = alpha/lambda, log = TRUE)
  ld <- ld + log(lambda*beta) + (lambda-1) * lbe
  if(!log) ld <- exp(ld)
  ld
}

#' Random number generation.
#'
#' @param n Number of draws to generate.
#' @param alpha Vector of shape parameters.
#' @param beta Vector of rate parameters.
#' @param lambda Vector of power parameters.
#' @return Vector of random draws.
rggamma <- function(n, alpha, beta, lambda) {
  rgamma(n, shape = alpha/lambda)^(1/lambda)/beta
}

#' Cumulative density.
#'
#' @param q Vector of quantiles.
#' @param alpha Vector of shape parameters.
#' @param beta Vector of rate parameters.
#' @param lambda Vector of power parameters.
#' @param lower.tail Logical, whether to return \code{P(X < q)} or \code{P(X > q)}.  The latter uses a specialized algorithm which is much more accurate than \code{1-P(X < q)} for small values of \code{q}.
#' @param log.p Logical, whether to return the CDF or the its logarithm, whether the latter uses an algorithm which is much more accurate than simply calculating \code{log(CDF)}.
#' @return A vector of CDF evaluations.
pggamma <- function(q, alpha, beta, lambda,
                    lower.tail = TRUE, log.p = FALSE) {
  pgamma((q*beta)^lambda, shape = alpha/lambda, rate = 1,
         lower.tail = lower.tail, log.p = log.p)
}

#' Quantile function.
#'
#' @param p Vector of probalities.
#' @param alpha Vector of shape parameters.
#' @param beta Vector of rate parameters.
#' @param lambda Vector of power parameters.
#' @param lower.tail Logical, whether to return quantile at \code{p} or \code{1-p}.
#' @param log.p Logical, whether to return quantile at \code{p} or \code{log(p)}.
#' @return A vector of CDF evaluations.
qggamma <- function(p, alpha, beta, lambda,
                    lower.tail = TRUE, log.p = FALSE) {
  qgamma(p, shape = alpha/lambda, rate = 1,
         lower.tail = lower.tail, log.p = log.p)^(1/lambda)/beta
}
