#--- multivariate normal distribution ------------------------------------------

#' Random number generation.
#'
#' @param n Number of random draws
#' @param mu Mean vector
#' @param V Variance matrix
#' @return A matrix containing the normal draws where each row is an iid observation.
rmvn <- function(n, mu, V) {
  d <- length(mu) # number of dimensions
  C <- chol(V) # cholesky decomposition: C upper-triangular and C'C = V
  Z <- matrix(rnorm(n*d), n, d) # n x d matrix of iid N(0,1)
  X <- Z %*% C # each row of X is iid N(0, C'C = V)
  X <- t(t(X) + mu) # each row of X is N(mu, V)
  X
}

#' PDF of a single observation.
#'
#' @param x Observation vector
#' @param mu Mean vector
#' @param V Variance matrix
#' @param log Logical; whether or not log-pdf is desired
#' @return Scalar value of (log-)pdf.
dmvn <- function(x, mu, V, log = FALSE) {
  d <- length(mu) # dimension of problem
  xVx <- t(x - mu) %*% solve(V) %*% (x - mu) # (x-mu)'V^{-1}(x-mu)
  ldV <- log(det(V)) # log(|V|)
  ll <- -.5 * (xVx + ldV + d * log(2*pi))
  if(!log) ll <- exp(ll) # always calculate on log scale
  ll
}

#' Loglikelihood function.
#'
#' @param mu Mean vector
#' @param V Variance matrix
#' @param Xbar Sample mean
#' @param S Sample sum-of-squares
#' @param n Sample size
#' @return Scalar; loglikelihood evaluated at \code{mu} and \code{V}.
#' @details Data input is sufficient statistics \code{Xbar} and \code{S} to avoid recalculating these for each \code{mu} and \code{V}.
lmvn <- function(mu, V, Xbar, S, n) {
  VS <- solve(V, S) # V^{-1} S
  # (Xbar-mu)'V^{-1}(Xbar-mu)
  xVx <- crossprod(Xbar - mu, solve(V, Xbar - mu))
  ldV <- determinant(x = V, log = TRUE)$mod[1] # log(|V|)
  -.5 * (sum(diag(VS)) + n * xVx + n * ldV)
}

#' Solve method for variance matrices.
#'
#' @param V Variance matrix
#' @param x Optional vector or matrix for which to solve system of equations.  If missing calculates inverse matrix.
#' @param ldV Optionally compute log determinant as well.
#' @return Matrix solving system of equations and optionally the log-determinant.
#' @details This function is faster and more stable than \code{solve} when \code{V} is known to be positive-definite.
solveV <- function(V, x, ldV = FALSE) {
  C <- chol(V) # cholesky decomposition
  if(missing(x)) x <- diag(nrow(V))
  # solve is O(ncol(C) * ncol(x)) with triangular matrices
  # using backward subsitution
  ans <- backsolve(r = C, x = backsolve(r = C, x = x, transpose = TRUE))
  if(ldV) {
    ldV <- 2 * sum(log(diag(C)))
    ans <- list(y = ans, ldV = ldV)
  }
  ans
}
