#--- mle and simplified log-likelihood for linear regression ---------------

source("stat440-utils.R")

# model is: y ~ N(X beta , sigma^2 V)

# "sufficient statistics" for calculating the MLE and profile loglikelihood
# these are: beta_hat, sigma_hat, log |V|, and n.
# if V is diagonal, supply diagV = diag(V) to accelerate things

#' Sufficient statistics for calculating the conditional MLE and profile likelihood function.
#'
#' @param y Length \code{n} vector of responses.
#' @param X \code{n x p} design matrix.
#' @param V Error correlation matrix.  Can be either (i) a dense \code{n x n} positive definite matrix, (ii) a positive vector of length \code{n}, for which \code{V} is diagonal with the corresponding elements, (iii) missing, in which case \code{V} is the identity matrix.
#' @return A list of elements from which the (conditional) MLE of \code{(beta, sigma)} given \code{y, X, V} and the profile likelihood \code{p(y | X, V, beta_hat, sigma_hat)} can be reconstructed.  Thus the elements are:
#' \describe{
#'   \item{\code{beta}}{The conditional MLE \code{beta_hat}.}
#'   \item{\code{sigma}}{The conditional MLE \code{sigma_hat}.}
#'   \item{\code{T}}{The matrix \code{t(X) %*% V^{-1} %*% X}.}
#'   \item{\code{U}}{The matrix \code{t(X) %*% V^{-1} %*% y}.}
#'   \item{\code{ldV}}{The log-determinant \code{log(det(V))}.}
#'   \item{\code{n}}{The number of observations.}
#' }
lm.suff <- function(y, X, V) {
  X <- as.matrix(X) # design matrix
  p <- ncol(X) # number of beta's
  n <- length(y) # number of observations
  if(missing(V)) {
    # V = diag(n)
    iVXy <- cbind(X,y)
    ldV <- 0
  } else if(length(V) == n) {
    # V = diag(V)
    iVXy <- cbind(X,y)/V
    ldV <- sum(log(V))
  } else if(identical(dim(V), c(n,n))) {
    # V is a dense +ve def matrix
    slv <- solveV(V = V, x = cbind(X, y), ldV = TRUE) # cholesky solve method
    iVXy <- slv$y # V^{-1}[X y]
    ldV <- slv$ldV # log(det(V))
  }
  XiVX <- crossprod(X, iVXy[,1:p]) # X'V^{-1}X
  XiVy <- crossprod(X, iVXy[,p+1]) # X'V^{-1}y
  bhat <- solve(XiVX, XiVy) # beta_hat
  iVyy <- iVXy[,p+1] - iVXy[,1:p] %*% bhat # V^{-1}(y - X bhat)
  # sigma^2_hat
  s2hat <- crossprod(y - X %*% bhat, iVyy)/n
  list(beta = drop(bhat), sigma = sqrt(drop(s2hat)),
       T = XiVX, U = drop(XiVy), ldV = ldV, n = n)
}

#' Profile likelihood for regression models.
#'
#' @param suff Sufficient statistics as calculated by \code{lm.suff}.
#' @param sigma Model standard deviation can be provided if it is known.  Usually it is unknown in which case this argument is missing.
#' @return The profile likelihood the profile likelihood \code{p(y | X, V, beta_hat, sigma_hat)}, or \code{p(y | X, V, beta_hat, sigma)}, depending on whether or not \code{sigma} is known.
lm.profll <- function(suff, sigma) {
  s2hat <- suff$sigma^2
  ldV <- suff$ldV
  n <- suff$n
  if(missing(sigma)) {
    -.5 * (n + n*log(s2hat) + ldV)
  } else {
    -.5 * (n * s2hat/sigma^2 + n * log(sigma^2) + ldV)
  }
}

## lm.suff <- function(y, X, V, diagV) {
##   X <- as.matrix(X)
##   p <- ncol(X) # number of beta's
##   n <- length(y) # number of observations
##   if(missing(diagV) && missing(V)) {
##     iVXy <- cbind(X,y)
##     ldV <- 0
##   } else if(!missing(diagV)) {
##     iVXy <- cbind(X,y)/diagV
##     ldV <- sum(log(diagV))
##   } else {
##     iVXy <- solve(V, cbind(X, y)) # V^{-1}[X y]
##     ldV <- determinant(V, log = TRUE)$mod[1]
##   }
##   XiVX <- crossprod(X, iVXy[,1:p]) # X'V^{-1}X
##   XiVy <- crossprod(X, iVXy[,p+1]) # X'V^{-1}y
##   bhat <- solve(XiVX, XiVy) # beta_hat
##   iVyy <- iVXy[,p+1] - iVXy[,1:p] %*% bhat # V^{-1}(y - X bhat)
##   # sigma^2_hat
##   s2hat <- crossprod(y - X %*% bhat, iVyy)/n
##   list(beta = drop(bhat), sigma = sqrt(drop(s2hat)),
##        T = XiVX, U = drop(XiVy), ldV = ldV, n = n)
## }
