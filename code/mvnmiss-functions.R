#--- EM for a simple multivariate normal model -----------------------------

require(mvtnorm)
source("stat440-utils.R") # for solveV

# the model is:
# y = (x,z) ~ N(0, V)
# x is always observed, but z is sometimes missing

#' Loglikelihood function.
#'
#' @param Y1 Matrix of complete observations
#' @param X0 Matrix of observations with \code{z} missing.
#' @param V Variance matrix.
#' @return Observed data loglikelihood.
mvn0.loglik <- function(V, Y1, X0) {
  px <- ncol(X0)
  sum(dmvnorm(x = Y1, sigma = V, log = TRUE)) +
    sum(dmvnorm(x = X0, sigma = V[1:px,1:px,drop=FALSE], log = TRUE))
}

#' MLE of zero-mean multivariate normal with partially missing observations.
#'
#' @param Y1 matrix of complete observations
#' @param X0 matrix of observations with z missing.
#' @param tol Relative tolerance for algorithm stopping time.
#' @param niter Maximum number of iterations for algorithm stopping time.
#' @details Estimation is done with EM algorithm.
#' @return The MLE of \code{V}.
mvn0.fit <- function(Y1, X0, tol = 1e-6, niter = 100) {
  if(!is.matrix(Y1) | !is.matrix(X0)) stop("Y1 and X0 must be matrices.")
  tX0 <- t(X0)
  tYY1 <- crossprod(Y1)
  # problem dimensions
  px <- ncol(X0)
  py <- ncol(Y1)
  pz <- py-px
  if(pz <= 0)
    stop("Complete data must have higher dimension than incomplete data.")
  n0 <- nrow(X0)
  n1 <- nrow(Y1)
  N <- n0+n1
  # index sequences
  ix <- 1:px
  iz <- px+1:pz
  i0 <- 1:n0
  i1 <- n0+1:pz
  # initial estimate based on complete observations only
  Vo <- crossprod(Y1)/n1
  for(ii in 1:niter) {
    # --- E-step ---
    # Vxx^{-1} [X0' Vxz]
    iV <- solveV(Vo[ix,ix,drop=FALSE],
                 cbind(tX0, Vo[ix,iz,drop=FALSE]))
    ziV <- crossprod(Vo[ix,iz,drop=FALSE], iV)
    # (Vzx Vxx^{-1} X0')'
    mu <- t(ziV[,i0,drop=FALSE])
    # Vzz - Vzx Vxx^{-1} Vxz
    W <- Vo[iz,iz,drop=FALSE] - ziV[,i1,drop=FALSE]
    T <- crossprod(cbind(X0, mu))
    T[iz,iz] <- T[iz,iz] + n0*W
    # --- M-step ---
    Vn <- (tYY1 + T)/N
    # relative tolerance
    re <- 2*max(abs((Vo-Vn)/(Vo+Vn)))
    Vo <- Vn
    if(re < tol) break
  }
  if(ii == niter && re > tol) warning("EM algorithm did not converge.")
  Vn
}
