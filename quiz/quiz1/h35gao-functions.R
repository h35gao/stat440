source("mvn-functions.R")


#' Compute the conditional mean and variance
#'
#' @param mu Mean vector
#' @param Sigma Variance matrix
#' @param x1 X1
#' @param ind1 A vector of TRUE/FALSE values of length n with exactly p TRUE values.
#' @return A list with elements named cmu2 and cSigma2
cmvn <- function(mu, Sigma, x1, ind1) {
  mu1 <- mu[ind1]
  mu2 <- mu[!ind1]
  sigma11 <- Sigma[ind1,ind1,drop=FALSE]
  sigma12 <- Sigma[ind1,!ind1,drop=FALSE]
  sigma21 <- Sigma[!ind1,ind1,drop=FALSE]
  sigma22 <- Sigma[!ind1,!ind1,drop=FALSE]
  
  cmu2<- mu2 + sigma21%*%solveV(sigma11,x1-mu1)
  cSigma2 <- sigma22 - sigma21%*%solveV(sigma11)%*%sigma12
  list(cmu2,cSigma2)
}

#' Compute the variance matrix of Brownian motion Bt
#' 
#' @param tseq Time sequence
#' @return The N × N variance matrix
bmV <- function(tseq) {
  outer(tseq,tseq,"pmin")
}


#' Compute the conditional mean and standard deviations 
#'
#' @param sseq Intermediate time sequence
#' @param tseq Time sequence
#' @param Bt Brownian motion values
#' @return A list of the conditional mean and standard deviations 
cbm <- function(sseq, tseq, Bt) {
  # length of tseq
  n <- length(tseq)
  # initialize cmu and csigma
  cmu <- rep(0,n)
  csigma <- rep(0,n)
  sig2 <- (1/sseq[1]+1/(tseq[1]-sseq[1]))^(-1)
  cmu[1] <- sig2*(Bt[1]/(tseq[1]-sseq[1]))
  csigma[1] <- sqrt(sig2)
  for (i in 2:n) {
    sig2 <- (1/(sseq[i]-tseq[i-1])+1/(tseq[i]-sseq[i]))^(-1)
    cmu[i] <- sig2*(Bt[i-1]/(sseq[i]-tseq[i-1])+Bt[i]/(tseq[i]-sseq[i]))
    csigma[i] <- sqrt(sig2)
  }
  list(cmu,csigma)
}