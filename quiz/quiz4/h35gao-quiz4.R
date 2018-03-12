source('ggamma-dpqr_functions.R')

set.seed(2018)

# Q1
ggamma.loglik <- function(alpha, eta, lambda, logy) {
  n <- length(logy)
  S <- sum(logy)
  Tlambda <- sum(exp(lambda*logy))
  shape <- alpha/lambda
  result <- n * (log(lambda) - lgamma(shape) + shape*log(eta)) + alpha*S - eta*Tlambda
  result
}

# Q3
ggamma.logmarg <- function(alpha, lambda, logy, khat, ghat) {
  n <- length(logy)
  S <- sum(logy)
  shape <- alpha/lambda
  result <- n*(log(lambda)-lgamma(shape))+alpha*S-khat*log(ghat)+lgamma(khat)
  result
}

# Q4
# marginal prior distribution for (alpha, lambda)
logprior <- function(alpha, lambda) {
  dlnorm(alpha, meanlog = 0, sdlog = 2, log = TRUE) + 
    dlnorm(lambda, meanlog = 0, sdlog = 2, log = TRUE)
}

ggamma.post <- function(nsamples, y, alpha0, lambda0, burn, mwg.sd, acc.out = FALSE) {
  logy <- log(y)
  n <- length(logy)
  if(missing(burn)) burn <- min(floor(nsamples/10), 1e3) # burnin
  Theta <- matrix(NA, nsamples, 3)
  colnames(Theta) <- c("alpha", "lambda", "eta")
  ntheta <- ncol(Theta) # number of parameters
  # keep MH acceptance probability, one for each parameter
  paccept <- rep(0, ntheta)
  names(paccept) <- colnames(Theta)
  # initialize
  khat <- n*alpha0/lambda0+1
  ghat <- sum(exp(lambda0*logy))
  eta0 <- rgamma(1,khat,ghat)
  theta.curr <- cbind(alpha0,lambda0,eta0)
  # only evaluate the posterior once per acceptance rate calculation
  marg.curr <- ggamma.logmarg(alpha0,lambda0,logy,khat,ghat) + logprior(alpha0,lambda0)
  loglik.curr <- ggamma.loglik(alpha0,eta0,lambda0,logy)
  for(ii in (-burn+1):nsamples) {
    # update alpha and lambda
    for(jj in 1:ntheta-1) {
      theta.prop <- theta.curr
      theta.prop[jj] <- theta.prop[jj] + mwg.sd[jj] * rnorm(1)
      khat <- n*theta.prop[1]/theta.prop[2]+1
      ghat <- sum(exp(theta.prop[2]*logy))
      if(all(theta.prop > 0)) {
        # only bother updating if draw is "valid"
        marg.prop <- ggamma.logmarg(theta.prop[1], theta.prop[2], logy, khat, ghat) + logprior(theta.prop[1], theta.prop[2])
        lacc <- marg.prop - marg.curr # log acceptance rate
        if(lacc > 0 || runif(1) < exp(lacc)) {
          # automatic accept if acc = exp(lacc) > 1
          theta.curr <- theta.prop
          marg.curr <- marg.prop
          paccept[jj] <- paccept[jj]+1
        }
      }
    }
    # update eta
    theta.prop <- theta.curr
    khat <- n*theta.prop[1]/theta.prop[2]+1
    ghat <- sum(exp(theta.prop[2]*logy))
    theta.prop[3] <- rgamma(1,khat,ghat)
    if(all(theta.prop > 0)) {
      # only bother updating if draw is "valid"
      loglik.prop <- ggamma.loglik(theta.prop[1], theta.prop[3], theta.prop[2], logy)
      lacc <- loglik.prop - loglik.curr # log acceptance rate
      if(lacc > 0 || runif(1) < exp(lacc)) {
        # automatic accept if acc = exp(lacc) > 1
        theta.curr <- theta.prop
        loglik.curr <- loglik.prop
        paccept[3] <- paccept[3]+1
      }
    }
    # storage
    if(ii > 0) Theta[ii,] <- theta.curr
  }
  paccept <- paccept/(nsamples+burn) # acceptance rate
  message("MH Acceptance Rate:")
  message("gamma: ", round(paccept[1]*100), "%")
  message("lambda: ", round(paccept[2]*100), "%")
  message("eta: ", round(paccept[3]*100), "%")
  if(acc.out) {
    ans <- list(Theta = Theta, accept = paccept)
  } else ans <- Theta
  ans
}