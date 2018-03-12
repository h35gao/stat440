#--- test functions for Weibull distribution -----------------------------------

source("weibull-functions.R")

#--- check log-likelihood function ---------------------------------------------

ntest <- 20 # number of datasets
ntheta <- 10 # number of likelihood evaluations per dataset

replicate(ntest, expr = {
  # simulate some data
  gamma0 <- runif(1, .5, 2)
  lambda0 <- runif(1, .5, 2)
  n <- sample(50:200, 1)
  y <- rweibull(n, shape = gamma0, scale = lambda0)
  logy <- log(y) # for loglikelihood
  # generate parameter values
  theta.seq <- cbind(gamma = runif(ntheta, .5, 2),
                     lambda = runif(ntheta, .5, 2))
  # unsimplified likelihood
  llu <- apply(theta.seq, 1, function(theta) {
    sum(dweibull(y, shape = theta[1], scale = theta[2], log = TRUE))
  })
  # simplified likelihood
  lls <- apply(theta.seq, 1, function(theta) {
    weibull.loglik(gamma = theta[1], lambda = theta[2], logy = logy)
  })
  max(abs(diff(llu-lls))) # = 0 if llu - lls = CONST
})

#--- posterior factorization ---------------------------------------------------

# simulate data
gamma0 <- 1.19
lambda0 <- 2.61
n <- 1e2
y <- rweibull(n, shape = gamma0, scale = lambda0)

# generate parameters
nsim <- 20
Eta <- rexp(nsim)
Gamma <- rexp(nsim)
# prior parameters for pi(eta | gamma) with pi(gamma) \propto 1
alpha0 <- rexp(1)
beta0 <- rexp(1)
logy <- log(y) # precompute
# log-posterior first version (joint distribution = likelihood + prior)
lp1 <- sapply(1:nsim, function(ii) {
  lp <- weibull.loglik(gamma = Gamma[ii], lambda = Eta[ii]^(-1/Gamma[ii]),
                       logy = logy) # likelihood
  lp + dgamma(Eta[ii], alpha0, beta0, log = TRUE) # prior
})
# second version (marginal + conditional)
lp2 <- sapply(1:nsim, function(ii) {
  lp <- weibull.lpcond(eta = Eta[ii], gamma = Gamma[ii],
                       logy, alpha0, beta0) # conditional
  lp + weibull.lpmarg(gamma = Gamma[ii], logy, alpha0, beta0) # marginal
})
range(diff(lp1-lp2))
