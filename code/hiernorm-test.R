# test 8 schools functions

source("hiernorm-functions.R")

# test posterior factorization:
# log p(lambda, tau | x) = log p(lambda | tau, x) + log p(tau | x) + CONST
# for any value of theta = (lambda, tau)

# generate a dataset
K <- 8
x <- rnorm(K, sd = 5)
sig2 <- rpois(K, lambda = 10)

ntest <- 10 # number of tests
lambda <- rnorm(ntest) # generate parameters
tau <- rexp(ntest)
# joint log posterior (prior is pi(lambda, tau) \propto 1)
lp1 <- sapply(1:ntest, function(ii) {
  loglik.theta(lambda = lambda[ii], tau = tau[ii], x = x, sig2 = sig2)
})
# conditional + marginal
lp2 <- sapply(1:ntest, function(ii) {
  lpcond.lambda(lambda = lambda[ii], tau = tau[ii], x = x, sig2 = sig2) +
    lpmarg.tau(tau = tau[ii], x = x, sig2 = sig2)
})
lp1 - lp2 # should different by a constant
