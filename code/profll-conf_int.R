#--- profile likelihood confidence intervals -----------------------------------

require(numDeriv)
source("stat440-utils.R") # contains mle.check
source("weibull-functions.R")

# true parameter values
gamma0 <- runif(1, 1.5, 1.7) # shape
lambda0 <- runif(1, .9, 1.1) # scale

# simulate data
n <- 20
y <- rweibull(n, shape = gamma0, scale = lambda0)

# plot data and pdf
hist(y, freq = FALSE)
curve(dweibull(x, shape = gamma0, scale = lambda0), col = "red", add = TRUE)

#--- parameter inference by profile likelihood ---------------------------------

# conditional MLE given gamma
lambda.hat <- function(gamma, logy) mean(exp(gamma * logy))^(1/gamma)

# profile likelihood function
# data input is logy = log(y) for computational efficiency
weibull.profll <- function(gamma, logy) {
  weibull.loglik(gamma = gamma,
                 lambda = lambda.hat(gamma, logy), logy = logy)
}


# calculate MLE
logy <- log(y)
# 1d optimizer requires interval on which to find min/max
# in practice we would need to figure out how to do this...
gamma.rng <- c(.1, 5) # this should be sufficient given theta_true
gamma.mle <- optimize(f = weibull.profll, interval = c(.1, 5),
                      logy = logy,
                      # instead of f = function(gamma) profll(gamma, logy)
                      maximum = TRUE)$maximum
lambda.mle <- lambda.hat(gamma.mle, logy)
theta.mle <- c(gamma.mle, lambda.mle)

mle.check(loglik = function(theta) {
  if(any(theta < 0)) return(-Inf)
  weibull.loglik(theta[1], theta[2], logy)
}, theta.mle = theta.mle, theta.names = expression(gamma, lambda))

#--- confidence intervals ------------------------------------------------------

# full likelihood function with theta = c(gamma, lambda)
wloglik <- function(theta, logy) weibull.loglik(theta[1], theta[2], logy)


# full likelihood se for gamma
se.full <- hessian(func = wloglik, x = theta.mle, logy = logy)
se.full <- sqrt(diag(-solve(se.full))[1])

# profile likelihood se for gamma
se.prof <- hessian(func = weibull.profll, x = gamma.mle, logy = logy)
se.prof <- sqrt(-1/se.prof[1])

diff(c(full = se.full, prof = se.prof))

#--- scratch -------------------------------------------------------------------

loglik <- function(theta) {
  if(any(theta < 0)) return(-Inf)
  sum(dweibull(y, shape = theta[1], scale = theta[2], log = TRUE))
}

lambda.hat <- function(gamma) mean(y^gamma)

profll <- function(gamma) loglik(c(gamma, lambda.hat(gamma)))

gamma.seq <- seq(.5, 2, len = 200)
gamma.ll <- sapply(gamma.seq, profll)

plot(gamma.seq, gamma.ll, type = "l")

# calculate MLE
gamma.mle <- optimize(profll, interval = c(.5, 2),
                      maximum = TRUE)$maximum
lambda.mle <- lambda.hat(gamma = gamma.mle)
theta.mle <- c(gamma.mle, lambda.mle)
