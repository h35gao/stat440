#--- profiling parameters in regression models ----------------------------

source("proflm-functions.R")
source("stat440-utils.R") # for mle.check
source("mvn-functions.R") # for density of multivariate normal

# generate a matrix of iid N(0,1)
rMnorm <- function(n, p) {
  if(missing(p)) p <- n
  matrix(rnorm(n*p), n, p)
}

#--- check that MLE maximizes the loglikelihood --------------------------------

# loglikelihood function
lm.loglik <- function(beta, sigma, y, X, V) {
  if(missing(V)) V <- diag(length(y)) # V is identity matrix
  if(length(V) == length(y)) {
    # V is diagonal
    ans <- sum(dnorm(y, mean = X %*% beta, sd = sigma * sqrt(V), log = TRUE))
  } else {
    # V is dense +ve def matrix
    ans <- dmvn(x = y, mu = c(X %*% beta), V = sigma^2 * V, log = TRUE)
  }
  ans
}

# reduced version for mle.check
loglik <- function(theta) {
  if(theta[p+1] < 0) return(-Inf)
  lm.loglik(beta = theta[1:p], sigma = theta[p+1], y = y, X = X, V = V)
}

# simulate data (probably want to do this multiple times)
n <- sample(10:20, 1)
p <- sample(1:4, 1)
X <- rMnorm(n, p)
beta <- rnorm(p)
sigma <- rexp(1)
# parameter names
theta.names <- c(paste0("beta[", 1:p, "]"), "sigma")
theta.names <- parse(text = theta.names)

# case 1: V is identity matrix
V <- diag(n)
y <- rmvn(1, mu = X %*% beta, V = sigma^2 * V)[1,]
suff <- lm.suff(y = y, X = X)
mle.check(loglik = loglik, theta.mle = c(suff$beta, suff$sigma),
          theta.names = theta.names)

# case 2: V is diagonal matrix
V <- rexp(n)
y <- rmvn(1, mu = X %*% beta, V = sigma^2 * diag(V))[1,]
suff <- lm.suff(y = y, X = X, V = V)
mle.check(loglik = loglik, theta.mle = c(suff$beta, suff$sigma),
          theta.names = theta.names)

# case 3: V is dense +ve def matrix
V <- crossprod(rMnorm(n))
y <- rmvn(1, mu = X %*% beta, V = sigma^2 * V)[1,]
suff <- lm.suff(y = y, X = X, V = V)
mle.check(loglik = loglik, theta.mle = c(suff$beta, suff$sigma),
          theta.names = theta.names)

#--- check that profile likelihood is likelihood at (conditional) MLE ----------

# case 1: V is identity matrix
V <- diag(n)
y <- rmvn(1, mu = X %*% beta, V = sigma^2 * V)[1,]
suff <- lm.suff(y = y, X = X)
# sigma unknown
ll.full <- loglik(theta = c(suff$beta, suff$sigma))[1]
ll.prof <- lm.profll(suff = suff) - n/2*log(2*pi) # calc drops 2*pi term
ll.full - ll.prof
# sigma known
ll.full <- loglik(theta = c(suff$beta, sigma))[1]
ll.prof <- lm.profll(suff = suff, sigma = sigma) - n/2*log(2*pi)
ll.full - ll.prof

# case 2: V is diagonal matrix
V <- rexp(n)
y <- rmvn(1, mu = X %*% beta, V = sigma^2 * diag(V))[1,]
suff <- lm.suff(y = y, X = X, V = V)
# sigma unknown
ll.full <- loglik(theta = c(suff$beta, suff$sigma))[1]
ll.prof <- lm.profll(suff = suff) - n/2*log(2*pi)
ll.full - ll.prof
# sigma known
ll.full <- loglik(theta = c(suff$beta, sigma))[1]
ll.prof <- lm.profll(suff = suff, sigma = sigma) - n/2*log(2*pi)
ll.full - ll.prof

# case 3: V is dense +ve def matrix
V <- crossprod(rMnorm(n))
y <- rmvn(1, mu = X %*% beta, V = sigma^2 * V)[1,]
suff <- lm.suff(y = y, X = X, V = V)
# sigma unknown
ll.full <- loglik(theta = c(suff$beta, suff$sigma))[1]
ll.prof <- lm.profll(suff = suff) - n/2*log(2*pi)
ll.full - ll.prof
# sigma known
ll.full <- loglik(theta = c(suff$beta, sigma))[1]
ll.prof <- lm.profll(suff = suff, sigma = sigma) - n/2*log(2*pi)
ll.full - ll.prof
