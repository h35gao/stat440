#--- check code for hlm fitting --------------------------------------------

source("hlm-functions.R")
source("stat440-utils.R")

# simulate data
n <- 100 # number of observations
px <- 2 # number of beta's
pw <- 3 # number of gamma's
X <- cbind(1, matrix(rnorm(n*(px-1)), n, px-1))
W <- cbind(1, matrix(rnorm(n*(pw-1)), n, pw-1))
# true parameter values
beta0 <- rnorm(px)
gamma0 <- rnorm(pw)
y <- rnorm(n, mean = X %*% beta0, sd = exp(W %*% (gamma0/2)))

theta.hat <- hlm.fit(y = y, X = X, W = W)

# quick check
rbind(true = beta0, est = theta.hat$beta) # beta
rbind(true = gamma0, est = theta.hat$gamma) # gamma

loglik <- function(theta) {
  hlm.loglik(beta = theta[1:px], gamma = theta[px+1:pw],
             y = y, X = X, W = W)
}

theta.names <- c(paste0("beta[", 1:px-1, "]"),
                 paste0("gamma[", 1:pw-1, "]"))
theta.names <- parse(text = theta.names) # convert to greek symbols

# thorough check
mle.check(loglik = loglik,
          theta.mle = c(theta.hat$beta, theta.hat$gamma),
          theta.names = theta.names, layout = c(2,3))

