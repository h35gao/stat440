#--- maximum likelihood for stochastic volatility model ------------------------

# model:
# dXt = (alpha - (tau*Vt)/2)dt + (tau*Vt)^(1/2) dBt
# dVt = -gamma(Vt - mu)dt + sigma Vt^lambda dBt

# Yt = (Xt, Vt)
# theta = (alpha, gamma, mu, sigma, lambda, tau)

source("sv-functions.R")
source("stat440-utils.R")

#--- simulate data --------------------------------------------------------

# model parameters
alpha <- .1
gamma <- 5
mu <- .2
sigma <- .3
lambda <- .5
tau <- 2
theta <- c(alpha = alpha, gamma = gamma, mu = mu,
           sigma = sigma, lambda = lambda, tau = tau)

dt <- 1/252 # 252 trading days in a year
# initial values
X0 <- 6.5
V0 <- .2

n <- 500 # number of observations

# simulated data
Yt <- sv.sim(n = n, Y0 = c(X0, V0), dt = dt, theta = theta)

par(mfrow = c(1,2), # 1 row, 2 columns
    mar = c(4,2.5,2,.1)+.1)
ylab <- expression(X[t], V[t])
cex <- 1
for(ii in 1:2) {
plot(Yt[,ii], pch = 16, cex = cex,
     xlab = "Date", ylab = "", main = ylab[ii])
}

#--- maximum likelihood estimation -----------------------------------------

# done with "optim" function, which takes two main arguments:
# (1) loglik(theta), the function to maximize
# (2) theta0, the starting value

loglik <- function(theta) {
  # parameter restrictions: gamma, mu, sigma, lambda, tau > 0
  if(any(theta[-1] <= 0)) return(-Inf)
  sv.loglik(theta = theta, Yt = Yt, dt = dt)
}

# wouldn't know these in practice, but good enough for testing
theta0 <- theta

# find mle
# last argument means maximize instead of minimize
mle.fit <- optim(fn = loglik, par = theta0,
                 control = list(fnscale = -1, maxit = 1000))

# any number except 0 means failed
mle.fit$convergence

# check results
theta.mle <- mle.fit$par
rbind(true = theta0, est = theta.mle)

# check if estimate is a mode of loglikelihood
# <=> each component is a mode of univariate function with all other
#     components fixed

theta.names <- c("alpha", "gamma", "mu", "sigma", "lambda", "tau")
theta.names <- parse(text = theta.names) # use greek letters

mle.check(loglik = loglik,
          theta.mle = theta.mle,
          theta.names = theta.names,
          layout = c(2,3))


# works well when starting value is close to true parameter value,
# but not so well otherwise...
mle.fit2 <- optim(fn = loglik,
                  par = 10*theta0, # bad starting value
                  control = list(fnscale = -1, maxit = 1000))
mle.fit2$convergence # optim thinks it has  converged...

# ...but maximum check fails miserably
mle.check(loglik = loglik,
          theta.mle = mle.fit2$par,
          theta.names = theta.names,
          layout = c(2,3))


#--- profile likelihood ----------------------------------------------------

# plot the 2-d profile likelihood of lambda and tau
npts <- 100
lambda.seq <- seq(from = 0, to = 2, len = npts)
tau.seq <- seq(from = 1, to = 2.5, len = npts)
eta.mat <- expand.grid(lambda.seq, tau.seq)
system.time({
  eta.ll <- apply(eta.mat, 1, function(eta) {
    sv.profll(eta = eta, Yt = Yt, dt = dt)
  })
  eta.ll <- matrix(eta.ll, npts, npts)
})

par(mfrow = c(1,1), mar = c(4, 4, 2, .1)+.1)
contour(lambda.seq, tau.seq, exp(eta.ll - max(eta.ll)),
        xlab = expression(lambda), ylab = expression(tau),
        drawlabels = FALSE,
        main = expression("\u2113"[prof](lambda, tau*" | "*bold(Y))))
points(lambda, tau, pch = 3, lwd = 3, cex = 3, col = "red")
legend("topright", legend = expression((list(lambda["true"], tau["true"]))),
       col = "red", pch = 3, pt.lwd = 2)

#--- find mle using profile likelihood -----------------------------------------

# step 1: optimize over reduced parameters eta = (lambda, tau)
eta.fit <- optim(par = c(1.3, .0015),
                fn = function(eta) {
                  if(any(eta < 0)) return(-Inf)
                  sv.profll(eta = eta, Yt = Yt, dt = dt)
                }, control = list(fnscale = -1, maxit = 1000))

# step 2: obtain remaining parameters analytically via
# phi_hat(eta) = (alpha, gamma, mu, sigma)
theta.mle <- sv.profmle(eta = eta.fit$par, Yt = Yt, dt = dt)

# check that this is indeed the mode
mle.check(loglik = loglik,
          theta.mle = theta.mle,
          theta.names = theta.names,
          layout = c(2,3))
