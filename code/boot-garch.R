#--- bootstrapping for GARCH forecasts ------------------------------------

# get apple data
require(quantmod)
getSymbols(Symbols = "AAPL",
           src = "yahoo", from = "1990-01-01")

# mar and oma can be used to remove unwanted margins
# and give plots common axes
par(mfrow = c(1,2),
    mar = c(3,4,2,.2)+.5, oma = c(2,0,0,1))
# full dataset
plot(AAPL$AAPL.Close, type = "l",
     xlab = "", main = "Apple Inc. Stock Prices 1990-2016",
     ylab = "AAPL")
# restricted range for inference
dates <- c("2006-01-01", "2012-01-01")
St <- AAPL[paste0(dates, collapse = "/")]$AAPL.Close
plot(diff(log(St)), type = "l",
     xlab = "",
     main = "Log-Returns between 2006-2012",
     ylab = expression(Delta*log("AAPL")))
mtext(text = "Date", cex = 1.2,
      outer = TRUE, side = 1, line = .5, font = 2)

#--- parameter estimation -------------------------------------------------

# this uses the garch11 package
# to install it:
# (1) download the file garch11_0.0.0.9001.tar.gz
# (2) open R in the folder where you saved it
# (3) run the command:
install.packages("garch11_0.0.0.9001.tar.gz",
                 repos = NULL, type = "source")

require(garch11)
source("stat440-utils.R")

# convert time series object to ordinary vector
# typically assume daily data even if wknds/holidays are excluded
St[1:10,]
# xts format doesn't shorten when you diff...
dXt <- diff(log(St))[-1] # log returns
# convert to simple numeric format for fitting functions
dXt <- as.numeric(dXt)
# often good to check with time series (as sometimes wknds are NA's)
anyNA(dXt)

mu.hat <- mean(dXt) # subtract the mean (and add it back in later)
eps <- dXt-mu.hat

# MLE
# automatic initialization
theta.mle <- garch.fit(eps = eps)

# check
# loglikelihood function
loglik <- function(theta) {
  garch.loglik(omega = theta[1], alpha = theta[2], beta = theta[3],
               eps = eps,
               sig20 = theta[1])
}
tnames <- expression(omega, alpha, beta)

mle.check(loglik = loglik, theta.mle = theta.mle,
          refit = TRUE, theta.names = tnames) # pass

# variance of stationary distribution (want this to be positive)
theta.mle["omega"]/(1 - theta.mle["alpha"] - theta.mle["beta"])

#--- confidence intervals for parameters ----------------------------------

require(numDeriv)

# observed Fisher information
FI.obs <- -hessian(x = theta.mle,
                   func = function(theta) {
                     omega <- theta[1]
                     alpha <- theta[2]
                     beta <- theta[3]
                     garch.loglik(omega, alpha, beta, eps, sig20 = omega)
                   })

FI.obs
# numerical 2nd derivative algorithm is unstable
# as omega.hat close to boundary = 0
# so instead calculate on the log-scale

FI.obs <- -hessian(x = log(theta.mle),
                   func = function(ltheta) {
                     omega <- exp(ltheta[1])
                     alpha <- exp(ltheta[2])
                     beta <- exp(ltheta[3])
                     garch.loglik(omega, alpha, beta, eps, sig20 = omega)
                   })

FI.obs # now works nicely

# confidence interval calculation
# first do it on log scale
ltheta.se <- sqrt(diag(solve(FI.obs))) # standard errors of log(theta)
ltheta.CI <- rbind(est = log(theta.mle),
                   lwr = log(theta.mle) - 1.96 * ltheta.se,
                   upr = log(theta.mle) + 1.96 * ltheta.se)
# now report on regular scale
theta.CI <- exp(ltheta.CI)
signif(theta.CI, 2)

# let's see how asymmetric these CI's are:
# calculate (est - L)/(est - U)
signif((theta.CI["est",]-theta.CI["lwr",])/
       (theta.CI["upr",]-theta.CI["lwr",]), 2)
# the CI for omega is somewhat left-skewed due to the boundary effect

# confidence interval for the probability of an
# gamma x 100% _drop_ in asset value tomorrow
garch.risk <- function(gamma, omega, alpha, beta, eps0, sig0) {
  # log return required for alpha change, i.e. (S1-S0)/S0 = alpha
  eps1 <- log(1 - gamma)
  sig1 <- sqrt(omega + alpha * eps0^2 + beta * sig0^2)
  pnorm(eps1, sd = sig1)
}

gamma <- .02 # 2% drop
# parameter estimates
omega.hat <- theta.mle["omega"]
alpha.hat <- theta.mle["alpha"]
beta.hat <- theta.mle["beta"]
# to forecast, need values of (sig_N, eps_N) of original time series
N <- length(eps) # length of time series
sig2.hat <- garch.sig2(omega = omega.hat, alpha = alpha.hat,
                       beta = beta.hat, eps = eps)
epsN <- eps[N]
sigN <- sqrt(sig2.hat[N])

Jac <- jacobian(function(ltheta) {
  garch.risk(gamma, exp(ltheta[1]), exp(ltheta[2]), exp(ltheta[3]),
             eps0 = epsN, sig0 = sigN)
}, x = log(theta.mle))

tau.mle <- garch.risk(gamma, omega.hat, alpha.hat, beta.hat, epsN, sigN)
tau.se <- sqrt(Jac %*% solveV(FI.obs, t(Jac)))[1]
signif(c(est = tau.mle, se = tau.se), 2)



#--- bootstrapping the uncertainty of value-at-risk -----------------------

# value at risk estimates
nfwd <- 5 # number of days to forecast
nsim <- 1e5 # number of fwd simulations used to estimate quantile
VaR.tau <- c(.01, .05, .1) # VaR levels to calculate

# parameter estimates
omega.hat <- theta.mle["omega"]
alpha.hat <- theta.mle["alpha"]
beta.hat <- theta.mle["beta"]
# to forecast, need values of (sig_N, eps_N) of original time series
N <- length(eps) # length of time series
sig2.hat <- garch.sig2(omega = omega.hat, alpha = alpha.hat,
                       beta = beta.hat, eps = eps)
epsN <- eps[N]
sigN <- sqrt(sig2.hat[N])
# also need last asset value
## SN <- St[N+1] # N returns -> N+1 values

# using z~iid N(0,1)
VaR.norm <- garch.VaR(nfwd = nfwd, nsim = nsim, tau = VaR.tau,
                      omega = omega.hat, alpha = alpha.hat,
                      beta = beta.hat,
                      eps0 = epsN, sig0 = sigN)

# re-introduce mu
garch.addmu <- function(VaR, nfwd, mu) {
  Y <- log(VaR + 1) # eps_1 + ... + eps_nfwd
  Y <- Y + nfwd * mu # add mu to each log-return
  exp(Y) - 1 # convert back to (S_nfwd - S_0)/S_0
}

VaR.norm <- garch.addmu(VaR.norm, nfwd = nfwd, mu = mu.hat)
signif(-VaR.norm, 2) # VaR (+ve means drop in value)

#--- relaxed model assumptions -------------------------------------------------

# original GARCH model:
#
# eps_t = sig_t * z_t,   z_t ~iid N(0,1)
# sig^2_t^2 = omega + alpha * eps^2_t-1 + beta * sig^2_t-1
#
# relaxed model:  z_t ~iid F(z), for some unspecified distribution F.
#
# VaR estimate obtained by resampling z.hat from its empirical distribution

# using empirical distribution of z
z.hat <- eps/sqrt(sig2.hat) # estimated innovations

hist(z.hat, breaks = 50, freq = FALSE,
     main = "Empirical vs Theoretical Innovations",
     xlab = expression(z))
curve(dnorm, col = "red", add = TRUE, lwd = 2)
legend("topleft", legend = expression(hat(z), N(0,1)),
       fill = c("black", "red"))

# VaR Calculation

# resample innovations
z.emp <- sample(z.hat, nfwd*nsim, replace = TRUE)
z.emp <- matrix(z.emp, nrow = nfwd, ncol = nsim)
# var calculation
VaR.emp <- garch.VaR(nfwd = nfwd, nsim = nsim, tau = VaR.tau,
                     omega = omega.hat, alpha = alpha.hat,
                     beta = beta.hat,
                     eps0 = epsN, sig0 = sigN, z = z.emp)
VaR.emp <- garch.addmu(VaR.emp, nfwd = nfwd, mu = mu.hat)

# the 1% empirical VaR is lower, but the 5% and 10% are higher.
# this reflects the fact that the empirical distribution is almost N(0,1),
# but with a couple of negative outliers.
VaR.est <- rbind(norm = VaR.norm, emp = VaR.emp)
signif(-VaR.est, 2)


#--- parallelized bootstrap on multicore system ---------------------------

# set up cluster
require(doParallel) # simplest way of parallelizing on all platforms
cl <- makeCluster(spec = 4) # 4 cores
registerDoParallel(cl) # synch with doParallel package
clusterSetRNGStream(cl = cl) # RNG for multiple cores


nboot <- 1e3 # number of bootstrap samples

# parametric bootstrap, z: N(0,1)
# need to tell each R worker to load garch11 package
system.time({
  VaR.bootP <- foreach(ii=1:nboot,
                            .combine = "cbind",
                            .packages = "garch11") %dopar% {
    # simulate time series
    eps.boot <- garch.sim(nobs = N, nts = 1,
                          omega = omega.hat, alpha = alpha.hat,
                          beta = beta.hat,
                          eps0 = eps[1])
    Yt.boot <- eps.boot + mu.hat # don't forget mu.hat!
    # parameter estimates (starting value for optim is MLE)
    mu.boot <- mean(Yt.boot)
    theta.boot <- garch.fit(eps = Yt.boot - mu.boot,
                            eta0 = alpha.hat/omega.hat,
                            beta0 = beta.hat)
    # VaR calculation
    CI.boot <- garch.VaR(nfwd = nfwd, nsim = nsim,
                         tau = VaR.tau,
                         omega = theta.boot["omega"],
                         alpha = theta.boot["alpha"],
                         beta = theta.boot["beta"],
                         eps0 = epsN, sig0 = sigN)
    CI.boot <- garch.addmu(CI.boot, nfwd, mu.boot)
    CI.boot # return
  }
})

#--- residual bootstrap --------------------------------------------------------

# This is for z_t ~iid F(z) with unspecified F(z)

system.time({
  VaR.bootR <- foreach(ii=1:nboot,
                       .combine = "cbind",
                       .packages = "garch11") %dopar% {
    omega.hat <- theta.mle["omega"]
    alpha.hat <- theta.mle["alpha"]
    beta.hat <- theta.mle["beta"]
    # simulate time series
    # sample residuals with replacement
    z.boot <- sample(z.hat, replace = TRUE)
    eps.boot <- garch.sim(nobs = N, nts = 1, z = z.boot,
                          omega = omega.hat, alpha = alpha.hat,
                          beta = beta.hat,
                          eps0 = eps[1])
    Yt.boot <- eps.boot + mu.hat
    # parameter estimates
    mu.boot <- mean(Yt.boot)
    theta.boot <- garch.fit(Yt.boot-mu.boot,
                            eta0 = alpha.hat/omega.hat,
                            beta0 = beta.hat)
    # VaR calculation
    z.emp <- sample(z.hat, size = nfwd*nsim, replace = TRUE)
    z.emp <- matrix(z.emp, nrow = nfwd, ncol = nsim)
    CI.boot <- garch.VaR(nfwd = nfwd, nsim = nsim, tau = VaR.tau,
                         omega = theta.boot["omega"],
                         alpha = theta.boot["alpha"],
                         beta = theta.boot["beta"],
                         eps0 = epsN, sig0 = sigN, z = z.emp)
    CI.boot <- garch.addmu(CI.boot, nfwd, mu.boot)
    CI.boot
  }
})

# NOTE: NEVER save the whole work space
# otherwise, you will restore the random seed
## save(VaR.norm, VaR.emp, VaR.bootP, VaR.bootR,
##      file = "garch-boot_sim3.RData")

# NOTE: ALWAYS remember to de-allocate cluster resources!
stopCluster(cl)

#--- display bootstrap confidence intervals ---------------------------------

load("garch-boot_sim3.RData")

# basic bootstrap intervals
boot.CI <- function(theta.hat, theta.boot, conf = .95) {
  conf <- 1-conf
  conf <- c(conf/2, 1-conf/2)
  quantile(theta.hat - theta.boot, probs = conf) + theta.hat
}

# in general, try to avoid as much copy-pasting as possible
VaR.CI.disp <- function(VaR.boot, VaR.est) {
  VaR.CI <- sapply(1:length(VaR.tau), function(ii) {
    boot.CI(VaR.est[ii], VaR.boot[ii,])
  })
  VaR.CI <- rbind(pt_est = VaR.est, VaR.CI,
                  CI_width = apply(VaR.CI, 2, diff))
  rownames(VaR.CI)[2:3] <- c("L_2.5%", "U_97.5%")
  names(dimnames(VaR.CI)) <- c("", "VaR_Level")
  VaR.CI
}

# parametric bootstrap
VaR.CI.P <- VaR.CI.disp(VaR.bootP, VaR.norm)

# residual bootstrap
VaR.CI.R <- VaR.CI.disp(VaR.bootR, VaR.emp)

# compare the two
lapply(list(Parametric = VaR.CI.P, Residual = VaR.CI.R),
       signif, digits = 2)

