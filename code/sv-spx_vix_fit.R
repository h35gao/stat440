#--- stochastic volatility modeling of financial assets -------------------

# download and plot GSPC and VIX data
# install.packages("quantmod")

require(quantmod)

# S&P500, VIX, and Apple
getSymbols(Symbols = c("^GSPC", "^VIX", "AAPL"),
           src = "yahoo", from = "1990-01-01")

# full ticker data
head(GSPC)

# only keep the closing prices
all(time(GSPC) == time(VIX)) # make sure dates are the same
all(time(GSPC) == time(AAPL))
# store date and strip closing price of fancy formatting
finData <- data.frame(Date = time(GSPC),
                      GSPC = as.numeric(GSPC$GSPC.Close),
                      VIX = as.numeric(VIX$VIX.Close),
                      AAPL = as.numeric(AAPL$AAPL.Close))
rownames(finData) <- NULL

head(finData)

# plot log-returns of spx and vix
N <- nrow(finData)
par(mfrow = c(1,2), mar = c(4,4,2,.5)+.1)
plot(finData$Date[1:(N-1)], diff(log(finData$GSPC)), pch = 16, cex = .5,
     xlab = "Date", ylab = expression(Delta*X[n]),
     main = expression("GSPC: Log-Returns "*(Delta*X[n])))
plot(finData$Date[1:N], finData$VIX, pch = 16, cex = .5,
     xlab = "Date", ylab = expression(V[n]),
     main = expression("VIX "*(V[n])))


#--- fit gbm model and plot goodness-of-fit test --------------------------

# plot gBM residuals and log-return time series
# the second plot shows that volatility is stochastic,
# which is why gBM model with constant volatility fails

dXt <- diff(log(finData$GSPC)) # log-returns of spx
par(mfrow = c(1,2), mar = c(4,4,2,.5)+.1)
# recall that gBM states that dXt ~iid N(mu, tau^2)
hist((dXt - mean(dXt))/sd(dXt), breaks = 100, freq = FALSE,
     main = "Goodness-of-Fit", xlab = "")
curve(dnorm, add = TRUE, col = "red")
legend("topright",
       legend = expression((Delta*X[n] - hat(mu))/hat(tau), N(0,1)),
       pch = c(22, NA), col = c("black", "red"), lty = c(NA, 1),
       seg.len = 1) # shorten line
# log-return time series
plot(finData$Date[1:(N-1)], diff(log(finData$GSPC)), pch = 16, cex = .5,
     xlab = "Date", ylab = expression(Delta*X[n]),
     main = expression("GSPC: Log-Returns "*(Delta*X[n])))


#--- estimate parameters of the SV model ---------------------------------

source("sv-functions.R")
# param names (useful for displaying things later)
theta.names <- c("alpha", "gamma", "mu", "sigma", "lambda", "tau")

# use data between 2005-2010
ind <- with(finData, {
  Date > "2005-01-01" & Date < "2011-01-01"
})

# process data to simplified form
Yt <- finData[ind,c("GSPC", "VIX")]
Yt <- as.matrix(Yt)
Yt[,1] <- log(Yt[,1]) # log-spx
colnames(Yt) <- c("X", "V")
N <- nrow(Yt)
dT <- 1/252

# plot the 2-d profile likelihood of lambda and tau
npts <- 50
lambda.seq <- seq(from = 1.1, to = 1.5, len = npts)
tau.seq <- seq(from = .001, to = .003, len = npts)
eta.mat <- expand.grid(lambda.seq, tau.seq)
system.time({
  eta.ll <- apply(eta.mat, 1, function(eta) {
    sv.profll(eta = eta, Yt = Yt, dt = dT)
  })
  eta.ll <- matrix(eta.ll, npts, npts)
})

# contour plot
par(mfrow = c(1,1), mar = c(4, 4, 2, .1)+.1)
contour(lambda.seq, tau.seq, exp(eta.ll - max(eta.ll)),
        xlab = expression(lambda), ylab = expression(tau),
#        ylim= c(.0015, .0020),
        main = expression(L[prof](lambda, tau*" | "*Y)))

#--- find mle using profile likelihood -----------------------------------------

loglik.prof <- function(eta) {
  if(any(eta < 0)) return(-Inf)
  sv.profll(eta = eta, Yt = Yt, dt = dT)
}

eta.fit <- optim(par = c(1.2, .002),
                fn = loglik.prof, control = list(fnscale = -1, maxit = 1000))

# sanity check
source("stat440-utils.R")
mle.check(loglik.prof, eta.fit$par, theta.names = expression(lambda, tau))

# full MLE
theta.mle <- sv.profmle(eta = eta.fit$par, Yt = Yt, dt = dT)
signif(theta.mle, 2) # display nicely


#--- confidence intervals -------------------------------------------------

# numerical derivatives:
# f'(x) ~ (f(x+dx) - f(x))/dx
# dx = 1/10, 1/100, 1/1000, 1/1e4

require(numDeriv) # a package for numerical derivatives

# full loglikelihood function
loglik <- function(theta) {
  # parameter restrictions: gamma, mu, sigma, lambda, tau > 0
  if(any(theta[-1] <= 0)) return(-Inf)
  sv.loglik(theta = theta, Yt = Yt, dt = dT)
}


FIobs <- -hessian(func = loglik, x = theta.mle)
# check that it's +ve def
# (1) eigen values all > 0
# Eigen decomp: V = G'DG, where D = diag, GG' = I
eigen(FIobs)$val
# (2) diagonal elements of choleski factor > 0
# Cholesky decomp: V = LL', where L is lower triangular
diag(chol(FIobs))


# full likelihood standard errors
theta.se <- sqrt(diag(solve(FIobs))) # standard errors: var(theta.mle) \approx solve(FIobs)

theta.CI <- rbind(mle = theta.mle,
                  se = theta.se,
                  L = theta.mle - 1.96 * theta.se,
                  U = theta.mle + 1.96 * theta.se)
colnames(theta.CI) <- theta.names
signif(t(theta.CI),2)

#--- profile likelihood for nuisance parameters ---------------------------

# suppose we only care about (lambda, tau)
# i.e., (alpha, gamma, mu, sigma) are "nuisance parameters"

# 1. use profile likelihood to obtain lambda.mle, tau.mle
# 2. use profile likelihood (instead of full likelihood)
#    to get confidence intervals

# profile likelihood CI
FIobs.eta <- -hessian(func = loglik.prof, x = eta.fit$par)
se.eta <- sqrt(diag(solve(FIobs.eta)))

# display standard errors
disp <- rbind(full = theta.se[5:6], prof = se.eta)
colnames(disp) <- c("lambda", "tau")
disp

#--- goodness-of-fit -----------------------------------------------------------

# look at residuals: histograms and qq-plots


# compare sv residuals to Black-Scholes (= gBm) residuals
Zbs <- diff(Yt)
Zbs <- (Zbs - mean(Zbs))/sd(Zbs) # BS residuals
Zsv <- sv.resid(Yt = Yt, dt = dT, theta = theta.mle) # SV residuals

# same plot every time
plot.res <- function(Z1, Z2, main1, main2) {
  par(mfrow = c(2,2), mar = c(3.5,3.5,1.5,.5)+.1, oma = c(0,0,0,0))
  # histograms
  hist(Z1, breaks = 100, freq = FALSE,
       xlab = "", ylab = "", main = main1)
  title(xlab = "Residuals", ylab = "Density", line = 2.5)
  curve(dnorm, col = "red", add = TRUE)
  hist(Z2, breaks = 100, freq = FALSE,
       xlab = "", ylab = "", main = main2)
  title(xlab = "Residuals", ylab = "Density", line = 2.5)
  curve(dnorm, col = "red", add = TRUE)
  # qq-plots
  qrng <- range(Z1, Z2)
  qpts <- qqnorm(Z1, type = "n", ylim = qrng, xlab = "", ylab = "")
  title(xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", line = 2.5)
  abline(a = 0, b = 1, lty = 2, col = "red")
  points(qpts, pch = 16, cex = .5)
  qpts <- qqnorm(Z2, type = "n", ylim = qrng, xlab = "", ylab = "")
  title(xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", line = 2.5)
  abline(a = 0, b = 1, lty = 2, col = "red")
  points(qpts, pch = 16, cex = .5)
}

plot.res(Z1 = Zbs, main1 = "gBm Model", Z2 = Zsv[,1], main2 = "SV Model")

# compare SV residuals for GSPC and VIX
plot.res(Z1 = Zsv[,1], main1 = "GSPC", Z2 = Zsv[,2], main2 = "VIX")

# compare GSPC residuals to AAPL residuals
# typically the VIX is a much better proxy for the volatility of SPX
# than for the volatility of individual assets

Zspx <- Zsv # save GSPC residuals

# fit SV model to AAPL
Yt <- finData[ind,c("AAPL", "VIX")]
Yt <- as.matrix(Yt)
Yt[,1] <- log(Yt[,1])
Yt[,2] <- Yt[,2]
colnames(Yt) <- c("X", "V")

# function to do all the fitting
theta.apl <- sv.fit(Yt = Yt, dt = dT, var.calc = FALSE)
Zapl <- sv.resid(Yt, dT, theta.apl)

plot.res(Z1 = Zspx[,1], main1 = "GSPC", Z2 = Zapl[,2], main2 = "AAPL")

# in this case AAPL residuals are significantly right-skewed, whereas
# SPX residuals are heavy-tailed.

# let's see how well variance of AAPL log-returns correlates with VIX dynamics
# definitely more correlation with GSPC log-returns than with those of AAPL
# also, turns out huge drop in AAPL log-ret corresponds to a stock split :)
# this is why it's critically important to include exogenous covariates
# (e.g., date of quarterly return disclosure, company acquisitions, etc) into
# the modeling framework.
N <- nrow(finData)
par(mfrow = c(1,3), mar = c(4,4,2,.5)+.1)
plot(finData$Date[1:(N-1)], diff(log(finData$GSPC)), pch = 16, cex = .5,
     xlab = "Date", ylab = expression(Delta*X[n]),
     main = expression("GSPC: Log-Returns "*(Delta*X[n])))
plot(finData$Date[1:(N-1)], diff(log(finData$AAPL)), pch = 16, cex = .5,
     xlab = "Date", ylab = expression(Delta*X[n]),
     main = expression("AAPL: Log-Returns "*(Delta*X[n])))
plot(finData$Date[1:N], finData$VIX, pch = 16, cex = .5,
     xlab = "Date", ylab = expression(V[n]),
     main = expression("VIX "*(V[n])))
