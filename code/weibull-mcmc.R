#--- illustration of basic MCMC techniques using the Weibull model -------------

source("weibull-mcmc_functions.R")
require(RColorBrewer) # nice color palette

#--- some examples of weibull distribtuion -------------------------------------

# various shape parameters
Shape <- c(.8, 1, 2, 4)
npar <- length(Shape)
Scale <- rep(1, npar)

# plot the pdfs
x <- seq(1e-4, 4, len = 1e3)
Y <- sapply(1:npar, function(ii) {
  y <- dweibull(x = x, shape = Shape[ii], scale = Scale[ii], log = TRUE)
  exp(y - max(y)) # each has max 1
})
par(mfrow = c(1,1), mar = c(4.5,4.5,2,.5))
plot(0, type = "n", xlim = range(x), ylim = range(Y),
     xlab = expression(y), ylab = expression(f(y)),
     main = expression(Y %~% "Weibull"*(list(gamma,1))))
clrs <- brewer.pal(n = npar, name = "Dark2")
# without invisible command, returns a list of NULLs
invisible({
  sapply(1:npar, function(ii) {
    lines(x, Y[,ii], col = clrs[ii], lwd = 3)
  })
})
legend("topright",
       legend = parse(text = paste0("gamma ==", Shape)),
       fill = clrs)


#--- posterior inference -------------------------------------------------------

# pick a specific set of parameters and simulate some data
# these parameters correspond to fitted wind speed velocities
# at an observatory in Thailand.  Source: Waewsak et al (2011).
gamma0 <- 1.19
lambda0 <- 2.61
n <- 1e2
y <- rweibull(n, shape = gamma0, scale = lambda0)

# plot
par(mfrow = c(1,1), mar = c(4.2,4,1.5,1))
# mode of the weibull distribution (to determine y-axis)
md <- lambda0 * ((gamma0-1)/gamma0)^(1/gamma0)
# histogram height
hh <- max(hist(y, breaks = 25, plot = FALSE)$density)
mx <- max(dweibull(md, shape = gamma0, scale = lambda0), hh)
main <- paste0("list(gamma == ", gamma0, ", lambda == ", lambda0, ")")
hist(y, breaks = 25, freq = FALSE,
     main = parse(text = main),
     xlab = "y", ylab = "Density", ylim = c(0, mx))
curve(dweibull(x, shape = gamma0, scale = lambda0),
      add = TRUE, col = "red")


#--- analytic posterior --------------------------------------------------------

# for 2D posterior distributions, grid-based approximation is very fast and
# has almost no error, so we call this the analytic posterior.

# create the grid on which to evaluate the posterior
npts <- 100 # number of points on each axis
gseq <- seq(.75, 1.75, len = npts) # gamma grid
lseq <- seq(1.5, 3.5, len = npts) # lambda grid
# stores each pair of coordinates on the npts x npts grid as
# a matrix of dimension  (npts^2 x 2)
Theta <- as.matrix(expand.grid(gseq, lseq))

# log-posterior evaluation: flat prior on (gamma, lambda)
logy <- log(y) # precompute observations on log scale
lpmat <- apply(Theta, 1, function(theta) {
  weibull.loglik(theta[1], theta[2], logy = logy)
})
lpmat <- matrix(lpmat, npts, npts) # store as a grid of size npts x npts

# plot the contour along with each marginal distribution

# precompute some plot statistics
wstat <- weibull.plot.stat(gseq = gseq, lseq = lseq, lpmat = lpmat)

clrs <- c("black", "red") # plot colors
cex <- 1.5 # size of labels
# legend
lgd <- expression(group("(", list(gamma[0], lambda[0]), ")"),
                  E*group("[", list(gamma,lambda)*" | "*bold(y), "]"))
par(mfrow = c(1,3), mar = c(4,4.2,2.5,1), xpd = FALSE)
# 2D contour plot
# the power is a cheap way of changing the contour levels
contour(wstat$gseq, wstat$lseq, wstat$gldens^.3, drawlabels = FALSE,
        cex.lab = cex, cex.main = cex,
        main = expression(p(list(gamma,lambda)*" | "*bold(y))),
        xlab = expression(gamma), ylab = expression(lambda))
# add mean
points(x = c(wstat$gmean, gamma0), y = c(wstat$lmean, lambda0),
       pch = "+", cex = 4, col = clrs)
legend("topright", legend = lgd, fill = clrs, cex = 1.2)
# marginal posterior in gamma
plot(wstat$gseq, wstat$gdens, type = "l",
     cex.lab = cex, cex.main = cex,
     xlab = expression(gamma), ylab = "Density",
     main = expression(p(gamma*" | "*bold(y))))
# add mean and credible interval
abline(v = c(wstat$gmean, gamma0), col = clrs, lty = 2, lwd = 2)
# marginal posterior in lambda
plot(wstat$lseq, wstat$ldens, type = "l",
     cex.lab = cex, cex.main = cex,
     xlab = expression(lambda), ylab = "Density",
     main = expression(p(lambda*" | "*bold(y))))
# mean and credible interval
abline(v = c(wstat$lmean, lambda0), col = clrs, lty = 2, lwd = 2)

#--- mcmc sampling -------------------------------------------------------------

# random walk metropolis hastings
# tune acceptance rate to ~25%
# often this is done by trial and error with e.g., nsamples = 100
nsamples <- 1e4
rwsd <- .3*c(1,1)
theta0 <- c(1,1)
Theta.rw <- weibull.rwm(nsamples, y, theta0 = theta0,
                        rwsd = rwsd, acc.out = TRUE)
acc.rw <- Theta.rw$accept
Theta.rw <- Theta.rw$Theta

# trace plot diagnostic shows what happens when acceptance rate
# is too high or too low
rwsd.low <- c(1,1) # jumps are too big
Theta.rw.low <- weibull.rwm(nsamples, y, theta0 = theta0,
                              rwsd = rwsd.low, acc.out = TRUE)
acc.rw.low <- Theta.rw.low$accept
Theta.rw.low <- Theta.rw.low$Theta
rwsd.high <- .01*c(1,1) # jumps are too small
Theta.rw.high <- weibull.rwm(nsamples, y, theta0 = theta0,
                              rwsd = rwsd.high, acc.out = TRUE)
acc.rw.high <- Theta.rw.high$accept
Theta.rw.high <- Theta.rw.high$Theta

# diagnostics:
# 1: trace plot
# 2: acf plot
# 3: effective sample size

# trace plots
par(mfcol = c(2,3), mar = c(4,3,1,.5)+.1, oma = c(0,2,1,1),
    xpd = NA)
tnames <- c("gamma", "lambda")
clrs <- c("red", "blue", "black")
acc <- round(c(low = acc.rw.low, high = acc.rw.high, opt = acc.rw)*100)
Th <- array(cbind(Theta.rw.low, Theta.rw.high, Theta.rw),
            dim = c(nsamples, 2, 3),
            dimnames = list(iter = NULL, theta = tnames,
              acc = names(acc)))
CI <- cbind(wstat$gCI, wstat$lCI)
nkeep <- 1e3
ylim <- apply(Th[1:nkeep,,], 2, range)
ylim <- apply(rbind(ylim, CI), 2, range)
for(jj in 1:length(acc)) {
  for(ii in 1:length(tnames)) {
    mcmc.trace(x = Th[1:nkeep,ii,jj], lwd = 2,
               ylim = ylim[,ii], col = clrs[jj])
    segments(x0 = rep(par("usr")[1], 2),
             x1 = rep(par("usr")[2], 2),
             y0 = CI[,ii], lty = 2)
    if(ii == 1) {
      title(main = paste0("Acceptance Rate: ", acc[jj], "%"),
            line = .5, cex.main = 1.5)
    }
    if(jj == 1) {
      title(ylab = parse(text = paste0("p(",tnames[ii],
                           "*\" | \"*bold(y))")),
            cex.lab = 1.5)
    }
  }
}

# acf plots
par(mfrow = c(1,2), xpd = FALSE,
    mar = c(4,4,2,0)+.1, oma = c(0, 0, 0, 1))
for(ii in 1:length(tnames)) {
  for(jj in 1:length(acc)) {
    mcmc.acf(x = Th[,ii,jj], type = "o", pch = 16,
             add = jj != 1, col = clrs[jj], ylim = c(0,1))
  }
  title(main = parse(text = paste0("p(",tnames[ii],
                       "*\" | \"*bold(y))")))
  if(ii == 2) {
    legend("bottomleft", legend = paste0(acc, "%"),
           title = "Accept. Rate", fill = clrs)
  }
}

# effective sample size (ESS)
disp <- apply(Th, 2:3, effect.size)
colnames(disp) <- paste0(acc, "%")
round(disp)

# code testing
# in general, this is a very difficult problem
# (e.g., compared to checking that MLE finding algorithm has converged
# to a local mode.)
# but for a 2D problem can compare MCMC to analytic marginal posteriors

# rerun with large sample size
nsamples <- 1e5
rwsd <- .3*c(1,1)
theta0 <- c(1,1)
Theta.rw <- weibull.rwm(nsamples, y, theta0 = theta0,
                        rwsd = rwsd, acc.out = FALSE)


Theta.seq <- cbind(gamma = wstat$gseq, lambda = wstat$lseq)
Theta.dens <- cbind(gamma = wstat$gdens, lambda = wstat$ldens)

par(mfrow = c(1,2), mar = c(4,4,2,0)+.1)
tnames <- parse(text = paste0("p(", c("gamma", "lambda"), "*\" | \"*bold(y))"))
# wrap sapply in "invisible" to prevent from printing a list of NULLs
invisible(sapply(1:2, function(ii) {
  mcmc.density(x = Theta.rw[,ii], theta.name = tnames[ii],
               type = "hist", breaks = 50, freq = FALSE)
  lines(Theta.seq[,ii], Theta.dens[,ii], col = "red", lwd = 2)
}))


#--- RWM vs MWG ----------------------------------------------------------------

# metropolis-within-gibbs
## rwsd <- .25 * c(1,1) # exactly same jump size as Random Walk
rwsd <- c(.2, .5) # optimal: 45%
Theta.mwg <- weibull.mwg(nsamples, y, theta0 = theta0, rwsd = rwsd,
                         acc.out = TRUE)
acc.mwg <- Theta.mwg$accept
Theta.mwg <- Theta.mwg$Theta

# ACF plots + ESS

tnames <- c("gamma", "lambda")
clrs <- c("red", "blue")
anames <- c("RWM", "MWG")
Th <- array(cbind(Theta.rw, Theta.mwg),
            dim = c(nsamples, 2, 2),
            dimnames = list(iter = NULL, theta = tnames,
                            alg = anames))
ess <- round(apply(Th, 2:3, effect.size))

par(mfrow = c(1,2), xpd = FALSE,
    mar = c(4,4,2,0)+.1, oma = c(0, 0, 0, 1))
for(ii in 1:length(tnames)) {
  for(jj in 1:length(anames)) {
    mcmc.acf(x = Th[,ii,jj], type = "o", pch = 16,
             add = jj != 1, col = clrs[jj], ylim = c(0,1))
  }
  title(main = parse(text = paste0("p(",tnames[ii],
                       "*\" | \"*bold(y))")))
  #if(ii == 1) {
    legend("topright", legend = paste0(anames, ": ", ess[ii,]),
           title = "Effective Samp. Size", fill = clrs)
  #}
}

#--- MIID algorithm ------------------------------------------------------------

# metropolized iid draws
# mode-quadrature proposal distribution:
# theta ~ N(theta_mode, theta_quad^{-1})
prop.mq <- weibull.mq(y, theta0)

# compare contours to true posterior
glprop <- dmvnorm(as.matrix(expand.grid(gseq, lseq)),
                  mean = prop.mq$mode, sigma = solve(prop.mq$quad))
glprop <- matrix(glprop, sqrt(length(glprop)))
gldens <- wstat$gldens
par(mfrow = c(1,1), mar = c(4.1,4.1,.5,.5)+.1, oma = c(0,0,0,0))
contour(gseq, lseq, gldens^.3, drawlabels = FALSE,
        cex.lab = cex, cex.main = cex, lwd = 2,
        #main = expression(p(list(gamma,lambda)*" | "*bold(y))),
        xlab = expression(gamma), ylab = expression(lambda))
contour(gseq, lseq, glprop^.3, add = TRUE, lwd = 2,
        col = "red", drawlabels = FALSE)
legend("topright", legend = c("True Posterior", "Mode-Quad Approx"),
       fill = c("black", "red"), cex = 1)

# metropolized iid draws -- no tuning parameters
nsamples <- 1e5
Theta.mi <- weibull.miid(nsamples, y, theta0 = theta0)

# compare effective sample sizes for all samplers
disp <- sapply(list(rwm = Theta.rw, mwg = Theta.mwg, miid = Theta.mi),
               function(Th) apply(Th, 2, effect.size))
round(disp)

#--- Marginal MCMC -------------------------------------------------------------

# let eta = 1/lambda^gamma
# then for pi(gamma, eta) = pi(gamma) * eta^(alpha-1) * exp(-beta * eta)
# posterior p(gamma | y) available analytically

# *** IMPORTANT *** change-of-variables:

# note that pi(gamma, eta) \propto 1 <=>
# pi(gamma, lambda) \propto gamma/lambda^(gamma+1)

# similarly, p(gamma, lambda | y) = p(gamma, eta | y) * gamma/lambda^(gamma+1)

# compare posterior for both priors:
# 1. pi(gamma, lambda) \propto 1
# 2. pi(gamma, eta) \propto 1
npts <- 100
gseq <- seq(.75, 1.75, len = npts)
lseq <- seq(1.25, 3.5, len = npts)
glseq <- as.matrix(expand.grid(gseq, lseq))
logy <- log(y)

# pi(gamma, lambda) \propto 1
lpmat1 <- apply(glseq, 1, function(theta) {
  weibull.loglik(gamma = theta[1], lambda = theta[2], logy = logy)
})
lpmat1 <- matrix(lpmat1, npts, npts)

# pi(gamma, eta) \propto 1 <=>
# pi(gamma, lambda) \propto gamma/lambda^(gamma+1)
lpmat2 <- apply(glseq, 1, function(theta) {
  # transform to (gamma, eta)
  gamma <- theta[1]
  lambda <- theta[2]
  eta <- 1/lambda^gamma
   # precompute quantity used both in marginal and conditional
  Tgamma <- sum(exp(gamma * logy))
  # p(gamma | y)
  lp <- weibull.lpmarg(gamma, logy,
                       alpha0 = 1, beta0 = 0, Tgamma = Tgamma)
  # p(eta | gamma, y)
  lp <- lp + weibull.lpcond(eta, gamma, logy,
                            alpha0 = 1, beta0 = 0, Tgamma = Tgamma)
  # change-of-variables to get p(gamma, lambda)
  lp + log(gamma) - (gamma+1) * log(lambda)
})
lpmat2 <- matrix(lpmat2, npts, npts)

wstat1 <- weibull.plot.stat(gseq, lseq, lpmat1)
wstat2 <- weibull.plot.stat(gseq, lseq, lpmat2)

# plot
legend.names <- expression(pi(gamma,lambda) %prop% 1,
                           pi(gamma,eta) %prop% 1,
                           "True Parameter Value")
par(mfrow = c(1,3), mar = c(4,4.2,2.5,1))
clrs <- c("black", "blue", "red")
cex <- 1.5
# contours
contour(wstat1$gseq, wstat1$lseq, wstat1$gldens^.3, drawlabels = FALSE,
        cex.lab = cex, cex.main = cex, col = clrs[1],
        main = expression(p(list(gamma,lambda)*" | "*bold(y))),
        xlab = expression(gamma), ylab = expression(lambda))
contour(wstat2$gseq, wstat2$lseq, wstat2$gldens^.3, drawlabels = FALSE,
        col = clrs[2], add = TRUE)
points(x = c(wstat1$gmean, wstat2$gmean, gamma0),
       y = c(wstat1$lmean, wstat2$lmean, lambda0),
       pch = "+", cex = 4, col = clrs)
legend("topright", legend = legend.names,
       fill = clrs, cex = 1.2)
# p(gamma | y)
plot(x = 0, type = "n",
     xlim = range(wstat1$gseq, wstat2$gseq),
     ylim = range(wstat1$gdens, wstat2$gdens),
     cex.lab = cex, cex.main = cex,
     xlab = expression(gamma), ylab = "Density",
     main = expression(p(gamma*" | "*bold(y))))
lines(wstat1$gseq, wstat1$gdens, col = clrs[1])
lines(wstat2$gseq, wstat2$gdens, col = clrs[2])
abline(v = c(wstat1$gmean, wstat2$gmean, gamma0),
       col = clrs, lty = 2, lwd = 2)
# p(lambda | y)
plot(x = 0, type = "n",
     xlim = range(wstat1$lseq, wstat2$lseq),
     ylim = range(wstat1$ldens, wstat2$ldens),
     cex.lab = cex, cex.main = cex,
     xlab = expression(lambda), ylab = "Density",
     main = expression(p(lambda*" | "*bold(y))))
lines(wstat1$lseq, wstat1$ldens, col = clrs[1])
lines(wstat2$lseq, wstat2$ldens, col = clrs[2])
abline(v = c(wstat1$lmean, wstat2$lmean, lambda0),
       col = clrs, lty = 2, lwd = 2)

# 1D Marginal MIID
Theta.mmi <- weibull.mmi(nsamples, y, gamma.rng = c(0, 10))

# compare to 2-d MIID
Theta.mi <- weibull.miid(nsamples, y, theta0 = c(1,1))

# effective sample size
sapply(list(MI_2D = Theta.mi, MI_1D = Theta.mmi),
       function(Theta) apply(Theta, 2, effect.size))

# check that posteriors are correct
par(mfrow = c(1,2))
hist(Theta.mmi[,"gamma"], breaks = 100, freq= FALSE,
     main = expression(p(gamma*" | "*bold(y))),
     xlab = expression(gamma))
lines(wstat2$gseq, wstat2$gdens, col = "red")
hist(Theta.mmi[,"lambda"], breaks = 100, freq= FALSE,
     main = expression(p(lambda*" | "*bold(y))),
     xlab = expression(lambda))
lines(wstat2$lseq, wstat2$ldens, col = "red")
