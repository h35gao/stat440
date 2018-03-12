#--- hierarchical modeling of the 8 schools data --------------------------

source("hiernorm-functions.R")

schools <- data.frame(school = LETTERS[1:8],
                      estimate = c(28, 8, -3, 7, -1, 1, 18, 12),
                      sd = c(15, 10, 16, 11, 9, 11, 10, 18))
schools

# hierarchical model is:
# Xi | mui ~ind N(mui, sigi^2)
# mui ~iid N(lambda, tau^2)

schX <- schools$estimate
schV <- schools$sd^2
K <- length(schX)

# plot Likelihood
tau.seq <- seq(0, 20, len = 1e3)
tau.ll <- sapply(tau.seq, lprof.tau, x = schX, sig2 = schV)
par(mfrow = c(1,1), mar = c(3.5, 3.8, .5, .5))
plot(tau.seq, exp(tau.ll-max(tau.ll)), type = "l",
     xlab = "", ylab = "")
title(xlab = expression(tau),
      ylab = expression(L[prof](tau*" | "*bold(x))),
      line = 2.25)

# on certain graphical devices (e.g., quartz on osx) you can plot the "ell"
# symbol for log-likelihood as follows
plot(tau.seq, tau.ll, type = "l", xlab = "", ylab = "")
title(xlab = expression(tau),
      ylab = expression("\u2113"[prof](tau*" | "*bold(x))),
      line = 2.25)

#--- bootstrap ------------------------------------------------------------

# maximum likelihood of lambda
tau.mle <- 0
lambda.mle <- lambda.prof(tau = tau.mle, x = schX, sig2 = schV)

# bootstrap estimates of theta = c(lambda, tau)
nboot <- 1e4
theta.bootP <- matrix(NA, nboot, 2) # parametric
colnames(theta.bootP) <- c("lambda", "tau")
theta.bootNP <- theta.bootP # non-parametric

for(ii in 1:nboot) {
  if(ii %% 5000 == 0) message("ii = ", ii)
  # parametric bootstrap
  # simulate responses
  x.boot <- rnorm(K, mean = lambda.mle, sd = sqrt(schV + tau.mle^2))
  sig2.boot <- schV # variances remain the same
  # mle
  theta.bootP[ii,] <- theta.mle(x = x.boot, sig2 = sig2.boot)
  # nonparametric bootstrap
  # resample schools (both means and variances)
  ind <- sample(K, replace = TRUE)
  x.boot <- schX[ind]
  sig2.boot <- schV[ind]
  theta.bootNP[ii,] <- theta.mle(x = x.boot, sig2 = sig2.boot)
}

# fraction of tau_boot == 0
# or rather < 1e-4 for numerical tolerance
tau.hat0 <- colMeans(cbind(P = theta.bootP[,"tau"],
                           NP = theta.bootNP[,"tau"]) < 1e-4)


# plot
par(mfrow = c(1, 2), mar = c(4,4,2,1))
hist(theta.bootP[,"tau"], breaks = 100, freq= FALSE,
     main = expression("Parametric Bootstrap"),
     xlab = expression(tilde(tau)))
legend("topright",
       legend = parse(text = paste0("Pr(tilde(tau)==0)==",
                                    signif(tau.hat0["P"]*100, 2),
                                    "*\"%\"")))
hist(theta.bootNP[,"tau"], breaks = 100, freq= FALSE,
     main = expression("Non-Parametric Bootstrap"),
     xlab = expression(tilde(tau)))
legend("topright",
       legend = parse(text = paste0("Pr(tilde(tau)==0)==",
                                    signif(tau.hat0["NP"]*100, 2),
                                    "*\"%\"")))

# bootstrap distribution vs fully pooled estimator

nsamples <- 1e4
mu.bootP <- rcond.mu(n = nsamples, tau = theta.bootP[,"tau"],
                     lambda = theta.bootP[,"lambda"],
                     x = schX, sig2 = schV)

mu.bootNP <- rcond.mu(n = nsamples, tau = theta.bootNP[,"tau"],
                      lambda = theta.bootP[,"lambda"],
                      x = schX, sig2 = schV)

lambda.se <- sqrt(1/sum(1/(schV + tau.mle)))

ylab <- expression("Density")
xlab <- expression(tilde(mu)[i])
clrs <- rainbow(K+1)[1:K]
par(mfrow = c(1,2), mar = c(4,4,2,.5), xpd = FALSE)
# parametric bootstrap
main <- "Parametric Bootstrap"
multi.dens(t(mu.bootP), col = clrs, lwd = 2, xlim = c(-15, 30))
curve(dnorm(x, mean = lambda.mle, sd = lambda.se), add = TRUE,
      lty = 2, lwd = 2)
title(xlab = xlab, ylab = ylab, main = main)
legend("topleft", legend = LETTERS[1:K], fill = clrs, title = "School",
       ncol = 2)
legend("topright", legend = "Full Pooling",
       lty = 2, lwd = 2)
# nonparametric bootstrap
main <- "Nonparametric Bootstrap"
multi.dens(t(mu.bootNP), col = clrs, lwd = 2, xlim = c(-15, 30))
title(xlab = xlab, ylab = ylab, main = main)
curve(dnorm(x, mean = lambda.mle, sd = lambda.se), add = TRUE,
      lty = 2, lwd = 2)


#--- bayesian inference ---------------------------------------------------

# compare posteriors for different priors:

tau.seq <- seq(0, 40, len = 1e3)
# pi(tau) \propto 1
tau.lp1 <- sapply(tau.seq, lpmarg.tau, x = schX, sig2 = schV)
# pi(tau^2) \propto 1
tau.lp2 <- tau.lp1 + log(tau.seq)

# numerical calculation of posterior mean
tau.pmean <- function(tau.seq, tau.lp) {
  tau.post <- exp(tau.lp - max(tau.lp)) # unnormalized posterior
  sum(tau.seq * tau.post)/sum(tau.post)
}

clrs <- c("black", "red")
par(mfrow = c(1,1), mar = c(4, 4.5, .5, .5))
plot(tau.seq, exp(tau.lp1 - max(tau.lp1)),
     type = "l", col = clrs[1], lwd = 2,
     xlab = expression(tau),
     ylab = expression(p(tau*" | "*bold(x))))
lines(tau.seq, exp(tau.lp2 - max(tau.lp2)), col = clrs[2], lwd = 2)
abline(v = c(tau.pmean(tau.seq, tau.lp1), tau.pmean(tau.seq, tau.lp2)),
       col = clrs, lwd = 2, lty = 2)
legend("topright",
       legend = expression(pi(tau)%prop%1, pi(tau^2)%prop%1),
       fill = clrs)

#--- posterior inference for mu ------------------------------------------------

nsamples <- 1e4
tau.seq <- seq(0, 50, len = 1e5)
# pi(tau) \propto 1
tau.lp1 <- sapply(tau.seq, lpmarg.tau, x = schX, sig2 = schV)
# p(tau | x)
tau.post <- sample(x = tau.seq, size = nsamples,
                   replace = TRUE, prob = exp(tau.lp1 - max(tau.lp1)))
# p(lambda | tau, x)
system.time({
  lambda.post <- rcond.lambda(n = nsamples, tau = tau.post,
                              x = schX, sig2 = schV)
})
# p(mu | x)
mu.post <- rcond.mu(n = nsamples, tau = tau.post, lambda = lambda.post,
                    x = schX, sig2 = schV)


# compare p(mu | x) and parametric bootstrap
clrs <- rainbow(K+1)[1:K]
par(mfrow = c(1,2), mar = c(4,4.2,2,.3)+.1)
# parametric bootstrap
ylab <- expression(p[Boot](tilde(mu)[i]))
xlab <- expression(tilde(mu)[i])
## main <- parse(text = paste0("\"Parametric Bootstrap, \"*hat(tau)[ML]==",
##                             tau.mle))
main <- "Parametric Bootstrap"
multi.dens(t(mu.bootP), col = clrs, lwd = 2)
curve(dnorm(x, mean = lambda.mle, sd = lambda.se), add = TRUE,
      lty = 2, lwd = 2)
title(xlab = xlab, ylab = ylab, main = main)
legend("topleft", legend = LETTERS[1:K],
       fill = clrs, title = "School", ncol = 2)
legend("topright", legend = "Full Pooling",
       lty = 2, lwd = 2)
# bayesian
ylab <- expression(p(mu[i]*" | "*bold(x)))
xlab <- expression(mu[i])
## main <- parse(text = "\"Bayesian Posterior\"")
main <- "Bayesian Posterior"
multi.dens(t(mu.post), col = clrs, lwd = 2)
curve(dnorm(x, mean = lambda.mle, sd = lambda.se), add = TRUE,
      lty = 2, lwd = 2)
title(xlab = xlab, ylab = ylab, main = main)
#legend("topright", legend = LETTERS[1:K], fill = clrs, title = "School")


#--- now the actual question of interest, i.e., rankings -----------------------

# rank of each school per posterior sample
rank.post <- apply(-mu.post, 2, rank) # minus sign since highest score
rank.bootP <- apply(-mu.bootP, 2, rank) # is first ranked

clrs <- rainbow(9)[-9]
par(mfrow = c(2,4), mar = c(4, 4, .5, .5), oma = c(0,0,2,0))
for(ii in 1:K) {
  barplot(height = table(rank.post[ii,])/nsamples, col = clrs[ii],
          xlab = paste("School", LETTERS[ii]), names.arg = 1:K)
}
mtext(side = 3, line = .5, outer = TRUE,
      text = "Posterior Rank Distribution per School")

clrs <- rainbow(9)[-9]
par(mfrow = c(2,4), mar = c(4, 4, .5, .5), oma = c(0,0,2,0))
for(ii in 1:K) {
  barplot(height = table(rank.bootP[ii,])/nsamples, col = clrs[ii],
          xlab = paste("School", LETTERS[ii]), names.arg = 1:K)
}
mtext(side = 3, line = .5, outer = TRUE,
      text = "Bootstrap Rank Distribution per School")


# order of the schools within each rank
order.post <- apply(mu.post, 2, order, decreasing = TRUE)
order.bootP <- apply(mu.bootP, 2, order, decreasing = TRUE)

# bayesian
clrs <- rainbow(9)[-9]
par(mfrow = c(2,4), mar = c(4, 4, .5, .5), oma = c(0,0,2,0))
for(ii in 1:K) {
  barplot(height = table(order.post[ii,])/nsamples, col = clrs,
          xlab = paste("Rank", ii), names.arg = LETTERS[1:8])
}
mtext(side = 3, line = .5, outer = TRUE,
      text = "Posterior School Distribution per Rank")

# parametric bootstrap
clrs <- rainbow(9)[-9]
par(mfrow = c(2,4), mar = c(4, 4, .5, .5), oma = c(0,0,2,0))
for(ii in 1:K) {
  barplot(height = table(order.bootP[ii,])/nsamples, col = clrs,
          xlab = paste("Rank", ii), names.arg = LETTERS[1:8])
}
mtext(side = 3, line = .5, outer = TRUE,
      text = "Bootstrap School Distribution per Rank")

