#--- euler vs true transition density for geometric brownian motion ------------

source("sde-functions.R")

# simulation parameters

S0 <- 5.1 # initial value
alpha <- .1 # rate parameter
sigma <- .3 # volatility parameter
T <- 3.2 # total time

Nseq <- c(1, 2, 5, 20, 100) # number of Euler steps, i.e., dt = T/N
NN <- length(Nseq)
nsim <- 1e5 # number of time series to generate

St.sim <- matrix(NA, nsim, NN) # store output
for(ii in 1:NN) {
  N <- Nseq[ii]
  dT <- T/N
  # keep last timepoint only, i.e., N*dT = T
  St.sim[,ii] <- gbm.sim(nsim, alpha = alpha, sigma = sigma, dT = dT,
                         s0 = S0, N = N)[N,]
}

#--- plot smoothed histograms and compare to true transition density. ----------

# number of points in kernel density.
# for numerical reasons best to make this power of 2.
ndens <- 512
xdens <- matrix(NA, ndens, NN)
ydens <- xdens
for(ii in 1:NN) {
  dens <- density(St.sim[,ii], n = ndens)
  xdens[,ii] <- dens$x
  ydens[,ii] <- dens$y
}
# calculate true pdf
xtrue <- seq(min(xdens), max(xdens), len = ndens)
ytrue <- dgbm(s = xtrue, alpha = alpha, sigma = sigma,
              s0 = S0, t = T)


par(mfrow = c(1,1), mar = c(4,4,1,.5)+.1)
clrs <- c("red", "blue", "orange", "brown", "green3")
lwd <- 2
plot(0, type = "n", xlim = c(min(xdens, xtrue), 30),
     ylim = range(ydens, ytrue),
     xlab = expression(S[T]), ylab = expression("Density"))
for(ii in 1:NN) {
  # euler approximations
  lines(xdens[,ii], ydens[,ii], col = clrs[ii], lwd = lwd)
}
lines(xtrue, ytrue, col = "black", lwd = lwd) # true density
# legend
lgd <- paste0("p(tilde(S)[T]^{(", Nseq, ")}*\" | \"*S[0])")
lgd <- c("p(S[T]*\" | \"*S[0])", lgd)
lgd <- parse(text = lgd)
legend("topright", legend = parse(text = lgd),
       fill = c("black", clrs))
