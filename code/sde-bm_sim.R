#--- simulate discrete-time skeleton of Brownian motion ------------------------

source("sde-functions.R")


# initial skeleton

T <- 2.5 # endpoint
nsteps <- 6 # dt = T/nsteps
nsim <- 5 # number of simulations

# simulate data
tseq <- seq(0, T, len = nsteps + 1) # simulation timepoints
# replicate repeats a calculation "n" times and concatenates the results
Bt <- replicate(n = nsim, expr = {
  # initial value is always 0
  c(0, rbm(tseq[-1]))
})

# plot trajectories
clrs <- rainbow(nsim+1)[1:nsim] # colors to use
# trick for plotting math symbols with variable inputs
main <- paste0("Delta*t==", signif(tseq[2] - tseq[1],2))
main <- parse(text = main)
par(mfrow = c(1,1), mar = c(4.5, 4.5, 2, .5)+.1)
# create an empty plot big enough to contain all the curves
plot(0, type = "n", xlim = range(tseq), ylim = range(Bt),
     xlab = expression(t), ylab = expression(B[t]),
     main = main)
for(ii in 1:nsim) {
  lines(tseq, Bt[,ii], lty = 3, pch = 21, type = "o", bg = clrs[ii])
}

#--- infill current skeleton ---------------------------------------------------

# given a time sequence and Bm skeletons,
# infill each by halving time interval between each observation
bm.infill <- function(tseq, Bt) {
  N <- length(tseq) - 1 # number of points to infill
  n <- ncol(Bt) # number of Bm time series
  # lower endpoints
  tL <- tseq[1:N]
  BL <- Bt[1:N,]
  # upper endpoints
  tU <- tseq[1+1:N]
  BU <- Bt[1+1:N,]
  # midpoints
  tM <- (tL + tU)/2
  # infill skeleton
  Bt2 <- matrix(NA, 2*N+1, n) # storage for new skeleton
  # T/F vector corresponding to the infill timepoints
  ifill <- rep(FALSE, 2*N+1)
  ifill[seq(2, 2*N+1, by = 2)] <- TRUE
  for(ii in 1:n) {
    # generate brownian bridge for each skeleton
    Bt2[ifill,ii] <- rbb(tseq = tM, tL = tL, tU = tU,
                         BL = BL[,ii], BU = BU[,ii])
  }
  Bt2[!ifill,] <- Bt # previous skeleton
  tseq2 <- rep(NA, 2*N+1) # storage for new timepoints
  tseq2[ifill] <- tM # new timepoints
  tseq2[!ifill] <- tseq # previous timepoints
  list(tseq = tseq2, Bt = Bt2)
}

Bt2 <- bm.infill(tseq = tseq, Bt = Bt)

clrs <- rainbow(nsim+1)[1:nsim] # colors to use
main <- paste0("Delta*t==", signif(Bt2$tseq[2] - Bt2$tseq[1],2))
main <- parse(text = main)
par(mfrow = c(1,1), mar = c(4.5, 4.5, 2, .5)+.1)
# create an empty plot big enough to contain all the curves
plot(0, type = "n", xlim = range(Bt2$tseq), ylim = range(Bt2$Bt),
     xlab = expression(t), ylab = expression(B[t]),
     main = main)
for(ii in 1:nsim) {
  lines(Bt2$tseq, Bt2$Bt[,ii], lty = 3, pch = 21, type = "o", bg = clrs[ii])
}

#--- infill many times to get a very fine skeleton -----------------------------

# cut interval in half nfill times
# so total number of timepoints is nsteps * 2^nfill + 1
nfill <- 10

BtN <- list(tseq = tseq, Bt = Bt) # initial value
for(ii in 1:nfill) {
  BtN <- bm.infill(tseq = BtN$tseq, Bt = BtN$Bt)
}

clrs <- rainbow(nsim+1)[1:nsim]
main <- paste0("Delta*t==", signif(BtN$tseq[2] - BtN$tseq[1],2))
main <- parse(text = main)
par(mfrow = c(1,1), mar = c(4.5, 4.5, 2, .5)+.1)
plot(0, type = "n", xlim = range(BtN$tseq), ylim = range(BtN$Bt),
     xlab = expression(t), ylab = expression(B[t]),
     main = main)
for(ii in 1:nsim) {
  lines(BtN$tseq, BtN$Bt[,ii])
  lines(BtN$tseq, BtN$Bt[,ii], lty = 3, col = clrs[ii])
}

