#--- comparison of some bootstrap methods ---------------------------------
#
# Inference Model: Ui ~iid Unif(0, theta)
#
# Data-Generating Process: Ui ~iid theta * Beta(alpha, alpha)
#   Note that Beta(1,1) = Unif(0,1)
#
# Estimand: theta
#
# Estimators: theta1.hat = max(U), theta2.hat = 2*mean(U)
#
# Sample size: n = 100, n = 1000
#
# Bootstrap techniques:
# * Non-parametric, i.e., resample U.boot from U without replacement
# * Parametric, i.e., U.boot ~iid Unif(0, theta.hat)
#
# Confidence intervals:
# * Basic Bootstrap:
#     CI = theta.hat + quantiles(theta.hat - theta.boot)
# * Percentile Bootstrap:
#     CI = quantiles(theta.boot)
#
#--------------------------------------------------------------------------

#' Generate a dataset, estimates, and bootstrap confidence intervals.
#'
#' @param theta0 true parameter value
#' @param alpha parameter specifying the data-generating process, which is Beta(alpha, alpha).
#' @param n sample size
#' @param nboot number of bootstrap samples to use for calculating CI's.
#' @return A \code{2 x 2 x 2 x 2} array with the following dimensions:
#' \describe{
#'   \item{\code{est}}{Parameter estimator: (i) sample max (max), (ii) 2x sample mean (mean2)}
#'   \item{\code{samp}}{Bootstrap sampling method: (i) nonparametric (NP), (ii) parametric (P)}
#'   \item{\code{int}}{Bootstrap CI method: (i) basic (basic), (ii) percentile (pct)}
#'   \item{\code{CI}}{Limits of confidence interval: lower (L = 2.5%), upper (U = 97.5%)}
#' }
#' @details \itemize{
#'   \item Simulations are "vectorized" for both bootstrap methods, i.e., all (n x nboot) bootstrap samples are generated in a single call.  This is much faster than drawing each of the nboot bootstrap samples separately.
#'   \item We have two estimators theta1.hat and theta2.hat, for which we can use the same non-parametric bootstrap samples.  For parametric bootstrap, in general we would have to generate separate samples for theta1.hat and theta2.hat.  However, in this case we can simple draw U.boot ~iid Unif(0,1), and multiply by theta.hat to get U.boot ~iid Unif(0, theta.hat).  Not only does this cut down on the computation time, but it also reduces the "Monte Carlo error" in comparing CI's for theta1.hat vs theta2.hat.
#' }
boot.sim <- function(theta0 = 1, alpha = 1, n, nboot) {
  # generate data
  U <- theta0 * rbeta(n = n, shape1 = alpha, shape2 = alpha)
  # estimates
  theta1.hat <- max(U) # estimator 1
  theta2.hat <- 2*mean(U) # estimator 2
  # bootstrap sampling
  # non-parametric
  U.bootNP <- sample(U, size = n*nboot, replace = TRUE)
  U.bootNP <- matrix(U.bootNP, n, nboot)
  # parametric
  # multiply by theta.hat to get the right distribution
  U.bootP <- matrix(runif(n*nboot), n, nboot)
  # estimators
  Theta.boot <- array(NA, dim = c(2,2,nboot),
                      dimnames = list(c("max", "mean2"),
                        c("NP", "P"), NULL))
  Theta.boot[1,1,] <- apply(U.bootNP, 2, max) # max/np
  Theta.boot[1,2,] <- apply(U.bootP * theta1.hat, 2, max) # max/p
  Theta.boot[2,1,] <- 2*colMeans(U.bootNP) # mean2/np
  Theta.boot[2,2,] <- 2*colMeans(U.bootP * theta2.hat) # mean2/p
  # confidence intervals
  CI <- array(NA, dim = c(2,2,2,2))
  dimnames(CI) <- list(est = c("max", "mean2"), samp = c("NP", "P"),
                       CI = c("L", "U"), int = c("basic", "pct"))
  # basic bootstrap
  Theta.hat <- c(max = theta1.hat, mean2 = theta2.hat)
  T.boot <- Theta.hat - Theta.boot # emulate T = theta0 - theta.hat
  T.CI <- apply(T.boot, 1:2,
                quantile, probs = c(.025, .975)) # quantiles of T.boot

  T.CI <- aperm(T.CI,                   # permute dimensions
                perm = c(2,3,1))        # quantiles are now last dim
  CI[,,,"basic"] <- T.CI + Theta.hat
  # percentile bootstrap
  CI.pct <- apply(Theta.boot, 1:2, quantile, probs = c(.025, .975))
  CI[,,,"pct"] <- aperm(CI.pct, perm = c(2,3,1))
  aperm(CI, perm = c(1,2,4,3))
}


#--- simulation setup -----------------------------------------------------

require(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl = cl)
clusterSetRNGStream(cl) # need parallel pseudo-rng's

theta0 <- 1 # true parameter value
# note that the problem doesn't really depend on theta0 since
# we can always rescale things such that theta0 = 1
N <- c(100, 1e3, 1e4) # sample size
Alpha <- c(1, 2) # DGP: U ~ theta0 * Beta(alpha, alpha)
nboot <- 1e3 # number of bootstrap samples
# number of times to repeat the whole experiment
# i.e., generate data and calculate various CI's
nsim <- 1e3

# the following method is useful for going through
# contingency tables with a single for-loop
sim.par <- expand.grid(n = N, alpha = Alpha) # simulation parameters
npar <- nrow(sim.par)
CI.Boot <- vector("list", npar) # list of outputs

system.time({
  for(ii in 1:npar) {
    message("n = ", sim.par$n[ii], ", alpha = ", sim.par$alpha[ii])
    # each element of ci.boot contains nsim bootstrap intervals
    # of each type
    ci.boot <- foreach(jj=1:nsim) %dopar% {
      boot.sim(theta0 = theta0, nboot = nboot,
               alpha = sim.par$alpha[ii],
               n = sim.par$n[ii])
    }
    # format list output as array and store
    CI.Boot[[ii]] <- array(unlist(ci.boot),
                           dim = c(dim(ci.boot[[1]]), nsim),
                           dimnames = c(dimnames(ci.boot[[1]]),
                             list(NULL)))
  }
})

stopCluster(cl) # deallocate resources when you're done

# store simulations as contingency table
dim(CI.Boot) <- c(length(N),length(Alpha))
dimnames(CI.Boot) <- list(n = N, alpha = Alpha)
CI.Boot

# always a good idea to save time-consumming simulations
# also, i like to comment out "save" commands to avoid over-writing
# things by accidentally pressing shortcut keys
# (this only has to happen to you once
# for you to never want it to happen again)
# save(CI.Boot, file = "unif-boot_sim5.RData")

#--- coverage and interval width calculations -----------------------------

load("unif-boot_sim5.RData")

# true coverage
CI.cover <- sapply(CI.Boot, function(ci.boot) {
  apply(ci.boot, 1:3, function(ci) {
    mean(ci["L",] <= theta0 & theta0 <= ci["U",])
  })
})

# save as "flat" table
CI.cover <- cbind(expand.grid(stat = c("max", "mean2"),
                              sim = c("NP", "P"),
                              int = c("basic", "pct"),
                              n = N, alpha = Alpha,
                              stringsAsFactors = FALSE),
                  cover = c(CI.cover))

head(CI.cover)

# rearrange/merge levels for simple display
CI.disp <- function(CI, responseName) {
  frm <- paste0(responseName, " ~ .")
  disp <- xtabs(formula = as.formula(frm),
                data = CI) # convert to regular table
  disp <- aperm(disp, c(4, 1, 2, 3, 5)) # permute dimensions for display
  # back to flat table
  disp <- as.data.frame.table(disp, responseName = responseName,
                              stringAsFactors = FALSE)
  # merge a few factors
  disp$sim.stat <- paste0(disp$sim, "_", disp$stat)
  disp$int.n <- paste0(disp$int, "_n=", disp$n)
  # back to contingency table
  frm <- paste0(responseName, " ~ int.n + sim.stat + alpha")
  disp <- xtabs(as.formula(frm), data = disp)
  names(dimnames(disp))[1:2] <- ""
  disp
}

round(CI.disp(CI.cover, "cover"), 2)

# mean interval width
CI.width <- sapply(CI.Boot, function(ci.boot) {
  apply(ci.boot, 1:3, function(ci) {
    mean(ci["U",] - ci["L",])
  })
})
CI.width <- cbind(expand.grid(stat = c("max", "mean2"),
                              sim = c("NP", "P"),
                              int = c("basic", "pct"),
                              n = N, alpha = Alpha),
                  width = c(CI.width))

round(CI.disp(CI.width, "width"), 2)


#--- scratch --------------------------------------------------------------


# true coverage
CI.cover <- apply(CI.Boot[[1,1]], c(1,2,4), function(ci) {
  mean(ci["2.5%",] <= theta0 & theta0 <= ci["97.5%",])
})

# mean interval width
CI.width <- apply(CI.Boot, c(1,2,4), function(ci) {
  mean(ci["97.5%",] - ci["2.5%",])
})



# display latex table
ltable1 <- function(CI.basic, CI.pct, wd1 = "1cm", wd2 = "7em",
                    digits = 3) {
  tmp <- paste0("K{", rep(wd1, 2), "}", collapse = "")
  message("\\begin{tabular}{r", tmp, "c", tmp, "}")
  tmp <- paste0("\\multicolumn{2}{K{",wd2,"}}")
  message("& ", tmp, "{Non-Parametric} & & ", tmp, "{Parametric} \\\\")
  tmp <- "& Max & Mean"
  message(tmp, " & ", tmp, "\\\\ \\cline{2-3}\\cline{5-6}")
  ci <- c(apply(signif(CI.basic, digits), 1:2,
                function(ci) paste0("(",ci[1],", ",ci[2],")")))
  ci <- paste0(c(ci[1:2], "", ci[3:4]), collapse = " & ")
  message("Basic & ", ci, " \\\\")
  ci <- c(apply(signif(CI.pct, digits), 1:2,
                function(ci) paste0("(",ci[1],", ",ci[2],")")))
  ci <- paste0(c(ci[1:2], "", ci[3:4]), collapse = " & ")
  message("Percentile & ", ci, " \\\\")
  message("\\end{tabular}")
}

ltable1(CI.basic = CI.basic, CI.pct = CI.pct,
        wd1 = "6em", wd2 = "12em", digits = 3)
