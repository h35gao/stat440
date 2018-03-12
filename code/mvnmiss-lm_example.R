# missing data in multivariate normal models

source("mvnmiss-functions.R")

#--- missing data example -------------------------------------------------

# regression model: y = alpha x + beta z + sigma eps

# extreme example
alpha <- 1
beta <- 1
sigma <- 1
# NOTE: don't use "nice" numbers for debugging code, because sometimes
# with nice numbers two errors can cancel each other out.
## alpha <- runif(1, .9, 1.2)
## beta <- runif(1, .5, 1.7)
## sigma <- runif(1, 1, 2)
theta0 <- c(alpha = alpha, beta = beta, sigma = sigma)

n <- 1e6
x <- rnorm(n)
z <- rnorm(n)
y <- alpha*x + beta*z + sigma*rnorm(n)

# missingness mechanism: P(missing) depends on y
C <- 2
p <- 1-c(.05, .9)
delta <- rbinom(n = n, size = 1, prob = ifelse(y < C, p[1], p[2]))

# regression comparison
M0 <- lm(y ~ x + z - 1) # without missing data
M1 <- lm(y ~ x + z - 1, subset = delta == 1) # with complete cases only

# beta's
coef(summary(M0))
coef(summary(M1)) # both estimate and 95% CI way off
# sigma's
c(M0 = sigma(M0), M1 = sigma(M1)) # M1 estimate way off

#--- missing data model for multivariate normals ---------------------------

# the model is:
#
# y = (x, z) ~ N(0, V)
#
# where x = (x1, ..., xp), z = (z1, ..., zq).
#
# for n iid observations, let S1 denote the complete cases, i.e., those
# for which both x and z are observed.
# let S0 denote the observations for which only x is observed.
#
# in other words, the observed data is Y1, X0, i.e.,
# the matrix of complete cases, and the matrix of incomplete cases for
# which only x is observed.

# Note: this model can be used for estimating the parameters of
# the regression problem by letting
# x_mvn <- c(y_reg, x_reg) and z_mvn <- z_reg.

# test code for the EM algorithm for fitting the mvn-miss model

source("stat440-utils.R")

# generate data
px <- 3
pz <- 2
n1 <- 100
n0 <- 100
py <- px+pz
V <- crossprod(matrix(rnorm(py^2),py))

Y1 <- rmvnorm(n = n1, sigma = V)
X0 <- rmvnorm(n = n0, sigma = V[1:px,1:px,drop=FALSE])

V.mle <- mvn0.fit(Y1 = Y1, X0 = X0, tol = 1e-8, niter = 1e3)

# take lower triangular part only
loglik <- function(L) {
  V <- matrix(NA,py,py)
  V[lower.tri(V,diag=TRUE)] <- L
  V <- t(V)
  V[lower.tri(V,diag=TRUE)] <- L
  if(any(eigen(V)$val <= 0)) return(-Inf)
  mvn0.loglik(V, Y1, X0)
}

tnames <- apply(expand.grid(1:py,1:py), 1, paste0, collapse = "")
tnames <- matrix(paste0("V[", tnames, "]"),py)
tnames <- parse(text = tnames[lower.tri(tnames,diag=TRUE)])
mle.check(loglik = loglik, theta.names = tnames,
          theta.mle = V.mle[lower.tri(V.mle,diag=TRUE)])

#--- application to regression model ---------------------------------------

# inference model: (y,z,x) ~ N(0, V)
# true model: p(x) ~ student-t
#             (y,z) | x ~ normal

# since the regression parameters don't appear in the marginal distribution
# of x, the inference model for x doesn't need to be correctly specified

# generate data
n <- 1e4 # use smaller sample size
x <- 1 + rt(n, df = 8) # not normal
# (y,z) | x ~ normal
z <- rnorm(n) # normal
y <- alpha*x + beta*z + sigma*rnorm(n) # normal
# missingness mechanism
C <- 2
p <- 1-c(.05, .9)
delta <- rbinom(n = n, size = 1, prob = ifelse(y < C, p[1], p[2]))

# observed data
D <- cbind(y = y, x = x, z = z, delta = delta)
D[D[,"delta"] == 0,"z"] <- NA # delete z when it's missing

# prepare arguments for EM calculation
Y1 <- D[D[,"delta"] == 1,c("x","y","z")]
X0 <- D[D[,"delta"] == 0, c("x","y")]

system.time({
  V.mle <- mvn0.fit(Y1 = Y1, X0 = X0)
})

# extract theta = (alpha, beta, sigma) from V
V2theta <- function(V) {
  theta <- solve(V[c("x", "z"),c("x","z")], V["y",c("x","z")])
  theta <- c(theta,
             sqrt(V["y","y"] - crossprod(V["y",c("x","z")], theta)))
  names(theta) <- c("alpha", "beta", "sigma")
  theta
}

theta.mle <- V2theta(V.mle)

# extremely close
rbind(true = c(alpha = alpha, beta = beta, sigma = sigma),
      mle = theta.mle)

# confidence intervals for theta:
# asymptotic: takes a bit of human effort
# bootstrap: lazy, make the computer work

nboot <- 1e4
system.time({
  theta.boot <- replicate(nboot, {
    # non-parametric bootstrap
    D.boot <- D[sample(n, replace = TRUE),]
    Y1.boot <- D.boot[D.boot[,"delta"] == 1,c("x","y","z")]
    X0.boot <- D.boot[D.boot[,"delta"] == 0, c("x","y")]
    V.boot <- mvn0.fit(Y1 = Y1.boot, X0 = X0.boot)
    V2theta(V.boot)
  })
})

# for each element of theta.hat, calculate the bootstrap confidence
# interval.
boot.CI <- function(theta.hat, theta.boot, conf = .95) {
  conf <- 1-conf
  conf <- c(conf/2, 1-conf/2)
  quantile(theta.hat - theta.boot, probs = conf) + theta.hat
}

# bootstrap intervals cover the true parameter value
# contrast this to throwing out missing data with 1e6 (i.e., 100x more)
# observations, for which CI's were way off.
theta.CI <- rbind(est = theta.mle,
                  sapply(1:length(theta.mle), function(ii) {
                    boot.CI(theta.mle[ii], theta.boot[ii,])
                  }))
signif(rbind(true = theta0, theta.CI), 3)
