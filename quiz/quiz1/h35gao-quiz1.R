source("mvn-functions.R")
source("h35gao-functions.R")

### 1(b)
# randomly generate dimension n in range [0,100]
n <- ceiling(runif(1)*100)
mu <- rnorm(n)
# guaranteed +ve def matrix
V <- crossprod(matrix(rnorm(n^2), n, n))

# generate 1000 observations
N <- 1000
X <- rmvn(n = N, mu = mu, V = V)
# randomly generate indicator vector
ind1 <- sample(c(TRUE,FALSE),n,TRUE)

X1 <- X[,ind1,drop=FALSE]
X2 <- X[,!ind1,drop=FALSE]

result <- cmvn(mu,V,t(X1),ind1)
cmu2 <- result[[1]]
cSigma2 <- result[[2]]

bool <- TRUE
# tolerance
tol <- 1e-10
# verify result
# loop through each observation
for (i in 1:N) {
  xmvn <- dmvn(X[i,],mu,V,log=TRUE)
  x1mvn <- dmvn(X1[i,],mu[ind1],V[ind1,ind1,drop=FALSE],log=TRUE)
  x2mvn <- dmvn(X2[i,],cmu2[,i],cSigma2,log=TRUE)
  bool <- bool & (abs(xmvn-x1mvn-x2mvn) <= tol)
}
cat("The result is:", bool)


### 2(c)
N <- sample(2:10, 1) # number of timepoints
# generate consistent sseq and tseq
stseq <- cumsum(rexp(2*N)) # generate both simultaneously
sseq <- stseq[seq(from = 1, to = 2*N, by = 2)] # odd elements
tseq <- stseq[seq(from = 2, to = 2*N, by = 2)] # even elements
# generate the variance matrix of (Bt,Bs)
seq <- c(tseq,sseq)
Sigma <- bmV(seq)
# use cvmn from Q1 to calculate the conditional mean and variance
mu <- rep(0,2*N)
Sigmat <- bmV(tseq)
Bt <- rmvn(n = 1, mu = rep(0,N), V = Sigmat)
ind <- c(rep(TRUE,N),rep(FALSE,N))
cmvnr <- cmvn(mu,Sigma,t(Bt),ind)
# use cbm to calculate the conditional mean and standard deviations
cbmr <- cbm(sseq,tseq,Bt)

# verify results
tol <- 1e-10
bool <- (norm(cmvnr[[1]]-cbmr[[1]]) <= tol)
bool <- bool & (norm(cmvnr[[2]]-diag(cbmr[[2]]^2)) <= tol)
cat("The result is:", bool)


