#--- heteroscedastic linear modeling of SENIC dataset --------------------------

source("hlm-functions.R")
require(numDeriv)
senic <- read.csv("senic.csv")

# homoscedastic error model (usual linear regression)
# regress length on all other variables in dataset _except_ region
# see ?formula for meaning of ./+/-/:/*/1 in "formula" objects.
M0 <- lm(length ~ . - region, data = senic)

# heteroscedastic error model
y <- senic$length # response variable
X <- model.matrix(M0) # mean covariate matrix
W <- X[,1:5] # log-variance covariate matrix
M1 <- hlm.fit(y = y, X = X, W = W)

#--- compare beta's -------------------------------------------------------

signif(rbind(M0 = coef(M0), M1 = M1$beta), 2)

# let's see if differences are significant, i.e.,
# outside confidence intervals

# standard errors
# M0: pre-calculated by lm
se0 <- sqrt(diag(vcov(M0)))
# M1: use Observed Fisher Information
se1 <- hessian(func = function(theta) {
  hlm.loglik(beta = theta[1:ncol(X)], gamma = theta[ncol(X) + 1:ncol(W)],
             y = y, X = X, W = W)
}, x = c(M1$beta, M1$gamma))
se1 <- sqrt(diag(solve(-se1))[1:ncol(X)])

# display nicely
disp <- rbind(rbind(M0 = coef(M0), `2*se` = 2*se0,
                    M1 = M1$beta, `2*se` = 2*se1))
signif(disp, 2) # estimate of each model pretty much within 2*se of the other

#--- compare prediction intervals -------------------------------------------------

M0.pred <- predict(M0,
                   interval = "prediction", level = .95)
M1.pred <- hlm.predict(beta = M1$beta, gamma = M1$gamma,
                       X = X, W = W, level = .95)

# plot by sorted M0 estimate
ord <- order(M0.pred[,"fit"])
n <- nrow(M0.pred)
x <- 1:n
Y <- cbind(M0.pred, M1.pred)[ord,]

# empty plot to set x/y limits
par(mfrow = c(1,1), mar = c(4, 4.5, 1, .5) + .1)
plot(0, type = "n", xlim = range(x), ylim = range(Y),
     xlab = "Sorted Observations", ylab = "Length of Stay (days)")
# colors, line width and line styles
lpars <- data.frame(col = rep(c("black", "red"), each = 3),
               lty = c(1,2,2),
               lwd = 4, # notice recycling of elements
               stringsAsFactors = FALSE) # keep col as string
for(ii in 1:ncol(Y)) {
  lines(x = x, y = Y[,ii],
        col = lpars$col[ii], lty = lpars$lty[ii], lwd = lpars$lwd[ii])
}
points(x, senic$length[ord], pch = 16, col = "blue") # observed data
lgd <- c("\"LM\"", "\"HLM\"", "y[obs]",
         "E*group(\"[\",y*\" | \"*list(bold(x),bold(w)),\"]\")", "\"95% pred. int.\"")
legend("topleft", legend = parse(text = lgd),
       pt.bg = c("black", "red", "black", "black", "black"),
       col = c("black", "black", "blue", "black", "black"),
       pch = c(22, 22, 16, NA, NA),
       lty = c(NA, NA, NA, 1, 2),
       lwd = c(NA, NA, NA, 2, 2), seg.len = 1.5)


# interval lengths
M0.IL <- (M0.pred[,"upr"]-M0.pred[,"lwr"])/2
M1.IL <- (M1.pred[,"upr"]-M1.pred[,"lwr"])/2

signif(rbind(M0 = c(mean = mean(M0.IL), sd = sd(M0.IL)),
             M1 = c(mean = mean(M1.IL), sd = sd(M1.IL))),2)

# fraction of observations for which heteroscedastic model
# had smaller prediction interval
mean(M1.IL < M0.IL)

#--- scratch ---------------------------------------------------------------

# get senic data
# senic <-scan("http://stat.ethz.ch/Teaching/Datasets/senic.dat",
# what=list(id=0,length=0,age=0,inf=0,cult=0,xray=0,beds=0,school=0,
# region=0,pat=0,nurs=0,serv=0))
# senic <- as.data.frame(senic)[-1]
# write.csv(senic, file = "senic.csv", row.names = FALSE)
