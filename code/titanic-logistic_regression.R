#--- analysis of the Titanic dataset ---------------------------------------

dimnames(Titanic) # included by default with R

# 4-way contingency table
mosaicplot(Titanic, main = "Survival on the Titanic",
           color = c("red", "blue"), cex = 2)

# 3-way contingency table (ignoring gender)
mosaicplot(~ Class + Age + Survived, data = Titanic,
           main = "Survival on the Titanic", color = c("red", "blue"), cex = 2)

#--- logistic regression ---------------------------------------------------

# convert the dataset to the appropriate format
# i.e., a data-frame with 5 columns:
# Class, Age, Sex (covariates)
# Total, Survived (response)

# data.frame of passengers per group
total <- as.data.frame(ftable(Sex ~ Class + Age, data = Titanic))
names(total)[4] <- "Total"
total

# data.frame of survivors only
# notice that the order of variables is different
surv <- as.data.frame(ftable(Class ~ Sex + Age, data = Titanic[,,,"Yes"]))
names(surv)[4] <- "Survived"
surv

# merge the two, order-safe way
titanic <- merge(x = surv, y = total, by = c("Class", "Age", "Sex"))
titanic

# logistic regression; main effects only
M <- glm(cbind(Survived, Total-Survived) ~ Class + Age + Sex,
         data = titanic, family = "binomial")

# coefficients
coef(M) # short form

coef(summary(M)) # long form

# parameter interpretation:
# (Intercept): ilogit(beta0) = expected survival probability of
# ClassCrew: exp(beta_ClassCrew) = expected odds-ratio between

# other useful glm commands

signif(vcov(M),2) # Inverse of Observed Fisher Information

predict(M, type = "response") # estimated survival probabilities

# display these nicely

# empirical survival probabilities
surv.emp <- titanic$Survived / titanic$Total
# glm survival probabilities
surv.glm <- predict(M, type = "response", newdata = titanic)


surv.fit <- cbind(titanic[c("Class", "Age", "Sex", "Total")],
                  ProbEmp = signif(surv.emp, 2),
                  ProbGLM = signif(surv.glm, 2))
surv.fit

# Difference between empirical and fitted in cells with large Total indicates
# that glm is pooling information too aggressively,
# i.e., should include interaction term(s).
# on the other hand this will increase uncertainty in cells with small Total.
# This is an example of the so-called "Bias-Variance trade-off".
