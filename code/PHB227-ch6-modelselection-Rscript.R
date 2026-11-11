setwd("~/Desktop/myfiles/Teaching/PHB227/PHB227-2026Fall")
getwd()
library(ISLR2)
library(knitr)
library(MASS)
library(class)
library(glmnet)
library(methods)

options(width=60)  # make the printing fit on the page
set.seed(1121)   # make the results repeatable
#require(PASWR)
require(ggplot2)
#require(xtable)


#####DATASET
Heart <- read.csv("Heart.csv", header = TRUE, sep=",")
Heart$AHD <- as.factor(Heart$AHD)
xnams2 <- c("Age", "RestBP", "Chol", "MaxHR")
fmla2 <- as.formula(paste("AHD ~ ", paste(xnams2, collapse= "+")))
      


###RIDGE and LASSO fit
xnams2
xxx <- as.matrix(Heart[,xnams2])
yyy <- Heart$AHD

gridd <- exp(seq(2,-6,-0.5))  ##lambda values

##ridge fit
rdg.fit <- glmnet(xxx,yyy,family="binomial",alpha=0,lambda=gridd)

##lasso fit
lso.fit <- glmnet(xxx,yyy,family="binomial",alpha=1,lambda=gridd)

##lambda
gridd

rdg.fit$lambda

coef(rdg.fit)

coef(lso.fit)

##very shrunk
rdg.fit$lambda[1]
coef(rdg.fit)[,1][xnams2]
coef(lso.fit)[,1][xnams2]


##almost no shrinkage
rdg.fit$lambda[length(gridd)]
coef(rdg.fit)[,length(gridd)][xnams2]
coef(lso.fit)[,length(gridd)][xnams2]


##logistic coeff: no shrinkage!
glm(fmla2, data=Heart,family="binomial")$coef[xnams2]


####coef output for each lambda
names(rdg.fit) ###https://cran.r-project.org/web/packages/glmnet/glmnet.pdf

rdg.fit   ###note df column and contrast with lasso
rdg.fit$a0
rdg.fit$beta #coef without intercept

lso.fit
lso.fit$a0
lso.fit$beta

#help(glmnet)
####software chooses lambda
#rdg.fit2 <- glmnet(xxx,yyy,family="binomial",alpha=0)
#lso.fit2 <- glmnet(xxx,yyy,family="binomial",alpha=1)
##################PLOTS 
par(mfrow=c(1,2))
par(las = 1,mar = c(4, 4, 0.1, 0.1), cex=0.5)  # tick labels direction
xnams2

lbs_fun <- function(fitt, ...) {
  L <- length(fitt$lambda)
  xx <- log(fitt$lambda[L])
  yy <- fitt$beta[, L]
  llabs <- names(yy)
  text(xx, yy, labels=llabs, ...)
#  legend('topright', legend=llabs, col=1:length(llabs), lty=1) 
}

plot(rdg.fit,xvar="lambda",label=TRUE,col=1:dim(coef(rdg.fit))[1],lwd=3)
lbs_fun(rdg.fit)
abline(h=0)

plot(lso.fit,xvar="lambda",label=TRUE,col=1:dim(coef(lso.fit))[1],lwd=3)
lbs_fun(rdg.fit)
#plot(lso.fit,xvar="lambda",label=TRUE,lwd=3)
abline(h=0)
#############################################################

####PREDICTION
nn <- dim(Heart)[1]
set.seed(1996)
trn <- sample(1:nn,ceiling(nn/2))


##ridge train fit, then test
rdgTR.fit <- glmnet(xxx[trn,],yyy[trn],family="binomial",alpha=0,lambda=gridd)

rdg.TST <- predict(rdgTR.fit, s=c(gridd[1],gridd[7],gridd[17]), newx=xxx[-trn,],type="class")

gridd[1]
table(rdg.TST[,1],Heart[-trn,"AHD"])

gridd[7]
table(rdg.TST[,2],Heart[-trn,"AHD"])

gridd[17]
table(rdg.TST[,3],Heart[-trn,"AHD"])

rdg.TST

##extracting predicted prob
rdg.TST2 <- predict(rdgTR.fit, s=c(gridd[1],gridd[7],gridd[17]), newx=xxx[-trn,],type="response")

rdg.TST2 ##note that col 1 prob are mostly less than 0.5..very shrunk!
#################################
##lasso train fit, then test
lsoTR.fit <- glmnet(xxx[trn,],yyy[trn],family="binomial",alpha=1,lambda=gridd)

lso.TST <- predict(lsoTR.fit, s=c(gridd[1],gridd[7],gridd[17]), newx=xxx[-trn,],type="class")

gridd[1] # only intercept
table(lso.TST[,1],Heart[-trn,"AHD"])

gridd[7]
table(lso.TST[,2],Heart[-trn,"AHD"])

gridd[17]
table(lso.TST[,3],Heart[-trn,"AHD"])

lso.TST

##extracting predicted prob
lso.TST2 <- predict(lsoTR.fit, s=c(gridd[1],gridd[7],gridd[17]), newx=xxx[-trn,],type="response")

lso.TST2  ## col1 and 2 have all pred prob < 0.5; lasso shrinks to 0. 

coef(lsoTR.fit)[,c(1,7,17)]
coef(rdgTR.fit)[,c(1,7,17)]
################ CROSS-VALIDATION TO SELECT LAMBDA
##cross-validation to select lambda
set.seed(2446)
cv.rdgeg <- cv.glmnet(xxx,yyy,family="binomial",alpha=0,lambda=gridd,nfolds=10)

cv.lsoeg <- cv.glmnet(xxx,yyy,family="binomial",alpha=1,lambda=gridd,nfolds=10)
c(cv.rdgeg$lambda.min,cv.lsoeg$lambda.min)
log(c(cv.rdgeg$lambda.min,cv.lsoeg$lambda.min))



par(mfrow=c(1,2))
plot(cv.rdgeg, main= "Ridge")
plot(cv.lsoeg, main="Lasso")


################### OPTIMAL LAMBDA
####Variables (and coefficients) at optimal $\lambda$s for RIDGE
c(cv.rdgeg$lambda.min,cv.rdgeg$lambda.1se)
log(c(cv.rdgeg$lambda.min,cv.rdgeg$lambda.1se))

##Coeff at "best" lambda
coef(cv.rdgeg,s="lambda.min")


##Coeff at  lambda + 1SE
coef(cv.rdgeg,s="lambda.1se")

###Variables (and coefficients) at optimal $\lambda$s for LASSO:
##Coeff at "best" lambda
c(cv.lsoeg$lambda.min,cv.lsoeg$lambda.1se)
log(c(cv.lsoeg$lambda.min,cv.lsoeg$lambda.1se))

coef(cv.lsoeg,s="lambda.min")

##Coeff at  lambda + 1SE
coef(cv.lsoeg,s="lambda.1se")

#names(cv.lsoeg)
#######################################################################
###########
####Using misclassification as the loss function

##cross-validation to select lambda
set.seed(2446)
cv.rdgeg2 <- cv.glmnet(xxx,yyy,family="binomial",alpha=0,lambda=gridd,nfolds=10,type.measure="class")

cv.lsoeg2 <- cv.glmnet(xxx,yyy,family="binomial",alpha=1,lambda=gridd,nfolds=10,type.measure="class")


par(mfrow=c(1,2))
plot(cv.rdgeg2, main = "Ridge")
plot(cv.lsoeg2, main ="Lasso")

##min lambda 
log(c(cv.rdgeg2$lambda.min,cv.lsoeg2$lambda.min))

##compare deviance vs miclassification coefs: same 'coz lambda_min is same
#coef(cv.lsoeg,s="lambda.min")
#coef(cv.lsoeg2,s="lambda.min")



###Kathy's  qt
#log(c(cv.rdgeg2$lambda.min, cv.rdgeg2$lambda.1se))
#cbind(log(cv.rdgeg2$lambda),cv.rdgeg2$cvm, cv.rdgeg2$cvup, cv.rdgeg2$cvlo)
########################################
##SIMULATION: COMPARE LASSO to RIDGE
##http://www4.stat.ncsu.edu/~post/josh/LASSO_Ridge_Elastic_Net_-_Examples.html#generate-data

#################EXAMPLE 1
# Generate data
set.seed(1988)  # Set seed for reproducibility
n <- 100  # Number of observations
p <- 90  #Number of predictors included in model

x <- matrix(rnorm(n*p), nrow=n, ncol=p)

y <- 10 * apply(x[, 1:2], 1, sum) + 
  5 * apply(x[, 3:4], 1, sum) +
  apply(x[, 5:6], 1, sum) +
  rnorm(n)



# Split data into train (2/3) and test (1/3) sets
train_rows <- sample(1:n, .66*n)
x.train <- x[train_rows, ]
x.test <- x[-train_rows, ]

y.train <- y[train_rows]
y.test <- y[-train_rows]

# Fit models on training set

gridd2 <- exp(seq(1,-8,-0.5))  ##lambda values
fit.lasso <- glmnet(x.train, y.train, family="gaussian", alpha=1,lambda=gridd2)
fit.ridge <- glmnet(x.train, y.train, family="gaussian", alpha=0,lambda=gridd2)
#fit.elnet <- glmnet(x.train, y.train, family="gaussian", alpha=.5,lambda=gridd2)
	  

# 10-fold Cross validation for each alpha = 0, 0.5, 1.0
# 
for (i in c(0,5,10)) {
    assign(paste("fit", i, sep=""), cv.glmnet(x.train, y.train, type.measure="mse", 
                                              alpha=i/10,family="gaussian",lambda=gridd2))
}


# Plot
par(mfrow=c(2,2))
# For plotting options, type '?plot.glmnet' in R console

plot(fit0, main="Ridge")
plot(fit10, main="LASSO")
#plot(fit5,main="Elastic Net")


##extract lambda.1se
log(c(fit0$lambda.1se,fit10$lambda.1se))
(c(fit0$lambda.1se,fit10$lambda.1se))

#extract coeff at lambda.1se
cbind(coef(fit.ridge,s=fit0$lambda.1se),coef(fit.lasso,s=fit10$lambda.1se))

##extract lambda.min and coeff at lambda.min
#log(c(fit0$lambda.min,fit10$lambda.min))
#cbind(coef(fit.ridge,s=fit0$lambda.min), coef(fit.lasso,s=fit10$lambda.min))

#compare coeff at lambda.1se and lambda.min
#cbind(coef(fit.ridge,s=fit0$lambda.min), coef(fit.ridge,s=fit0$lambda.1se))
#cbind(coef(fit.lasso,s=fit10$lambda.min),coef(fit.lasso,s=fit10$lambda.1se))

###predict on test set
#lambda.1se
yhat0 <- predict(fit0, s=fit0$lambda.1se, newx=x.test)
yhat10 <- predict(fit10, s=fit10$lambda.1se, newx=x.test)
#yhat5 <- predict(fit5, s=fit5$lambda.1se, newx=x.test)

#lambda.min
#yhat0 <- predict(fit0, s=fit0$lambda.min, newx=x.test)
#yhat10 <- predict(fit10, s=fit10$lambda.min, newx=x.test)
#yhat5 <- predict(fit5, s=fit5$lambda.min, newx=x.test)

#mean sq errors on test set
mse0 <- mean((y.test - yhat0)^2)
mse10 <- mean((y.test - yhat10)^2)
#mse5 <- mean((y.test - yhat5)^2)

mse0
mse10
#mse5

###MSEs over lambda
mse.rdg <- numeric(0)
mse.lso <- numeric(0)

for (i in 1:length(gridd2)){
         yhatt0 <- predict(fit0, s=gridd2[i], newx=x.test)
	 mse.rdg[i] <- mean((y.test - yhatt0)^2)
	 yhatt10 <- predict(fit10, s=gridd2[i], newx=x.test)
	 mse.lso[i] <- mean((y.test - yhatt10)^2)
#cat(gridd2[i],mse.rdg[i],mse.lso[i],"\n")
	 rm(yhatt0,yhatt10)
	 }

##plot mses over lambda
plot(log(gridd2),mse.rdg,type="n")
lines(log(gridd2),mse.rdg,col="blue")

plot(log(gridd2),mse.lso,type="n")
lines(log(gridd2),mse.lso,col="red")

#log(c(fit0$lambda.min,fit10$lambda.min))
#log(c(fit0$lambda.1se,fit10$lambda.1se))

#plot(log(gridd2),mse.lso,type="n",xlim=c(-4,0.5),ylim=c(0,20))
#lines(log(gridd2),mse.lso,col="red")


##plot coefficient paths
par(mfrow=c(1,2))
plot(fit.ridge,xvar="lambda",label=TRUE)
plot(fit.lasso,xvar="lambda",label=TRUE)


################################
###EXAMPLE 2

# Generate data
set.seed(1996)  # Set seed for reproducibility
n <- 100  # Number of observations
p <- 90  # Number of predictors included in model
real_p <- 90  # Number of true predictors (rest have null association)
x <- matrix(rnorm(n*p), nrow=n, ncol=p)
y <- apply(x[,1:real_p], 1, sum) + rnorm(n)

# Split data into train (2/3) and test (1/3) sets
train_rows <- sample(1:n, .66*n)
x.train <- x[train_rows, ]
x.test <- x[-train_rows, ]

y.train <- y[train_rows]
y.test <- y[-train_rows]

# Fit models 
fit.lasso <- glmnet(x.train, y.train, family="gaussian", alpha=1,lambda=gridd2)
fit.ridge <- glmnet(x.train, y.train, family="gaussian", alpha=0,lambda=gridd2)
#fit.elnet <- glmnet(x.train, y.train, family="gaussian", alpha=.5,lambda=gridd2)


# 10-fold Cross validation for each alpha = 0, 1.0
# 
for (i in c(0,5,10)) {
    assign(paste("fit", i, sep=""), cv.glmnet(x.train, y.train, type.measure="mse", 
                                              alpha=i/10,family="gaussian",lambda=gridd2))
}


# Plot:
par(mfrow=c(2,2))
# For plotting options, type '?plot.glmnet' in R console
plot(fit0, main="Ridge")
plot(fit10, main="LASSO")
#plot(fit5, main="Elastic Net")


##lambda.1se and coeff at lambda.1se
log(c(fit0$lambda.1se,fit10$lambda.1se))
cbind(coef(fit.ridge,s=fit0$lambda.1se),coef(fit.lasso,s=fit10$lambda.1se))

##coeff at lambda.min
#log(c(fit0$lambda.min,fit10$lambda.min))
#cbind(coef(fit.ridge,s=fit0$lambda.min), coef(fit.lasso,s=fit10$lambda.min))


#cbind(coef(fit.ridge,s=fit0$lambda.min), coef(fit.ridge,s=fit0$lambda.1se))
#cbind(coef(fit.lasso,s=fit10$lambda.min),coef(fit.lasso,s=fit10$lambda.1se))

##predict on test set

#lambda.1se
yhat0 <- predict(fit0, s=fit0$lambda.1se, newx=x.test)
yhat10 <- predict(fit10, s=fit10$lambda.1se, newx=x.test)
#yhat5 <- predict(fit5, s=fit5$lambda.1se, newx=x.test)

#lambda.min
#yhat0 <- predict(fit0, s=fit0$lambda.min, newx=x.test)
#yhat10 <- predict(fit10, s=fit10$lambda.min, newx=x.test)
#yhat5 <- predict(fit5, s=fit5$lambda.min, newx=x.test)


mse0 <- mean((y.test - yhat0)^2)
mse10 <- mean((y.test - yhat10)^2)
#mse5 <- mean((y.test - yhat5)^2)

mse0
mse10
#mse5

###MSEs over lambda
mse.rdg <- numeric(0)
mse.lso <- numeric(0)
for (i in 1:length(gridd2)){
         yhatt0 <- predict(fit0, s=gridd2[i], newx=x.test)
	 mse.rdg[i] <- mean((y.test - yhatt0)^2)
	 yhatt10 <- predict(fit10, s=gridd2[i], newx=x.test)
	 mse.lso[i] <- mean((y.test - yhatt10)^2)
	 rm(yhatt0,yhatt10)
	 }



plot(log(gridd2),mse.rdg,type="n")
lines(log(gridd2),mse.rdg,col="blue")

plot(log(gridd2),mse.lso,type="n")
lines(log(gridd2),mse.lso,col="red")


##coeff paths
#par(mfrow=c(1,2))
#plot(fit.ridge,xvar="lambda",label=TRUE)
#plot(fit.lasso,xvar="lambda",label=TRUE)

############################################
##EXAMPLE 3

# Generate data
set.seed(19873)  # Set seed for reproducibility
n <- 100  # Number of observations
p <- 50  # Number of predictors included in model


#rhoo <- 0.9 ## corr between predictors
CovMatrixx <- outer(1:p, 1:p, function(x,y) {.7^abs(x-y)})
#CovMatrixx <- matrix(rhoo, nrow=p, ncol=p)
#diag(CovMatrixx) <- 1
#table(CovMatrixx)

x <- mvrnorm(n, rep(0,p), CovMatrixx)
y <- 10*apply(x[, 1:2], 1, sum) + 
      5*apply(x[, 3:4], 1, sum) +
    apply(x[, 5:14], 1, sum) +
    rnorm(n)
  
# Split data into train (2/3) and test (1/3) sets
set.seed(1988)
train_rows <- sample(1:n, .66*n)
x.train <- x[train_rows, ]
x.test <- x[-train_rows, ]

y.train <- y[train_rows]
y.test <- y[-train_rows]

# Fit models 

fit.lasso <- glmnet(x.train, y.train, family="gaussian", alpha=1)
fit.ridge <- glmnet(x.train, y.train, family="gaussian", alpha=0)
fit.elnet <- glmnet(x.train, y.train, family="gaussian", alpha=.5)


# 10-fold Cross validation for each alpha = 0, 1.0
# 
for (i in c(0,5,10)) {
    assign(paste("fit", i, sep=""), cv.glmnet(x.train, y.train, type.measure="mse", 
                                              alpha=i/10,family="gaussian",lambda=gridd2))
}


# Plot:
par(mfrow=c(2,3))
# For plotting options, type '?plot.glmnet' in R console
plot(fit0, main="Ridge")
plot(fit10, main="LASSO")
plot(fit5, main="Elastic Net")

#par(mfrow=c(1,1))
#plot(fit.ridge,xvar="lambda",label=TRUE)
#plot(fit.lasso,xvar="lambda",label=TRUE)

#cbind(coef(fit.ridge,s=fit0$lambda.1se),
#coef(fit.lasso,s=fit10$lambda.1se),
#coef(fit.elnet,s=fit5$lambda.1se)
#)



yhat0 <- predict(fit0, s=fit0$lambda.1se, newx=x.test)
yhat10 <- predict(fit10, s=fit10$lambda.1se, newx=x.test)
yhat5 <- predict(fit5, s=fit5$lambda.1se, newx=x.test)

#yhat0 <- predict(fit0, s=fit0$lambda.min, newx=x.test)
#yhat10 <- predict(fit10, s=fit10$lambda.min, newx=x.test)
#yhat5 <- predict(fit5, s=fit5$lambda.min, newx=x.test)


mse0 <- mean((y.test - yhat0)^2)
mse10 <- mean((y.test - yhat10)^2)
mse5 <- mean((y.test - yhat5)^2)

mse0
mse10
mse5

###MSEs over lambda
mse.rdg <- numeric(0)
mse.lso <- numeric(0)
mse.elnet <- numeric(0)

for (i in 1:length(gridd2)){
         yhatt0 <- predict(fit0, s=gridd2[i], newx=x.test)
	 mse.rdg[i] <- mean((y.test - yhatt0)^2)
	 yhatt10 <- predict(fit10, s=gridd2[i], newx=x.test)
	 mse.lso[i] <- mean((y.test - yhatt10)^2)
	 yhatt5 <- predict(fit5, s=gridd2[i], newx=x.test)
	 mse.elnet[i] <- mean((y.test - yhatt5)^2)
	 rm(yhatt0,yhatt5,yhatt10)
	 }



plot(log(gridd2),mse.rdg,type="n")
lines(log(gridd2),mse.rdg,col="blue")

plot(log(gridd2),mse.lso,type="n")
lines(log(gridd2),mse.lso,col="red")


plot(log(gridd2),mse.elnet,type="n")
lines(log(gridd2),mse.elnet,col="green")

#plot(log(gridd2),mse.lso,type="n",xlim=c(0,0.5),ylim=c(0,10))
#lines(log(gridd2),mse.lso,col="red")
