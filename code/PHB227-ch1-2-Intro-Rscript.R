setwd("~/Desktop/myfiles/Teaching/PHB227/PHB227-2026Fall")
library(ISLR2)
library(knitr)
library(MASS)
library(FNN)
#library(class)
#library(glmnet)
#library(methods)
#library(tree)
#library(randomForest)
library(ggplot2)


options(width=60)  # make the printing fit on the page
set.seed(1121)   # make the results repeatable
#require(PASWR)
require(ggplot2)
require(xtable)

##DATASET
##read in Heart dataset
Heart <- read.csv("Heart.csv", header = TRUE, sep=",")
names(Heart)[-1]
Heart$AHD <- as.factor(Heart$AHD)
summary(Heart[,c(2,6,9,15)])



#Outcome: MaxHR; Predictors: Age, Chol\\
#We want to approximate MaxHR as a function of Age and Chol.


##PLOTS
par(mfrow=c(1,2))
#par(las = 1,mar = c(4, 4, 0.1, 0.1), cex=0.5)  # tick labels direction
plot(Heart$Age, Heart$MaxHR)
plot(Heart$Chol,Heart$MaxHR)
#loess(Heart$Age, Heart$RestBP)
#boxplot(Chol~Sex,data=Heart)



##3-d plot 

par(mfrow=c(1,1))
#par(las = 1,mar = c(4, 4, 0.1, 0.1), cex=0.5)  # tick labels direction
library(scatterplot3d)
with(data=Heart, scatterplot3d(Age,Chol,MaxHR))

##PARAMETRIC MODEL
##fit linear model
hrfit <- lm(MaxHR~Age + Chol,data=Heart)
coef(hrfit)


##predict MaxHR for a 54-year old with cholesterol level of 171 mg/dl

##prediction in R
predict(hrfit, newdata = data.frame(Age=54, Chol=171))


##linear prediction: by hand
sum(coef(hrfit) *c(1,54,171))


##plot linearfit
x11 <- seq(min(Heart$Age),max(Heart$Age),(max(Heart$Age) - min(Heart$Age))/100)
x22 <- seq(min(Heart$Chol),max(Heart$Chol),(max(Heart$Chol) - min(Heart$Chol))/100)
fnn1 <- function(xx,yy){predict(hrfit, newdata = data.frame(Age=xx, Chol=yy))}
zz <- outer(x11,x22,fnn1)
pp <- persp(x11,x22,zz,theta = 30, phi = 20, expand = 0.6, shade = 0.1,col="lightblue",xlab="Age", ylab="Chol", zlab="MaxHR",xlim=c(min(Heart$Age),max(Heart$Age)),ylim=c(min(Heart$Chol),max(Heart$Chol)), zlim=c(min(Heart$MaxHR),max(Heart$MaxHR)),ticktype="detailed", nticks=10,cex=0.5,main="Linear")
oobs <- trans3d(Heart$Age,Heart$Chol,Heart$MaxHR,pp)
ppred <- trans3d(Heart$Age,Heart$Chol,fitted(hrfit),pp)
points(oobs, col="blue",pch=16)
segments(oobs$x, oobs$y, ppred$x, ppred$y)

##NONPARAMETRIC
##loess fit
hrfit2 <- loess(MaxHR~Age + Chol,data=Heart) ##loess
#hrfit2

##knn fit
library(FNN)
kknum <- 4
knn.pred2 <- knn.reg(cbind(Heart$Age,Heart$Chol), cbind(x11,x22),Heart$MaxHR, k=kknum) ##knn
x11[1:10]
x22[1:10]
knn.pred2$pred[1:10]


##plot loess and knn
#par(mfrow=c(1,2))
x11 <- seq(min(Heart$Age),max(Heart$Age),(max(Heart$Age) - min(Heart$Age))/100)
x22 <- seq(min(Heart$Chol),max(Heart$Chol),(max(Heart$Chol) - min(Heart$Chol))/100)

##loess plot
fnn2 <- function(xx,yy){predict(hrfit2, newdata = data.frame(Age=xx, Chol=yy))}
zz2 <- outer(x11,x22,fnn2)
pp2 <- persp(x11,x22,zz2,theta = 30, phi = 20, expand = 0.6, shade = 0.1,col="lightblue",xlab="Age", ylab="Chol", zlab="MaxHR",xlim=c(min(Heart$Age),max(Heart$Age)),ylim=c(min(Heart$Chol),max(Heart$Chol)), zlim=c(min(Heart$MaxHR),max(Heart$MaxHR)),ticktype="detailed", nticks=10,cex=0.3,main="Loess")
oobs2 <- trans3d(Heart$Age,Heart$Chol,Heart$MaxHR,pp2)
ppred2 <- trans3d(Heart$Age,Heart$Chol,fitted(hrfit2),pp2)
points(oobs2, col="blue",pch=16)
segments(oobs2$x, oobs2$y, ppred2$x, ppred2$y)

###knn plot


fnn3 <- function(xx,yy){knn.reg(cbind(Heart$Age,Heart$Chol), cbind(xx,yy),Heart$MaxHR, k=kknum)$pred}

zz3 <- outer(x11,x22,fnn3)
pp3 <- persp(x11,x22,zz3,theta = 30, phi = 20, expand = 0.6, shade = 0.1,col="lightblue",xlab="Age", ylab="Chol", zlab="MaxHR",xlim=c(min(Heart$Age),max(Heart$Age)),ylim=c(min(Heart$Chol),max(Heart$Chol)), zlim=c(min(Heart$MaxHR),max(Heart$MaxHR)),ticktype="detailed", nticks=10,cex=0.3,main="k-NN")
oobs3 <- trans3d(Heart$Age,Heart$Chol,Heart$MaxHR,pp3)
ppred3 <- trans3d(Heart$Age,Heart$Chol,fitted(hrfit2),pp3)
points(oobs3, col="blue",pch=16)
segments(oobs3$x, oobs3$y, ppred3$x, ppred3$y)


##knn prediction: by hand
##Find distance from (x0,y0) = (54,171) to all pairs (Age,Chol) in Heart dataset
##identify the kknum pairs of (Age,Chol) which  are closest to (x0,y0) = (54,171)
##find MaxHR corresponding to these closest kknum
##Average these kknum MaxHR values
# (Heart$Age-54)^2 + (Heart$Chol - 171)^2)
(sort((Heart$Age-54)^2 + (Heart$Chol - 171)^2))[1:kknum]

Heart$Age[order((Heart$Age-54)^2 + (Heart$Chol - 171)^2)][1:kknum]
Heart$Chol[order((Heart$Age-54)^2 + (Heart$Chol - 171)^2)][1:kknum]
Heart$MaxHR[order((Heart$Age-54)^2 + (Heart$Chol - 171)^2)][1:kknum]
mean(Heart$MaxHR[order((Heart$Age-54)^2 + (Heart$Chol - 171)^2)][1:kknum])


##Prediction via R functions:

##k-nn

knn.reg(cbind(Heart$Age,Heart$Chol), cbind(54,171),Heart$MaxHR, k=kknum)


##loess
predict(hrfit2, newdata = data.frame(Age=54, Chol=171))

####CLASSIFICATION

#Outcome: AHD
#Predictors: Age, Chol
#We want to approximate AHD as a function of Age, Chol.

#par(las = 1,mar = c(4, 4, 0.1, 0.1), cex=0.5)  # tick labels direction
par(mfrow=c(1, 2))
for (i in c(2,6)){
          boxplot(Heart[,i] ~ AHD, data =Heart, main=names(Heart)[i])
} 







