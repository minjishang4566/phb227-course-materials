##ctrl x 2; esc x for shell
setwd("~/Desktop/myfiles/Teaching/PHB227/PHB227-2026Fall")
#load(".RData")
options(width=60)  # make the printing fit on the page
set.seed(1121)   # make the results repeatable
 
require(ggplot2)
require(xtable)
library(ISLR2)
library(knitr)
library(MASS)
library(class)

##dataset 
Heart <- read.csv("Heart.csv", header = TRUE, sep=",")
names(Heart)[-1]
Heart$AHD <- as.factor(Heart$AHD)
summary(Heart[,c(2,3,5,6,9,15)])

##basic plots
par(mfrow=c(2,2))
for (i in c(2,5:6,9)){
          boxplot(Heart[,i] ~ AHD, data =Heart, 
main=names(Heart)[i])
} 

##########################################################
##LOGISTIC
xnams <- c("Age", "Sex", "RestBP", "Chol", "MaxHR")
fmla <- fmla <- as.formula(paste("as.factor(AHD) ~ ", paste(xnams, collapse= "+")))
logr.fit <- glm(fmla,data = Heart, family=binomial)
summary(logr.fit)$coef

##evaluating logistic model at a new datum
tmpdata <- c(1,54,1,100, 220, 150)

##calculated prob
exp(sum(logr.fit$coef*tmpdata))/(1+exp(sum(logr.fit$coef*tmpdata)))

tmpdata2 <- data.frame(Age=54,Sex=1,RestBP=100, Chol=220, MaxHR=150)
predict(logr.fit, tmpdata2,type = "response")

hrt.probs <- predict(logr.fit, type = "response")
summary(hrt.probs)

par(mfrow=c(1,1))
boxplot(hrt.probs ~ Heart$AHD, 
main="Logistic Reg-based Probability of AHD")

tbl <- table(Heart$AHD, hrt.probs>0.5)

crrcls <- (tbl[1,1] + tbl[2,2])/sum(tbl)

tbl

crrcls 

sens <- tbl[2,2]/(tbl[2,1] + tbl[2,2]) 
spec <- tbl[1,1]/(tbl[1,1] + tbl[1,2]) 
ppv <- tbl[2,2]/(tbl[1,2] + tbl[2,2]) 
npv <- tbl[1,1]/(tbl[1,1] + tbl[2,1]) 



Accuracy.logr <- c(crrcls,sens,spec,ppv,npv)
dd <- data.frame(Accuracy.logr, row.names = c("Correct Classification", "Sensitivity", "Specificity", "Positive Predictive Value", "Negative Predictive value"))
#names(dd) <- "Accuracy"
round(dd,3)

#table(Heart$AHD, hrt.probs>0.3)  ##improves sens, but higher misclassification
#table(Heart$AHD, hrt.probs>0.8)

###################################################################
##LDA
lda.fit <- lda(fmla,data = Heart)
round(lda.fit$prior,3)
round(lda.fit$means,3)


##lda checks
###pi_k
#table(Heart$AHD)/length(Heart$AHD)

##means
#tapply(Heart$Age, list(Heart$AHD), mean)
#tapply(Heart$Chol, list(Heart$AHD), mean)

#table(Heart$Sex, Heart$AHD)
#prop.table(table(Heart$Sex, Heart$AHD),2)  ##col %s


round(coef(lda.fit),3) ##compare to logistic

##evaluating lda model at a new datum
tmpdata2 <- data.frame(Age=54,Sex=1,RestBP=100, Chol=220, MaxHR=150)
predict(lda.fit,tmpdata2)

lda.pred <- predict(lda.fit,Heart)
#table(Heart$AHD, lda.pred$class)
summary(lda.pred$posterior[,2])  ##posterior prob of Y=1

#par(las = 1,mar = c(4, 4, 0.1, 0.1))  # tick labels direction
boxplot(lda.pred$posterior[,2] ~ Heart$AHD, 
main="LDA-derived Probability of AHD")

tbl <- table(Heart$AHD, lda.pred$posterior[,2]>0.5)

crrcls <- (tbl[1,1] + tbl[2,2])/sum(tbl)
sens <- tbl[2,2]/(tbl[2,1] + tbl[2,2]) 
spec <- tbl[1,1]/(tbl[1,1] + tbl[1,2]) 
ppv <- tbl[2,2]/(tbl[1,2] + tbl[2,2]) 
npv <- tbl[1,1]/(tbl[1,1] + tbl[2,1]) 

tbl

dd$Accuracy.lda <- c(crrcls,sens,spec,ppv,npv)

round(dd,3)

#######################################################################
##QDA

qda.fit <- qda(fmla,data = Heart)
round(qda.fit$prior,3)
round(qda.fit$means,3)

##evaluating QDA moldel at a new datum
tmpdata2 <- data.frame(Age=54,Sex=1,RestBP=100, Chol=220, MaxHR=150)
predict(qda.fit,tmpdata2)


qda.pred <- predict(qda.fit,Heart)
summary(qda.pred$posterior[,2])  ##posterior prob of Y=1

boxplot(qda.pred$posterior[,2] ~ Heart$AHD, 
main="QDA-derived Probability of AHD")

tbl <- table(Heart$AHD, qda.pred$posterior[,2]>0.5)

crrcls <- (tbl[1,1] + tbl[2,2])/sum(tbl)
sens <- tbl[2,2]/(tbl[2,1] + tbl[2,2]) 
spec <- tbl[1,1]/(tbl[1,1] + tbl[1,2]) 
ppv <- tbl[2,2]/(tbl[1,2] + tbl[2,2]) 
npv <- tbl[1,1]/(tbl[1,1] + tbl[2,1]) 

tbl

dd$Accuracy.qda <- c(crrcls,sens,spec,ppv,npv)
round(dd,3)

#########################################################
##k-NN
##k=5
knn.pred <- knn(Heart[,xnams], Heart[,xnams],Heart$AHD, k=5)
#tbl <- table(knn.pred,Heart$AHD)
tbl <- table(Heart$AHD, knn.pred)
round((tbl[1,1] + tbl[2,2])/sum(tbl),2)

crrcls <- (tbl[1,1] + tbl[2,2])/sum(tbl)
sens <- tbl[2,2]/(tbl[2,1] + tbl[2,2]) 
spec <- tbl[1,1]/(tbl[1,1] + tbl[1,2]) 
ppv <- tbl[2,2]/(tbl[1,2] + tbl[2,2]) 
npv <- tbl[1,1]/(tbl[1,1] + tbl[2,1]) 


dd$Accuracy.knn5 <- c(crrcls,sens,spec,ppv,npv)
round(dd,3)

##(Age=54,Sex=1,RestBP=100, Chol=220, MaxHR=150)
tmpdata3 <- c(54,1,100,220,150)
knn(Heart[,xnams],tmpdata3,Heart$AHD, k=5)
##################################


###$k=1:$
knn.pred <- knn(Heart[,xnams], Heart[,xnams],Heart$AHD, k=1)
tbl <- table(knn.pred,Heart$AHD)
(tbl[1,1] + tbl[2,2])/sum(tbl)


####$k=1$, predict on "new" set:
knn.pred <- knn(Heart[1:230,xnams], Heart[231:303,xnams],Heart[1:230,15], k=1)
tbl <- table(knn.pred,Heart$AHD[231:303])
round((tbl[1,1] + tbl[2,2])/sum(tbl),2)


####$k=5$, predict on "new" set:
knn.pred <- knn(Heart[1:230,xnams], Heart[231:303,xnams],Heart[1:230,15], k=5)
tbl <- table(knn.pred,Heart$AHD[231:303])
round((tbl[1,1] + tbl[2,2])/sum(tbl),2)


################################
##MULTI CLASS EXAMPLE

##create categorical MaxHR
Heart$HRtrtle <- cut(Heart$MaxHR, c(0,quantile(Heart$MaxHR,c(0.33,0.67)),max(Heart$MaxHR)+1))

summary(Heart$HRtrtle)

levels(Heart$HRtrtle) <- c("l","m","h")
summary(Heart$HRtrtle)

##I. multinomial logistic reg

library(nnet)
mlogist <- multinom(HRtrtle~Age+RestBP+Chol+Sex, data=Heart, na.actoin=na.omit)
summary(mlogist)
#exp(coef(mlogist))

##inference
#summary(mlogist)$coeff/summary(mlogist)$standard.error


##predicted probabilities
fitted(mlogist)[1:20,]
predict(mlogist)[1:20]


#table(predict(mlogist),Heart$HRtrtle)

##accuracy
#sum(diag(table(predict(mlogist),Heart$HRtrtle)))/sum(table(predict(mlogist),Heart$HRtrtle))

###################################################
##II. Multiclass LDA

mlda <- lda(HRtrtle~Age+RestBP+Chol+Sex,data = Heart)
mlda


predmlda <- predict(mlda,Heart)
predmlda$posterior[1:20,]
predmlda$class[1:20]


##Accuracy
#table(predmlda$class,Heart$HRtrtle)
#sum(diag(table(predmlda$class,Heart$HRtrtle)))/sum(table(predmlda$class,Heart$HRtrtle))

##########################
##III. multiclass QDA


mqda <- qda(HRtrtle~Age+RestBP+Chol+Sex,data = Heart)
#mqda


predmqda <- predict(mqda,Heart)
predmqda$posterior[1:20,]
predmqda$class[1:20]

##Accuracy
#table(predmqda$class,Heart$HRtrtle)
#sum(diag(table(predmqda$class,Heart$HRtrtle)))/sum(table(predmqda$class,Heart$HRtrtle))

