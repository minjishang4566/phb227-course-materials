library(ISLR2)
library(knitr)
library(MASS)
library(class)

setwd("~/Desktop/myfiles/Teaching/PHB227/PHB227-2026Fall")
set.seed(112)   # make the results repeatable
#require(PASWR)
require(ggplot2)
#require(xtable)
Heart <- read.csv("Heart.csv", header = TRUE, sep=",")
Heart$AHD <- as.factor(Heart$AHD)
xnams <- c("Age", "Sex", "RestBP", "Chol", "MaxHR")
fmla <- as.formula(paste("AHD ~ ", paste(xnams, collapse= "+")))


##### TRAINING SET-VALIDATION SET SPLIT 
nn <- dim(Heart)[1]
set.seed(1234)

trn <- sample(1:nn,ceiling(nn/2))
logrTR.fit <- glm(fmla,data = Heart, family=binomial,subset= trn)

# training error
trn.pred <- predict(logrTR.fit, Heart[trn,], type = "response")
trn.conf <- table(trn.pred>0.5,Heart[trn,"AHD"])

trn.error <- 1- (trn.conf[1,1] + trn.conf[2,2])/sum(trn.conf)

list(trn.conf, round(trn.error,2))


# testing error, compare with training error
tst.pred <- predict(logrTR.fit, Heart[-trn,], type = "response")

tst.conf <- table(tst.pred>0.5,Heart[-trn,"AHD"])
tst.error <- 1- (tst.conf[1,1] + tst.conf[2,2])/sum(tst.conf)

list(tst.conf, round(tst.error,2))

#Now let's repeat this $10$ times, i.e., randomly split training and validation.

##train/valid
nn <- dim(Heart)[1]
tst.error <- rep(0,10)
trn.error <- rep(0,10)

for (i in 1:10){
#set.seed(1234)
trn <- sample(1:nn,ceiling(nn/2))

logrTR.fit <- glm(fmla,data = Heart, family=binomial,subset= trn)

trn.pred <- predict(logrTR.fit, Heart[trn,], type = "response")
trn.conf <- table(trn.pred>0.5,Heart[trn,"AHD"])
trn.error[i] <- 1- (trn.conf[1,1] + trn.conf[2,2])/sum(trn.conf)


tst.pred <- predict(logrTR.fit, Heart[-trn,], type = "response")
tst.conf <- table(tst.pred>0.5,Heart[-trn,"AHD"])
tst.error[i] <- 1- (tst.conf[1,1] + tst.conf[2,2])/sum(tst.conf)
}

round(trn.error,2)
round(tst.error,2)

summary(logrTR.fit)
####################################################################
##LOOCV
library(boot)
logr.fit <- glm(fmla,data = Heart, family=binomial)
cost.fn <- function(r, pii = 0) mean(abs(r-pii) > 0.5)

cv.glm(Heart,logr.fit,cost.fn)$delta ##logistic model (including bias-corrected version)

########################################################

##Loki's CV function
##dataa is dataset
##methodd is classification methods: "ldaa"=lda, "logisticc"=logistic;"qdaa"=qda;
##stsdd = TRUE means setseed
##fmlaa is model formula
###responsee is numeric outcome as 0 -1 variable 
##costt is cost function: squared error loss or 0-1 loss or....
##numflds is #of folds in CV
##xnams is predictor list
##resp.coll is column in dataa with response
##knn.num = num of neighbors for knn

cost.fn <- function(r, pii = 0) mean(abs(r-pii) > 0.5)
##BEGIN CV FN
cv_loki.fn <- function(stsd,dataa, methodd, fmlaa, responsee, costt,numflds,xnamms, resp.coll, knn.numm){
               if (stsd == TRUE) {set.seed(1234)}
               nn <- dim(dataa)[1]
               kk <- floor(nn/numflds) ##size of each fold. May need to tweak
               fitt.pred <- rep(NA,nn)
#cat(nn,"\t",kk,"\t", mean(fitt.pred),"\n")

               fold <- list(NULL)
               fold[[1]] <- sample((1:nn),kk) ## sample to get first fold
               for (i in 2:(numflds-1)){ ##sample to get each  fold
                    fold[[i]] <- sample((1:nn)[-unlist(fold[1:(i-1)])],kk)
               }
               fold[[numflds]] <- (1:nn)[-unlist(fold[1:(numflds-1)])] ##last fold may not have size kk
 
              for (i in 1:length(fold)){ ##get held-out set, and training set
                  heldoutt <- fold[[i]]
                  trainn <- dataa[-heldoutt,] 
                  if (methodd == "ldaa"){ 
                      fitt <- lda(fmla,data = trainn,na.rm=T)
                      fp <- predict(fitt,dataa[heldoutt,])$posterior[,2]
                   }
                  else if (methodd=="logisticc"){
                       fitt <- glm(fmla,data = trainn, family=binomial,na.action=na.omit)
                       fp <- predict(fitt, dataa[heldoutt,], type = "response")
                   }
                   else if (methodd == "qdaa"){ 
                       fitt <- qda(fmla,data = trainn,na.rm=T)
                       fp <- predict(fitt,dataa[heldoutt,])$posterior[,2]
                   }
                   else if (methodd=="knnn"){
                        fitt <- knn(trainn[,xnamms],dataa[heldoutt,xnamms],trainn[,resp.coll],knn.numm)
                        fp <- as.numeric(fitt)-1
                   }
#cat(i, "\t",length(fp),"\t", length(heldoutt),"\t",dim(trainn),"\t", mean(fp),"\n")
#cat(round(fp,2))
                   fitt.pred[heldoutt]  <- fp
#cat(i, "\t",length(fp),"\t", length(heldoutt),"\t",dim(trainn),"\n")
#cat(mean(fp),"\t", mean(fitt.pred[!is.na(fitt.pred)]), "\t", sum(!is.na(fitt.pred)),"\n")
                   rm(heldoutt,trainn,fp)
              }
              errorr <- costt(responsee,fitt.pred)
#              acc.tbl <- table(dataa[,resp.coll],fitt.pred>0.5) 
#              acc <- sum(diag(acc.tbl))/sum(acc.tbl)
              outt <- list(numflds, methodd,fitt.pred, errorr)
              names(outt) <- c("number of folds", "fittting method", "predicted probabilities","error")
#cat(coef(fitt), "\n",length(fitt.pred), mean(fitt.pred),"\n")
#cat(outt[[4]],"\n",outt[[3]])
              outt
}

##END CV FN
##############################################################

###LOOCV misclassification error for Heart data for
##logistic
 jnkk <- cv_loki.fn(TRUE,Heart,"logisticc",fmla,as.numeric(Heart$AHD)-1,cost.fn,dim(Heart)[1])
jnkk[[4]]

#### LDA
 jnkk <- cv_loki.fn(TRUE, Heart,"ldaa",fmla,as.numeric(Heart$AHD)-1,cost.fn,dim(Heart)[1])
jnkk[[4]]



### QDA
 jnkk <- cv_loki.fn(TRUE, Heart,"qdaa",fmla,as.numeric(Heart$AHD)-1,cost.fn,dim(Heart)[1])
jnkk[[4]]

####k-NN with $k=5$
 jnkk <- cv_loki.fn(TRUE, Heart,"knnn",fmla,as.numeric(Heart$AHD)-1,cost.fn,dim(Heart)[1],xnams,15,5)
jnkk[[4]]


################
##FIVE-FOLD CV: one iteration for logistic
#logr.fit <- glm(fmla,data = Heart, family=binomial)
cv.glm(Heart,logr.fit,cost.fn,K=5)$delta


##TEN-FOLD: one iteration for logistic
cv.glm(Heart,logr.fit,cost.fn,K=10)$delta

###########################################################
###FIVE-FOLD cross validation misclassification error for Heart data for
##logistic, lda, qda, kNN

###Logistic
 jnkk <- cv_loki.fn(TRUE, Heart,"logisticc",fmla,as.numeric(Heart$AHD)-1,cost.fn,5)
jnkk[[4]]

##### LDA
 jnkk <- cv_loki.fn(TRUE,Heart,"ldaa",fmla,as.numeric(Heart$AHD)-1,cost.fn,5)
jnkk[[4]]

#### QDA
 jnkk <- cv_loki.fn(TRUE,Heart,"qdaa",fmla,as.numeric(Heart$AHD)-1,cost.fn,5)
jnkk[[4]]

######k-NN with $k=5$
 jnkk <- cv_loki.fn(TRUE,Heart,"knnn",fmla,as.numeric(Heart$AHD)-1,cost.fn,5,xnams,15,5)
jnkk[[4]]


################################################################
###REPEATED CV

### doing repeated CV
numcvs <- 100

lgg <- rep(-1,numcvs) #logistic
ldd <- rep(-1,numcvs) #lda
qdd <- rep(-1,numcvs) #qda
kkn <- rep(-1,numcvs) #knn
stsdd <- FALSE

for (j in 1:numcvs){
         cat(j)
         lgg[j] <- cv_loki.fn(stsdd,Heart,"logisticc",fmla,as.numeric(Heart$AHD)-1,cost.fn,5,FALSE)[[4]]
         ldd[j] <- cv_loki.fn(stsdd,Heart,"ldaa",fmla,as.numeric(Heart$AHD)-1,cost.fn,5,FALSE)[[4]]
         qdd[j] <- cv_loki.fn(stsdd,Heart,"qdaa",fmla,as.numeric(Heart$AHD)-1,cost.fn,5,0,15)[[4]]
         kkn[j] <- cv_loki.fn(stsdd,Heart,"knnn",fmla,as.numeric(Heart$AHD)-1,cost.fn,5,xnams,15,5)[[4]]
}

repcv.df <-data.frame(c(rep("logistic",numcvs),rep("lda",numcvs),rep("qda",numcvs),rep("knn",numcvs)), c(lgg,ldd,qdd,kkn),row.names=NULL)

names(repcv.df) <- c("Method", "Error")

boxplot(repcv.df$Error~repcv.df$Method)
######################
##BOOTSTRAP 
##bootstrap sample: 
nn <- dim(Heart)[1]
boott1 <- sample(nn,nn,replace=T)
sort(boott1)
table(boott1)


###Using bootstrap for estimation
##Recall the logistic regression model for AHD.
##model fit coefs
summary(glm(fmla,data = Heart, family = binomial))

#################
##bootstrap coefs: 1 iteration
nn <- dim(Heart)[1]
coef(glm(fmla,data = Heart, family = binomial, subset=sample(nn,nn,replace=T)))[c("RestBP","MaxHR")]

####
##repeat 10 times
nmbts <- 10
btcoefs <- matrix(0, ncol=2, nrow=nmbts)
for (i in 1:nmbts){
                 btcoefs[i,] <- coef(glm(fmla,data = Heart, family = binomial, subset=sample(nn,nn,replace=T)))[c("RestBP","MaxHR")]
	       }

btcoefs ##note coefs vary over bootstraps

##compare bootstrap SEs to aymptotic SEs
sqrt(diag(var(btcoefs)))

####repeat with nmbts <- 1000
nmbts <- 1000
btcoefs <- matrix(0, ncol=2, nrow=nmbts)
for (i in 1:nmbts){
                   btcoefs[i,] <- coef(glm(fmla,data = Heart, family = binomial, subset=sample(nn,nn,replace=T)))[c("RestBP","MaxHR")]
	          }
summary(btcoefs)

##compare bootstrap SEs to aymptotic SEs
sqrt(diag(var(btcoefs)))

##corr
cor(btcoefs[,1],btcoefs[,2])

###silly example
##are coefs of MaxHR and RestBP equal  (in abs value)?
summary(abs(btcoefs[,2]) - abs(btcoefs[,1]))
quantile(abs(btcoefs[,2]) - abs(btcoefs[,1]), c(0.025,0.975))

##########################################################
##Trouble shooting why LOOCV for knn gives different error estimates
ajnk <- cv_loki.fn(FALSE, Heart,"knnn",fmla,as.numeric(Heart$AHD)-1,cost.fn,dim(Heart)[1],xnams,15,5)

bjnk <- cv_loki.fn(FALSE, Heart,"knnn",fmla,as.numeric(Heart$AHD)-1,cost.fn,dim(Heart)[1],xnams,15,5)

#ajnk[[4]]
#bjnk[[4]]
#ajnk[[3]] - bjnk[[3]] ##subjects 12, 31, 233 have discrepant predictions
##xnams <- c("Age", "Sex", "RestBP", "Chol", "MaxHR")
rownn <- 233 #12 
kknum <- 5

##Note: in LOOCV the  neighbor with 0 distance is the datum itself, which is dropped.
sort(
  (Heart$Age-Heart$Age[rownn])^2 + 
  (Heart$Sex -Heart$Sex[rownn])^2+ 
  (Heart$RestBP - Heart$RestBP[rownn])^2 + 
  (Heart$Chol - Heart$Chol[rownn])^2 + 
  (Heart$MaxHR - Heart$MaxHR[rownn])^2
)[1:11]#[1:kknum]

#Nearest distances for row 12, kknum=5
####  0 140 262 268 306 307 307 342 367 387 393
###drop 0, amiguity between nbrs 5 and 6 both have 307 dist

##row 31 distances nbrs 5 and 6 have distance 251
##  0  99 162 183 239 251 251 307 358 366 377

##row 233 distances nbrs 5 and 6 have distance 753
##  0 318  522  613  653  753  753  834  860 1240 1681
