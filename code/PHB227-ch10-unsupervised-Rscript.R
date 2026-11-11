setwd("~/Desktop/myfiles/Teaching/PHB227/PHB227-2024Fall")
getwd()
library(ISLR2)
library(knitr)
#library(MASS)
library(class)
library(glmnet)
library(methods)
library(tree)
library(randomForest)
#library(ggplot2)

options(width=60)  # make the printing fit on the page
set.seed(1121)   # make the results repeatable
#require(PASWR)
require(ggplot2)
require(xtable)



########################### PRINCIPAL COMPONENTS

####USARRESTS DATASET

USArrests[1:10,]

apply(USArrests,2,mean) 
apply(USArrests,2,var)  ##need to scale!

states=row.names(USArrests) #states

##FIT PCA 
pr.out=prcomp(USArrests, scale=TRUE)
names(pr.out)

pr.out

##EXTRACT PCs
pr.out$rotation  ##columns are eigenvectors 


##VISUALIZING
par(mfrow=c(1,1))
biplot(pr.out,scale=0,cex=0.5)
abline(h=0,v=0)

##PC SCORES
dim(pr.out$x)
pr.out$x
##by hand PCi score for jth row 
#sum(pr.out$rotation[,i]*(scale(USArrests)[j,]))

pr.out$x[c("Virginia","California","Georgia","Vermont","Connecticut"),]

pr.out$center
USArrests[c("Virginia","California","Georgia","Vermont","Connecticut"),]


####Variance explained by each principal component
pr.var=pr.out$sdev^2
pve <- pr.var/sum(pr.var)

pr.var
pve




########### HOW MANY PCs to KEEP??

par(mfrow=c(1,2))
plot(pve, xlab="Principal Component", ylab="Proportion of Variance Explained ", ylim=c(0,1),type="b")
plot(cumsum(pve), xlab="Principal Component ", ylab=" Cumulative Proportion of Variance Explained", ylim=c(0,1), type = "b")



##### EFFECT OF (NOT) SCALING
par(mfrow=c(1,1))
pr.out2=prcomp(USArrests, scale=FALSE)
biplot(pr.out2,scale=0,xlabs=rep(".", nrow(USArrests)))

pr.out2$rotation
pr.out$scale



###OTHER PCs
biplot(pr.out,scale=0,choices=c(1,3),cex=0.5)
abline(h=0,v=0)

pr.out$center
pr.out$x[c("Alaska","Louisiana","Rhode Island"),c(1,3)]
USArrests[c("Alaska","Louisiana","Rhode Island"),]


#####SIGN FLIP
par(mfrow=c(1,1))
pr.out3 <- pr.out
pr.out3$rotation <- -pr.out3$rotation
pr.out3$x <- -pr.out3$x
biplot(pr.out3,scale=0,cex=0.5)
abline(h=0,v=0)

#########################################################
############### K MEANS CLUSTERING 

##DATASET: GENE EXPRESSION IN CANCER
nci.labs <-  NCI60$labs ##64 samples' cancer types
nci.data <-  NCI60$data ##64 samples with 6830 gene expressions

nci.labs
dim(nci.data)
table(nci.labs)

##scale variables
sd.data=scale(nci.data) ##scale 
set.seed (2024)

 
##FIT k-MEANS 
km.out <- kmeans(sd.data, 4, nstart=20) ##kmeans with 4 clusters
names(km.out)

##cluster centroids
dim(km.out$centers)

km.clusters <- km.out$cluster
km.clusters

table(km.clusters)
km.out$size


table(nci.labs,km.clusters)

c(km.out$totss,km.out$tot.withinss,km.out$betweenss,km.out$tot.withinss+km.out$betweenss)
km.out$withinss
sum(km.out$withinss)

####### K-means clustering: random restarts}

###random starts
set.seed(2024)
km.out1 <- kmeans(sd.data, 4, nstart=1)
km.out1$tot.withinss

km.out5 <- kmeans(sd.data, 4, nstart=5)
km.out5$tot.withinss

km.out20 <- kmeans(sd.data, 4, nstart=20)
km.out20$tot.withinss

km.out100 <- kmeans(sd.data, 4, nstart=100)
km.out100$tot.withinss
###How will the total sum of squares vary with random starts?
#########################

###CHOOSING K in K-Means
##within ss will decrease to 0, so look at plot and find elbow, but there may not always be a clear elbow

##choosing $K$
k.withinss <- numeric(0)
k.betwss <- numeric(0)
chh <- numeric(0)##Calinski-Harabasz
for (i in seq(2,10,1)) {
    cat(i,"\t")
    tmpp <- kmeans(sd.data, i, nstart=20)
    k.withinss[i] <- tmpp$tot.withinss 
    k.betwss[i] <- tmpp$betweenss
    tmpd <- tmpp$tot.withinss/(dim(sd.data)[1]-i)
    tmpn <- tmpp$betweenss/(i-1)
#cat(i,"\t",tmpn,"\t",dim(sd.data)[1]-i,"\t", tmpd,"\n")
    chh[i]<-  tmpn/tmpd
    rm(tmpp,tmpn,tmpd)
}

plot(seq(2,10,1), k.withinss[-1],xlab="Number of clusters", ylab="Within cluster variation")
lines(seq(2,10,1), k.withinss[-1])


plot(seq(2,10,1), chh[-1],xlab="Number of clusters", ylab="Calinksi-Harabasz Index")
lines(seq(2,10,1), chh[-1])



###################### nci gene data: 3 clusters
###3-cluster soln
km.out3cls <- kmeans(sd.data, 3, nstart=20)

km.out3cls$size

table(nci.labs, km.out3cls$cluster)
table(km.out$cluster, km.out3cls$cluster)

###############################################################
##### HIERARCHICAL CLUSTERING
nci.labs <-  NCI60$labs ##64 samples' cancer types
nci.data <-  NCI60$data ##64 samples with 6830 gene expressions


##scale variables
sd.data=scale(nci.data) ##scale 
data.dist=dist(sd.data) ##pairwise distances between the samples

##https://www.rdocumentation.org/packages/stats/versions/3.4.1/topics/hclust
hc.complete <- hclust(data.dist) ##hierarchical clustering with complete linkage


names(hc.complete)

hc.complete$method
hc.complete$dist.method


##plot dendrogram
plot(hclust(data.dist), labels=nci.labs, main="Complete
Linkage", xlab="", sub="",ylab="",cex=0.5)



###deconstructing output
#hc.complete$height ##height at which clusters fuse

#hc.complete$merge  ##which clusters fuse; -j refers to observation j; +j refers to cluster at step j

cbind(hc.complete$merge,hc.complete$height) 
nci.labs


### create clusters by cutting tree
hc.cplt.cls <- cutree(hc.complete,4)
hc.cplt.cls
table(hc.cplt.cls)
t(table(hc.cplt.cls, nci.labs))

#table(cutree(hc.complete,5))
#######Average linkage
plot(hclust(data.dist, method="average"), labels=nci.labs,
main="Average Linkage", xlab="", sub="",ylab="",cex=0.5)

##compare complete to average linkage 
hc.average <- hclust(data.dist, method="average")
hc.avg.cls <- cutree(hc.average,4)

table(hc.avg.cls,hc.cplt.cls)


##compare complete to kmeans
table(km.clusters,hc.cplt.cls)

########## Using distance/height to define clusters
hc.cls.ht139 <- cutree(hc.complete,h=139)
table(hc.cls.ht139,hc.cplt.cls)

hc.cls.ht120 <- cutree(hc.complete,h=120)
table(hc.cls.ht120)

##plot dendogram 
plot(hclust(data.dist), labels=nci.labs, main="Complete
Linkage", xlab="", sub="",ylab="",cex=0.5)

abline(h=139)
abline(h=120)
