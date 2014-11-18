fun <- function( j ) {
 library(mvtnorm);
 library(glmnet);
 n=1100;
 p=1000;
 q=100;
 beta<-rep(0,min(p,n));
 index<-1; #index of PCs associatd with outcome Y
 beta[index]<-rep(0.01,length(index));
 alpha=0; #0:ridge penalty,1:lasso#
 alpha2=0;#0:ridge penalty,1:lasso#
 n.train=50;
 nfolds=5;
 nperm=500;
 index.val<-(n.train+1):n;

 source("SimAugmPred_Functions.R");

#Correlation matrix of X1#
rho.Hub = c(.9,.5);
tau.Hub = c( (.9 - .5) / (10-2), (.7-.6)/(5-2));
eps.Hub = min(1-(rho.Hub) - 0.75*(tau.Hub) ) - .01;
cor1 <- simcor.H(k=4, size=rep(250,4),	rho=rbind( c(.9,.8),c(.5,.3),c(.9,.6),c(.05,.03) ), power=1,
	epsilon=eps.Hub, eidim=2);

#Correlation matrix of X2#
cor2<-diag(q);


#cor2 <- simcor.H(k=4, size=rep(25,4),	rho=rbind( c(.9,.8),c(.5,.3),c(.9,.6),c(.5,.3) ), power=1,
#	epsilon=eps.Hub, eidim=2);



#generate X1, X2 and Y
data<-Sim.Model.1(n=n,p=p,q=q,cor.mat.1=cor1,cor.mat.2=cor2,beta=beta);

#Create outer CV partition#
folds<-createFolds(1:n.train, k = nfolds, list = T)

# Estimates with 2CV#
index.train<-1:n.train;
Y.X1.2CV<-glmnet.2CV(X=data$X1[index.train,],Y=data$Y[index.train],alpha=alpha,folds=folds,nfolds=nfolds);
Q2.X1.2CV<-Y.X1.2CV$Q2;
p.X1.2CV<-Y.X1.2CV$p;

res<-data$Y[index.train]-p.X1.2CV;

X2.2CV<-Q2.test.2CV(X1=data$X1[index.train,],X2=data$X2[index.train,],res=res,Y=data$Y[index.train],alpha=alpha2,folds=folds,nfolds=nfolds,nperm=nperm);
Q2.X2.2CV<-X2.2CV$Q2
Q2.pvalue<-X2.2CV$pvalue1
Q2.GLOBAL<-X2.2CV$Q2.global


#Validation quantities#
X1.scale<-sapply(1:ncol(data$X1),function(i)scale(data$X1[index.train,i]))
fit1=glmnet(X1.scale,data$Y[index.train],family="gaussian",standardize=F,alpha=alpha);
cv.fit1=cv.glmnet(X1.scale,data$Y[index.train],family="gaussian",standardize=F,alpha=alpha);

X1.VAL.scale<-sapply(1:ncol(data$X1),function(i)scale(data$X1[index.val,i]));
p.cv.glmnet=predict(fit1,X1.VAL.scale,s=cv.fit1$lambda.min);
Q2.X1.VAL<-1-(sum((data$Y[index.val]-p.cv.glmnet)^2)/sum((data$Y[index.val]-mean(data$Y[index.val]))^2));

res0<-data$Y[index.val]-p.cv.glmnet;

X2.val1.scale<-sapply(1:ncol(data$X2),function(i)scale(data$X2[index.val[1:n.train],i]))
fit2=glmnet(X2.val1.scale,res0[1:n.train],family="gaussian",standardize=F,alpha=alpha);
cv.fit2=cv.glmnet(X2.val1.scale,res0[1:n.train],family="gaussian",standardize=F,alpha=alpha);

X2.val2.scale<-sapply(1:ncol(data$X2),function(i)scale(data$X2[(2*n.train+1):n,i]));
p.cv.glmnet=predict(fit2,X2.val2.scale,s=cv.fit2$lambda.min);
Q2.X2.VAL<-1-(sum((res0[(n.train+1):(n-n.train)]-p.cv.glmnet)^2)/sum((res0[(n.train+1):(n-n.train)]-mean(res0[(n.train+1):(n-n.train)]))^2));

Q2.GLOBAL.VAL<-1-(sum((res0[(n.train+1):(n-n.train)]-p.cv.glmnet)^2)/sum((data$Y[(2*n.train+1):n]-mean(data$Y[(2*n.train+1):n]))^2));



#FINAL OUTPUT FOR EACH SCENARIO:

out<-c(Q2.X1.VAL,Q2.X1.2CV,Q2.X2.VAL,Q2.X2.2CV,Q2.GLOBAL.VAL,Q2.GLOBAL,Q2.pvalue);
return(out);
}