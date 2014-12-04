fun <- function( j ) {
 n=100;
 p=1000;
 q=100;
 beta<-rep(0,min(p,n));
 index<-1; #index of PCs associatd with outcome Y
 beta[index]<-rep(0.01,length(index));
 alpha=0; #0:ridge penalty,1:lasso#
 alpha2=0.5;#0:ridge penalty,1:lasso#
 n.train=50;
 nfolds=5;
 nperm=500;
 index.val<-(n.train+1):n;

 source("SimAugmPred_Functionsv2.R");

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
seed=j	
set.seed(j)
data<-Sim.Model.1(n=n,p=p,q=q,cor.mat.1=cor1,cor.mat.2=cor2,beta=beta, seed);


#Create outer CV partition#
folds<-createFolds(1:n.train, k = nfolds, list = T)

# Estimates with 2CV#
index.train<-1:n.train;
Y.X1.2CV <- glmnet.2CV(X=data$X1[index.train,],Y=data$Y[index.train],alpha=alpha,folds=folds,nfolds=nfolds)
Q2.X1.2CV <- Y.X1.2CV$Q2
p.X1.2CV <- Y.X1.2CV$p
bhat.X1.2CV <- Y.X1.2CV$betahat[-1,] #added coefficients extraction
RIscoreu0v0 <- rowSums(matrix(as.logical(bhat.X1.2CV), ncol=ncol(bhat.X1.2CV), nrow=nrow(bhat.X1.2CV) ))

res<-data$Y[index.train]-p.X1.2CV;

X2.2CV<-Q2.test.2CV(X1=data$X1[index.train,],X2=data$X2[index.train,],res=res,Y=data$Y[index.train],alpha=alpha2,folds=folds,nfolds=nfolds,nperm=nperm);
alarm()
Q2.X2.2CV<-X2.2CV$Q2
Q2.pvalue<-X2.2CV$pvalue1
Q2.GLOBAL<-X2.2CV$Q2.global

#Validation quantities#
X1.scale= scale(data$X1) 
fit1=glmnet(X1.scale[index.train,], data$Y[index.train],family="gaussian",standardize=F,alpha=alpha);
cv.fit1=cv.glmnet(X1.scale[index.train,],data$Y[index.train],family="gaussian",standardize=F,alpha=alpha);

p.cv.glmnet=predict(fit1,X1.scale[index.val,],s=cv.fit1$lambda.min);
bhat.cv.glmnet=coef(fit1,s=cv.fit1$lambda.min)[-1];
Q2.X1.VAL<-1-(sum((data$Y[index.val]-p.cv.glmnet)^2)/sum((data$Y[index.val]-mean(data$Y[index.val]))^2));

res0<-data$Y[index.val]-p.cv.glmnet;

X2.scale<- scale(data$X2)
fit2=glmnet(X2.scale[index.val[1:n.train], ], res0[1:n.train],family="gaussian",standardize=F,alpha=alpha);
cv.fit2=cv.glmnet(X2.scale[index.val[1:n.train], ], res0[1:n.train],family="gaussian",standardize=F,alpha=alpha);
b2hat.cv.glmnet=coef(fit2,s=cv.fit2$lambda.min)[-1];
#p.cv.glmnet=predict(fit2,X2.scale[(2*n.train+1):n,],s=cv.fit2$lambda.min);
p.cv.glmnet=predict(fit2,X2.scale[(n.train+1):n,],s=cv.fit2$lambda.min);

Q2.X2.VAL <- 1-(sum((res0-p.cv.glmnet)^2)/sum((res0-mean(res0))^2));
Q2.GLOBAL.VAL <- 1-(sum((res0-p.cv.glmnet)^2)/sum((data$Y[(n.train+1):n]-mean(data$Y[(n.train+1):n]))^2));

#FINAL OUTPUT FOR EACH SCENARIO:

out<-c(Q2.X1.VAL,Q2.X1.2CV,Q2.X2.VAL,Q2.X2.2CV,Q2.GLOBAL.VAL,Q2.GLOBAL,Q2.pvalue);
return(out);
}