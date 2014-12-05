source("SimAugmPred_Functionsv2.R")

fun <- function(j) {
seed = j	
set.seed(j)
 n = 100
 p = 1000
 q = 100
 beta <- rep(0,min(p,n))
 index <- 1 #index of PCs associatd with outcome Y
 beta[index] <- rep(0.01,length(index))
 alphas = seq(0, 1, by = 0.1)
 n.train = 50
 nfolds = 5
 nperm = 100
 index.val<-(n.train+1):n

#Correlation matrix of X1#
rho.Hub = c(.9,.5)
tau.Hub = c( (.9 - .5) / (10-2), (.7-.6)/(5-2))
eps.Hub = min(1-(rho.Hub) - 0.75*(tau.Hub) ) - .01
cor1 <- simcor.H(k=4, size=rep(250,4), rho=rbind( c(.9,.8),c(.5,.3),c(.9,.6),c(.05,.03) ), power=1, epsilon=eps.Hub, eidim=2)

#Correlation matrix of X2#
cor2<-diag(q)

#cor2 <- simcor.H(k=4, size=rep(25,4),	rho=rbind( c(.9,.8),c(.5,.3),c(.9,.6),c(.5,.3) ), power=1,
#	epsilon=eps.Hub, eidim=2);

#generate X1, X2 and Y
data <- Sim.Model.1(n=n, p=p, q=q, cor.mat.1=cor1, cor.mat.2=cor2, beta=beta, seed)

#Create outer CV partition#
folds<-createFolds(1:n.train, k = nfolds, list = T)

# Estimates with 2CV#
index.train<-1:n.train
Y.X1.2CV <- glmnet.2CV(X=data$X1[index.train,],Y=data$Y[index.train],alphas=seq(0,1, by= 0.1),folds=folds,nfolds=nfolds)
Q2.X1.2CV <- Y.X1.2CV$Q2
p.X1.2CV <- Y.X1.2CV$p
bhat.X1.2CV <- Y.X1.2CV$betahat #added coefficients extraction
RIscoreu0v0 <- rowSums(matrix(as.logical(bhat.X1.2CV), ncol=ncol(bhat.X1.2CV), nrow=nrow(bhat.X1.2CV) ))
res<-data$Y[index.train]-p.X1.2CV
X2.2CV<-Q2.test.2CV(X1=data$X1[index.train,],X2=data$X2[index.train,],res=res,Y=data$Y[index.train],alpha=alpha2,folds=folds,nfolds=nfolds,nperm=nperm)

Q2.X2.2CV <- X2.2CV$Q2
Q2.pvalue <- X2.2CV$pvalue1
Q2.GLOBAL <- X2.2CV$Q2.global
betapValues <- X2.2CV$pvaluesBeta
###############################################
#         VALIDATION             SET          #
###############################################
#first data source:
X1.scale= scale(data$X1) 
cv2alpha=lapply(alphas, function(alpha) cv.glmnet(X1.scale[index.train,],data$Y[index.train],family="gaussian",standardize=FALSE,alpha=alpha, type.measure="mse"))
optimumPerAlpha<-sapply(seq_along(alphas), function(x){
				curcvs<-cv2alpha[[x]]
				indOfMin<-match(curcvs$lambda.min, curcvs$lambda)
				c(lam=curcvs$lambda.min, alph=alphas[x], cvu=curcvs$cvup[indOfMin],cvm=curcvs$cvm[indOfMin])})
posOfOptimum<-which.min(optimumPerAlpha["cvu",])
overall.lambda.min<-optimumPerAlpha["lam",posOfOptimum]
overall.alpha.min<-optimumPerAlpha["alph",posOfOptimum]
fit1=glmnet(X1.scale[index.train,], data$Y[index.train],family="gaussian",standardize=FALSE,alpha=overall.alpha.min, lambda=overall.lambda.min)
p.cv.glmnet=predict(fit1,X1.scale[index.val,])
bhat.cv.glmnet=coef(fit1)[-1]

Q2.X1.VAL<-1-(sum((data$Y[index.val]-p.cv.glmnet)^2)/sum((data$Y[index.val]-mean(data$Y[index.val]))^2))
res0<-data$Y[index.val]-p.cv.glmnet
#-------------------#
#second data source #
#-------------------#
X2.scale<- scale(data$X2)
cv2alpha2=lapply(alphas, function(alpha) cv.glmnet(X2.scale[index.train,],data$Y[index.train],family="gaussian",standardize=FALSE,alpha=alpha, type.measure="mse"))
optimumPerAlpha2<-sapply(seq_along(alphas), function(x){
			curcvs<-cv2alpha2[[x]]
			indOfMin<-match(curcvs$lambda.min, curcvs$lambda)
			c(lam=curcvs$lambda.min, alph=alphas[x], cvu=curcvs$cvup[indOfMin],cvm=curcvs$cvm[indOfMin])})
posOfOptimum2 <- which.min(optimumPerAlpha2["cvu",])
overall.lambda.min2 <- optimumPerAlpha2["lam",posOfOptimum2]
overall.alpha.min2 <-optimumPerAlpha2["alph",posOfOptimum2]
fit2=glmnet(X2.scale[index.train,], data$Y[index.train],family="gaussian",standardize=FALSE,alpha=overall.alpha.min2, lambda=overall.lambda.min2)
b2hat.cv.glmnet2=coef(fit2)[-1]
p.cv.glmnet2 = predict(fit2, X2.scale[index.val,])

Q2.X2.VAL <- 1-(sum((res0-p.cv.glmnet2)^2)/sum((res0-mean(res0))^2))
Q2.GLOBAL.VAL <- 1-(sum((res0-p.cv.glmnet2)^2)/sum((data$Y[index.val]-mean(data$Y[index.val]))^2))

#FINAL OUTPUT FOR EACH SCENARIO:
out<-list(statistics = c(Q2.X1.VAL, Q2.X1.2CV, Q2.X2.VAL, Q2.X2.2CV, Q2.GLOBAL.VAL, Q2.GLOBAL, Q2.pvalue),
		 variables = list(bhatX1=bhat.X1.2CV, RIvarX1=RIscoreu0v0, PvalX2var=betapValues, bhatX1val= bhat.cv.glmnet, bhatX2val=b2hat.cv.glmnet))
out
}