#Project:MIMOmics              #
#                              #
#Author: Mar Rodríguez Girondo #
#Date: November 2014           #
################################

###########################################################################
##  Simulating the Hub Matrix (entries filled in using Toeplitz structure)#
###########################################################################

# this function calculates a Toeplitz matrix with values descending
# from a user specified maximum to minimum.  The matrix has a 
# block diagonal structure.  The size and base correlation for each
# block is user specified.


# k is the number of groups
# size is a vector of length k specifying the size of each group 
# rho is a vector of length k specifying base correlation values
# epsilon <- (1-min(rho) - 0.75*min(tau) ) - .01
# tau_k = (max(rho_k) - min(rho_k) )/ (size_k -2) 
# eidim is the space from which the noise is generated, the smaller the more noise
# power = 2 makes the correlations stay high
# power = 0.5 makes the correlations descent rapidly


simcor.H <- function(k=6, size=c(10,5,8,7,15,50), 
	rho=rbind(c(.9,.7), c(.7,.7), c(.7,.2), c(.5,.3), c(.9,.85), c(.3,.2)), power=1,
	epsilon=.08, eidim=2){


	ndim <- sum(size)# dim of correlation matrix
	bigcor<- matrix(rep(0, ndim*ndim), ncol=ndim)

### generating the basic correlation matrix


	for (i in 1:(k) ){

	cor <- toeplitz(rho.func(rho[i,1],rho[i,2],power,size[i]) )

	if (i==1){bigcor[1:size[1], 1:size[1]] <- cor}
	if (i!=1){bigcor[(sum(size[1:(i-1)]) + 1):sum(size[1:i]),
		(sum(size[1:(i-1)]) + 1):sum(size[1:i])] <- cor}
	}
	diag(bigcor) <- 1 - epsilon


### adding noise to the correlation matrix

	eivect <- c( )
	for (i in 1:ndim) {
	ei <- runif(eidim, -1, 1)
	eivect <- cbind(eivect, sqrt(epsilon) * ei / sqrt(sum(ei^2) ) )
	}
 

	bigE <- t(eivect) %*% eivect
	cor.nz <- bigcor + bigE
	cor.nz
	}



# rho.func is needed for filling in the rest of the structure of the Hub
# correlation matrix
# r.max is the maximum user specified correlation
# r.min is the minimum user specified correlation
# power is the power at which the correlations descend
# p is the size of the correlation block

rho.func <- function(r.max, r.min, power,p){
	rhovec <-c()

	rhovec[1] <- 1
	for(i in 2:p){
	rhovec[i] <- r.max - ((i-2)/(p-2))^power*(r.max-r.min)
	}
	rhovec}





#Simulation of data #
Sim.Model.1<-function(n,p,q,cor.mat.1,cor.mat.2,beta){

#Generate X1 (first source of predictors)
X1<-rmvnorm(n=n,sigma=cor.mat.1)
s1<-svd(X1)
U01<-rmvnorm(n=n,sigma=diag(min(p,n)))
X1<-U01%*%diag(s1$d)%*%t(s1$v)


#Generate X2 (second source of predictors)
X02<-rmvnorm(n,sigma=cor.mat.2)
s2<-svd(X02)
U02<-rmvnorm(n,sigma=diag(min(q,n)))
U02[,1:min(p,q)]<-U01[,1:min(p,q)]
#U02[,index]<-rmvnorm(n,sigma=diag(length(index)))


V02<-s2$v
D02<-s2$d
X2<-U02%*%diag(s2$d)%*%t(V02)

#Generate response
betaX1<-s1$v%*%beta
Y<-X1%*%betaX1+rnorm(n,0,1)

return(list(X1=X1,X2=X2,Y=Y,betaX1=betaX1))
}

#Create folds for K-fold CV procedures#
createFolds<-function (y, k = 10, list = TRUE, returnTrain = FALSE) 
{
    if (is.numeric(y)) {
        cuts <- floor(length(y)/k)
        if (cuts < 2) 
            cuts <- 2
        if (cuts > 5) 
            cuts <- 5
        y <- cut(y, unique(quantile(y, probs = seq(0, 1, length = cuts))), 
            include.lowest = TRUE)
    }
    if (k < length(y)) {
        y <- factor(as.character(y))
        numInClass <- table(y)
        foldVector <- vector(mode = "integer", length(y))
        for (i in 1:length(numInClass)) {
            seqVector <- rep(1:k, numInClass[i]%/%k)
            if (numInClass[i]%%k > 0) 
                seqVector <- c(seqVector, sample(1:k, numInClass[i]%%k))
            foldVector[which(y == dimnames(numInClass)$y[i])] <- sample(seqVector)
        }
    }
    else foldVector <- seq(along = y)
    if (list) {
        out <- split(seq(along = y), foldVector)
        names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))), 
            sep = "")
        if (returnTrain) 
            out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
    }
    else out <- foldVector
    out
}


############################################################################
#glmnet wrapper to get double cv predictions and predictive performance Q^2#   
#number of folds (nfolds) and partitions (folds) are entered as input      #
############################################################################

glmnet.2CV<-function(X,Y,alpha=alpha,folds,nfolds){

if(nrow(X)!=length(Y)) {stop("Dimensions of X and Y do not match!")}
alpha=alpha
n=nrow(X)


fit.glmnet=lapply(1:nfolds,function(i)glmnet(X[-folds[[i]],],Y[-folds[[i]]],family="gaussian",standardize=T,alpha=alpha))
cv.fit.glmnet=lapply(1:nfolds,function(i)cv.glmnet(X[-folds[[i]],],Y[-folds[[i]]],family="gaussian",standardize=T,alpha=alpha))

p.cv.glmnet=unlist(lapply(1:nfolds,function(i)predict(fit.glmnet[[i]],matrix(X[folds[[i]],],ncol=ncol(X)),s=cv.fit.glmnet[[i]]$lambda.min)))
p.cv.glmnet<-p.cv.glmnet[order(unlist(folds))]

#Calculate CV mean of the outcome#

p0<-rep(NA,length(Y))
glm0<-lapply(1:nfolds,function(i)glm(Y[-folds[[i]]]~1,family="gaussian"))
for(i in 1:nfolds){
p0[folds[[i]]]=rep(coef(glm0[[i]]))}


Q2.glmnet<-1-(sum((Y-p.cv.glmnet)^2)/sum((Y-p0)^2))

return(list(Q2=Q2.glmnet,p=p.cv.glmnet))
}


############################################################################
#Added value testing                                                       #
############################################################################

Q2.test.2CV<-function(X1,X2,Y,res,alpha=alpha,folds,nfolds,nperm=100){

if(nrow(X1)!=length(Y)) {stop("Dimensions of X and Y do not match!")}
if(nrow(X2)!=length(Y)) {stop("Dimensions of X and Y do not match!")}


tmp<-glmnet.2CV(X2,res,alpha=alpha,folds=folds,nfolds=nfolds)
p=tmp$p
Q2=tmp$Q2


p0<-rep(NA,length(Y))
glm0<-lapply(1:nfolds,function(i)glm(Y[-folds[[i]]]~1,family="gaussian"))
for(i in 1:nfolds){
p0[folds[[i]]]=rep(coef(glm0[[i]]))}

Q2.global<-1-(sum((res-p)^2)/sum((Y-p0)^2))


p1.null<-sapply(1:nperm,function(i)glmnet.2CV(X1,Y,alpha=alpha,folds=folds,nfolds=nfolds)$p)
Q2.null<-sapply(1:nperm,function(i)glmnet.2CV(X2[sample(nrow(X2)),],Y-p1.null[,i],alpha=alpha,folds=folds,nfolds=nfolds)$Q2)
pvalue1 <- sum(Q2.null>Q2)/nperm 
pvalue2 <- sum(abs(Q2.null)>abs(Q2))/nperm 

return(list(Q2=Q2,pvalue1=pvalue1,pvalue2=pvalue2,Q2.global=Q2.global))
}

