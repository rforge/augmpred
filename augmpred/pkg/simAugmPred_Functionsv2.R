#Project:MIMOmics              #
#                              #
#Author: Mar Rodríguez Girondo  with additions of Tatsiana Khamiakova#
#Date: November 2014           #
################################
source("includeCcode.R")
library(glmnet)
###########################################################################
##  Simulating the Hub Matrix (entries filled in using Toeplitz structure)#
###########################################################################
# latets addition: instead of mvnorm function from R
#use C source function to generate multivariate normals
#syntax is the same, but it works better on higher dimensions
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
#seed parameter is added to pass to the C function
Sim.Model.1<-function(n, p,q,cor.mat.1,cor.mat.2,beta, seed){

#Generate X1 (first source of predictors)

X1<-mvrnormArma(n,mu=rep(0, dim(cor.mat.1)[1]), sigma=cor.mat.1, seed)
s1<-svd(X1)
U01<- mvrnormArma(n,mu=rep(0, min(p,n)), sigma=diag(min(p,n)), seed)

X1<-U01%*%diag(s1$d)%*%t(s1$v)

#Generate X2 (second source of predictors)
X02<- mvrnormArma(n,mu=rep(0, dim(cor.mat.2)[1]), sigma=cor.mat.2, seed)
s2<-svd(X02)
U02<- mvrnormArma(n,mu=rep(0, min(q,n)), sigma=diag(min(q,n)), seed)
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
#Tatsiana: variable selection 
#function to calculate relative importance of the features
#should be called only if it is desired, in functions
#glmnet.2CV and Q2.test.2CV

RelativeImportance <- function(betas, AccVec, u=0, v=0){
	if(v==0){weights <- matrix(as.numeric(as.logical(betas)), ncol=ncol(betas), nrow=nrow(betas) ) } 
	else if (v==1) {weights <- abs(betas)} else if (v==2){weights <- betas^2}
	RI <- AccVec^u%*% t(weights)
	return(list(RI=as.vector(RI), N= length(AccVec)))
}
RelativeImportanceu0v0 <- function(betas){
	
	rowSums(matrix(as.logical(betas), ncol=ncol(betas), nrow=nrow(betas) ))
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
#a function for cross-validation on alpha and lambda
#addition for the optimum alpha search
#see stackexchange: http://stats.stackexchange.com/questions/17609/cross-validation-with-two-parameters-elastic-net-case/17612#17612
cv2.glmnet <- function(x, y, weights, offset = NULL, lambda = NULL, type.measure = c("mse", 
				"deviance", "class", "auc", "mae"), nfolds = 10, foldid, grouped = TRUE, keep = FALSE, parallel = FALSE, alphas=seq(0,1,by=0.1),...)
{
	#preamble: checking for parameters and passing them from ... to the functioin
	if (missing(type.measure)) 
		type.measure = "default"
	else type.measure = match.arg(type.measure)
	if (!is.null(lambda) && length(lambda) < 2) 
		stop("Need more than one value of lambda for cv.glmnet")
	N = nrow(x)
	if (missing(weights)) 
		weights = rep(1, N)
	else weights = as.double(weights)
	y = drop(y)
	glmnet.call = match.call(expand.dots = TRUE)
	which = match(c("type.measure", "nfolds", "foldid", "grouped", 
					"keep", "alphas"), names(glmnet.call), F)
	if (any(which)) 
		glmnet.call = glmnet.call[-which]
	
	cv2alpha <- vector(mode="list", length=length(alphas))
	for(indAlpha in 1:length(alphas)){
	alpha=alphas[indAlpha]
	glmnet.call[[1]] = as.name("glmnet")
	glmnet.object = glmnet(x, y, weights = weights, offset = offset, 
			lambda = lambda, alpha=alpha, ...)
	#glmnet.object$call = glmnet.call #temporarily disable since it does not get alpha
	is.offset = glmnet.object$offset
	lambda = glmnet.object$lambda
	if (inherits(glmnet.object, "multnet")) {
		nz = predict(glmnet.object, type = "nonzero")
		nz = sapply(nz, function(x) sapply(x, length))
		nz = ceiling(apply(nz, 1, median))
	}
	else nz = sapply(predict(glmnet.object, type = "nonzero"), 
				length)
	if (missing(foldid)) 
		foldid = sample(rep(seq(nfolds), length = N))
	else nfolds = max(foldid)
	if (nfolds < 3) 
		stop("nfolds must be bigger than 3; nfolds=10 recommended")
	outlist = as.list(seq(nfolds))
	if (parallel && require(foreach)) {
		outlist = foreach(i = seq(nfolds), .packages = c("glmnet")) %dopar% 
				{
					which = foldid == i
					if (is.matrix(y)) 
						y_sub = y[!which, ]
					else y_sub = y[!which]
					if (is.offset) 
						offset_sub = as.matrix(offset)[!which, ]
					else offset_sub = NULL
					glmnet(x[!which, , drop = FALSE], y_sub, lambda = lambda, 
							offset = offset_sub, weights = weights[!which], alpha=alpha, 
							...)
				}
	}
	else {
		for (i in seq(nfolds)) {
			which = foldid == i
			if (is.matrix(y)) 
				y_sub = y[!which, ]
			else y_sub = y[!which]
			if (is.offset) 
				offset_sub = as.matrix(offset)[!which, ]
			else offset_sub = NULL
			outlist[[i]] = glmnet(x[!which, , drop = FALSE], alpha=alpha,
					y_sub, lambda = lambda, offset = offset_sub, 
					weights = weights[!which], ...)
		}
	}
	fun = paste("cv", class(glmnet.object)[[1]], sep = ".")
	cvstuff = do.call(fun, list(outlist, lambda, x, y, weights, 
					offset, foldid, type.measure, grouped, keep))
	cvm = cvstuff$cvm
	cvsd = cvstuff$cvsd
	cvname = cvstuff$name
	out = list(alpha=alpha, lambda = lambda, cvm = cvm, cvsd = cvsd, cvup = cvm + 
					cvsd, cvlo = cvm - cvsd, nzero = nz, name = cvname, glmnet.fit = glmnet.object)
	if (keep) 
		out = c(out, list(fit.preval = cvstuff$fit.preval, foldid = foldid))
	lamin = if (type.measure == "auc") 
				getmin(lambda, -cvm, cvsd)
			else getmin(lambda, cvm, cvsd)
	obj = c(out, as.list(lamin))
	class(obj) = "cv.glmnet"
	cv2alpha[[indAlpha]] = obj	
	}
#cv2alpha
#either leave it here and process in another function
#or put extra stuff to check the min (alpha, lambda)
#go for option 2, so that next to results for ecah alpha we have the min alpha and lambda
#collect the optimum lambda for each alpha
}

############################################################################
#glmnet wrapper to get double cv predictions and predictive performance Q^2#   
#number of folds (nfolds) and partitions (folds) are entered as input      #
############################################################################
glmnet.2CV<-function(X,Y,alphas,folds,nfolds){

if(nrow(X)!=length(Y)) {stop("Dimensions of X and Y do not match!")}
n=nrow(X)
fit.glmnet=lapply(1:nfolds,function(i) glmnet(X[-folds[[i]],],Y[-folds[[i]]],family="gaussian",standardize=T,alpha=alpha))

cv.fit.glmnet<- vector(mode="list", length=nfolds)
for(i in 1:nfolds){
			cv2alpha=lapply(alphas, function(alpha) cv.glmnet(X[-folds[[i]],],Y[-folds[[i]]],family="gaussian",standardize=T,alpha=alpha, type.measure="mse"))
			optimumPerAlpha<-sapply(seq_along(alphas), function(x){
						curcvs<-cv2alpha[[x]]
						indOfMin<-match(curcvs$lambda.min, curcvs$lambda)
						c(lam=curcvs$lambda.min, alph=alphas[x], cvu=curcvs$cvup[indOfMin],cvm=curcvs$cvm[indOfMin])
					})
#step 3: find the overall optimum
			posOfOptimum<-which.min(optimumPerAlpha["cvu",])
			overall.lambda.min<-optimumPerAlpha["lam",posOfOptimum]
			overall.alpha.min<-optimumPerAlpha["alph",posOfOptimum]
			overall.criterionthreshold<-optimumPerAlpha["cvu",posOfOptimum]
			finalmod = glmnet(X[-folds[[i]],],Y[-folds[[i]]], family="gaussian", standardize=T, alpha=overall.alpha.min, lambda=overall.lambda.min)
			cv.fit.glmnet[[i]] <- list(finalmod=finalmod, alphamin=overall.alpha.min, lambdamin=overall.lambda.min)
}

#cv2.fit.glmnet=lapply(1:nfolds,function(i)cv2.glmnet(X[-folds[[i]],],Y[-folds[[i]]],family="gaussian",standardize=T,alphasOfInterest=seq(0,1,by=0.1)))

#added: extracted beta-coefficients to monitor feature selection for Elastic net and Lasso
coef.glmnet = sapply(1:nfolds, function(i) as.matrix(coef(cv.fit.glmnet[[i]]$finalmod)))
# 
p.cv.glmnet=unlist(lapply(1:nfolds,function(i)predict(cv.fit.glmnet[[i]]$finalmod, matrix(X[folds[[i]],],ncol=ncol(X)))))
p.cv.glmnet<-p.cv.glmnet[order(unlist(folds))]

#Calculate CV mean of the outcome#

p0<-rep(NA,length(Y))
glm0<-lapply(1:nfolds,function(i)glm(Y[-folds[[i]]]~1,family="gaussian"))
for(i in 1:nfolds){
p0[folds[[i]]]=rep(coef(glm0[[i]]))}
Q2.glmnet<-1-(sum((Y-p.cv.glmnet)^2)/sum((Y-p0)^2))

return(list(Q2=Q2.glmnet,p=p.cv.glmnet, betahat=coef.glmnet))
}


############################################################################
#Added value testing                                                       #
#+ testing for the average number of selection counts on the permuted data
#which variables get selected by chance
############################################################################

Q2.test.2CV<-function(X1,X2,Y,res,alpha=alpha,folds,nfolds,nperm=100){

if(nrow(X1)!=length(Y)) {stop("Dimensions of X and Y do not match!")}
if(nrow(X2)!=length(Y)) {stop("Dimensions of X and Y do not match!")}


tmp<-glmnet.2CV(X2,res,alpha=alpha,folds=folds,nfolds=nfolds)
p=tmp$p
Q2=tmp$Q2
bhat= tmp$betahat[-1,]
RIscoreu0v0 <- rowSums(matrix(as.logical(bhat), ncol=ncol(bhat), nrow=nrow(bhat) ))

p0<-rep(NA,length(Y))
glm0<-lapply(1:nfolds,function(i)glm(Y[-folds[[i]]]~1,family="gaussian"))
for(i in 1:nfolds){
p0[folds[[i]]]=rep(coef(glm0[[i]]))}

Q2.global<-1-(sum((res-p)^2)/sum((Y-p0)^2))

tmpOutput <- lapply(1:nperm,function(i)glmnet.2CV(X1,Y,alpha=alpha,folds=folds,nfolds=nfolds))
p1.null <- sapply(1:nperm, function(i) tmpOutput[[i]]$p)
RIscoreu0v0.null  <- sapply(1:nperm, function(i) RelativeImportanceu0v0(tmpOutput[[i]]$betahat[-1,]))

#p1.null<-sapply(1:nperm,function(i)glmnet.2CV(X1,Y,alpha=alpha,folds=folds,nfolds=nfolds)$p)
#Q2.null<-sapply(1:nperm,function(i)glmnet.2CV(X2[sample(nrow(X2)),],Y-p1.null[,i],alpha=alpha,folds=folds,nfolds=nfolds)$Q2)
tmpOutput2 <- lapply(1:nperm,function(i)glmnet.2CV(X2[sample(nrow(X2)),],Y-p1.null[,i],alpha=alpha,folds=folds,nfolds=nfolds))
Q2.null <- sapply(1:nperm, function(i) tmpOutput2[[i]]$Q2)
RIscore2u0v0.null  <- sapply(1:nperm, function(i) RelativeImportanceu0v0(tmpOutput2[[i]]$betahat[-1,]))

pvalue1 <- sum(Q2.null>Q2)/nperm 
pvalue2 <- sum(abs(Q2.null)>abs(Q2))/nperm 
pvalueBetas2 <- rowSums(RIscore2u0v0.null>RIscoreu0v0)/nperm 
#getting RI score is computationally almost for free- we are there anyway

return(list(Q2=Q2,pvalue1=pvalue1,pvalue2=pvalue2,Q2.global=Q2.global,RI.X1=RIscoreu0v0.null, RI.X2=RIscore2u0v0.null, pvaluesBeta=pvalueBetas2))
}

