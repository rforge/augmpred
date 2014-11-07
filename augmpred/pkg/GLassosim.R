# Project: Mimomics
# 
# Author: lucp2149
###############################################################################
SimulateNoise<- function(seed=123, NSamples=100, NGenes=1000, sdNoise=1){
		Sigma <- diag(NGenes)
		mu <- rep(0, NGenes)
		set.seed(seed)
		predictorsGrp <- mvrnormArma(NSamples, mu, Sigma, seed=seed)
		colnames(predictorsGrp) <- paste("gene", c(1:NGenes), sep="")
		outcome <- rnorm(NSamples, 0, sdNoise)
		SimData <- list(x=predictorsGrp, y=scale(outcome))
		return(SimData)
}


#make the simulation dataset following Whitten's model
#y=X*beta+eps
#X - matrix of predictors from a multivariate distribution N(0,Sigma)
#Sigma- block-diagonal matrix
SimulateData<- function(rho1 = 0.5, rho2 = 0.6, seed=123, NSamples=100, NGenes=1000, NGenesGroup=50,
			betasMargins=c(0.5, 0.7, -0.7, -0.5), sdNoise=1){
		Sigma <- diag(NGenes)
		Sigma[1:NGenesGroup, 1:NGenesGroup] <- rho1
		Sigma[(NGenesGroup+1):(NGenesGroup*2), (NGenesGroup+1):(NGenesGroup*2)]  <- rho2
		for(i in 1:(NGenesGroup*2)){
				Sigma[i,i] <- 1
		}
		set.seed(seed)
		mu <- rep(0, NGenes)
		predictorsGrp <- mvrnormArma(NSamples, mu, Sigma, seed=seed)
		colnames(predictorsGrp) <- paste("gene", c(1:NGenes), sep="")
		#simulate the outcome accoridng to the linear model 
		betas <- rep(0, NGenes)
		betas[1:floor(NGenesGroup/2)] <- runif(floor(NGenesGroup/2), betasMargins[1], betasMargins[2])
		betas[(NGenesGroup+1):(NGenesGroup+ floor(NGenesGroup/2))] <- runif(floor(NGenesGroup/2), betasMargins[3], betasMargins[4])
		outcome <- predictorsGrp %*% betas + rnorm(NSamples, 0, sdNoise)
		SimData <- list(x=predictorsGrp, y=scale(outcome))
		return(SimData)
}

FitElnetOrLassokFoldCV<-function(X,y,varnames,ind.train=1:nrow(X),ind.test=1:nrow(X),alpha.value=1,STDx=TRUE, penalty.factor){
	options(warn=-1)
	fit1 <- glmnet(X[ind.train,],y[ind.train],alpha = alpha.value,standardize = STDx, penalty.factor=penalty.factor)
	pred.fit1 <- predict(fit1,matrix(X[ind.test,],ncol=ncol(X)))#,s=(1:100)*0.01) - changed, since I would like to use the same lambdas as were chosen by glmnet to fit the model
	PRMSS<-colMeans((pred.fit1 - y[ind.test])^2)
	OptIndex<-which.min(PRMSS)
	MSEi<-min(PRMSS)
	CORR<-as.vector(cor(predict(fit1,matrix(X[ind.test,],ncol=ncol(X))),y[ind.test]))[OptIndex]
	#temp<-coef(fit1,s=OptIndex)[-1] # here we keep the intercept
	coefEstim<-coef(fit1)[-1,OptIndex]
	selVar<-varnames[which(abs(coefEstim) > 0)]
	finalPred<- pred.fit1[ , OptIndex][-1]
	return(list(selVar=selVar,MSEi=MSEi,CORR=CORR,coefEstim=coefEstim,finalPred=finalPred))
}

#cross-validation (default 3fold-cv)
CvELnetOrLasso<-function(y,x,alpha,STDx,ncv=10,nFold=3, penalty.factor)
{
	
	CORRcv<-MSEcv<-rep(NA,ncv)
	seliter <- numeric()
	noFeaSelected<-rep(NA,ncv)
	FeaSelected <- matrix(0, ncol=ncv, nrow=ncol(x))
	nTest<-floor(length(y)/nFold)
	nTrain<-length(y)-nTest
	
	for (j in 1:ncv){ 
		set.seed(j)
		ind.train <-as.vector(sort(sample(1:length(y),nTrain,replace=F) ) )
		ind.test <- as.vector(c(1:length(y))[-ind.train])
		res3fCv<-NA  
		res3fCv<-FitElnetOrLassokFoldCV(X=x,y=y,varnames=colnames(x),ind.train=ind.train,ind.test=ind.test,alpha.value=alpha,STDx=STDx,penalty.factor)
		seliter <- c(seliter,res3fCv$selVar)
		MSEcv[j]<-res3fCv$MSEi
		CORRcv[j]<-res3fCv$CORR
		noFeaSelected[j]<-length(res3fCv$selVar)
		FeaSelected[, j] <- res3fCv$coefEstim
	}
	feaNames<-names(table(seliter))
	FeaSelFreq<-as.vector(table(seliter))
	return(list(FeaSelFreq=FeaSelFreq,MSEcv=MSEcv,CORRcv=CORRcv,ncv=ncv,feaNames=feaNames, FeaSelected=FeaSelected))
}