# Project: Mimomics
# this script has funtions from Pushpike's package, without using qSTAR data utils
# which means it just fits elastic net/LASSO and performs cross-validation
# and makes the final prediction with the features' selection 
# Author: lucp2149
###############################################################################
setwd("C://Mimomics//FinRisk")
library(glmnet)
#fpnames is replaced with varnames; output selGene is replaced with selVar
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
	#HERE, ALL FEATURES WERE COLLAPSED; 
	#INSTEAD CONCATENATE THEM A VECTOR TO SEE WHEN THEY WERE DISCOVERED JOINTLY
	feaNames<-names(table(seliter))
	FeaSelFreq<-as.vector(table(seliter))
#return(no.fp.cv)
	return(list(FeaSelFreq=FeaSelFreq,MSEcv=MSEcv,CORRcv=CORRcv,ncv=ncv,feaNames=feaNames, FeaSelected=FeaSelected))
}
#innerouterC, without output of an object, use list instead
InnerOuterCV<-function (y, x, alpha, STDx, ncv = 10, nFold = 3, innerCV=100, AtLeast=0.05, penalty.factor) 
{
	#initialize some vectors
	CORRcv <- MSEcv <- rep(NA, ncv)
	seliter <- numeric()
	seliterOutBag<-numeric()
	noFeaSelected <- rep(NA, ncv)
	MCorrCvInner<-rep(NA, ncv)
	MMseCvInner<-rep(NA, ncv)
	
	#deving the data into outbag and inbag
	nTest <- floor(length(y)/nFold)
	nTrain <- length(y) - nTest
	for (j in 1:ncv) {
		indbag <- as.vector(sort(sample(1:length(y), nTrain, 
								replace = F)))
		outbag <- as.vector(c(1:length(y))[-indbag])
		res3fCv <- NA
		res3fCv <- CvELnetOrLasso(y[indbag],x[indbag,],alpha,STDx=STDx,ncv=innerCV,nFold=nFold,penalty.factor)
		
		FreqFea <- res3fCv$FeaSelFreq
		names(FreqFea) <- res3fCv$feaNames
		Temp <- sort(FreqFea, decreasing = TRUE)
		if ((max(Temp)) < AtLeast * res3fCv$ncv) 
			stop("AtLeast: minimum  selection frequency is higher than maximum selection frequency")
		TopListLassoCV <- Temp[Temp > AtLeast * res3fCv$ncv]
		
		#mostly selected gene during the inner cross validations
		seliter <- c(seliter, names(TopListLassoCV))
		
		#median correlations of the inner cross valdiations
		MCorrCvInner[j]<-median(res3fCv$CORRcv,na.rm=T)
		MMseCvInner[j]<-median(res3fCv$MSEcv,na.rm=T)
		
		#Evaluate the model based on outofbag data
		OutbagEval<-FitElnetOrLassokFoldCV(X = x[, is.element(colnames(x), 
								names(TopListLassoCV))], y = y, varnames = names(TopListLassoCV), 
				ind.train = indbag, ind.test = outbag, alpha.value = alpha, 
				STDx = STDx, penalty.factor)
		
		MSEcv[j] <- OutbagEval$MSEi
		CORRcv[j] <- OutbagEval$CORR
		seliterOutBag<-c(seliterOutBag,OutbagEval$selVar)
		
	}
	feaNamesInner <- names(table(seliter))
	InnerFeaSelFreq <- as.vector(table(seliter))
	
	feaNamesOutbag <- names(table(seliterOutBag))
	OuterFeaSelFreq <- as.vector(table(seliterOutBag))
	
	return(list(InnerFeaSelFreq = InnerFeaSelFreq,OuterFeaSelFreq=OuterFeaSelFreq, MSEcv = MSEcv, 
					CORRcv = CORRcv,MCorrCvInner=MCorrCvInner,MMseCvInner=MMseCvInner, ncv = ncv, feaNamesInner = feaNamesInner,feaNamesOutbag=feaNamesOutbag))
}

#final cross validation and feature selection

FinalEvaluation<-function(x,y,CvObject,AtLeast=0.10,MixPara=1,MinCorr=0.8,STDx=TRUE,penalty.factor){
	
	#  if (class(CvObject)!="InOutCV") stop("Invalid Class object")
	FreqFea<-CvObject$OuterFeaSelFreq
	names(FreqFea)<-CvObject$feaNamesOutbag
	Temp<-sort(FreqFea,decreasing = TRUE)
	#Temp[Temp>50]
	
	
	if ((max(Temp))<AtLeast*CvObject$ncv) stop("reduce the percentage of selection frequency")
	TopListLassoCV<-Temp[Temp>AtLeast*CvObject$ncv]
	
	EstCoef<-corrtop<-nSelected<-rep(NA,length(TopListLassoCV)-1)
	for (jj in 2:length(corrtop)){
		
		LassoOnCvSelected<-FitElnetOrLassokFoldCV(X=x[,is.element(colnames(x),names(TopListLassoCV)[1:jj])],
				y=y,varnames=names(TopListLassoCV),ind.train=1:nrow(x) ,ind.test=1:nrow(x) ,alpha.value=MixPara,STDx=STDx, 
				penalty.factor=rep(1, sum(is.element(colnames(x), names(TopListLassoCV)[1:jj]))))
		nSelected[jj]<-jj
		corrtop[jj]<-LassoOnCvSelected$CORR
	}
	
	if (max(corrtop, na.rm = T) < MinCorr) {
		MinCorr<-quantile(corrtop, na.rm = T,probs=0.75)
		warning("Value specified in argument MinCorr is higher than the maximum correlation obtained with top K features.")
	}
	corrtop<-sort(corrtop)
	HowMany<-min(nSelected[corrtop>MinCorr],na.rm=T)
	
	
	LassoOnTopKSelected<-FitElnetOrLassokFoldCV(X=x[,is.element(colnames(x),names(TopListLassoCV)[1:HowMany])],
			y=y,varnames=names(TopListLassoCV)[1:HowMany],ind.train=1:nrow(x) ,ind.test=1:nrow(x) ,alpha.value=MixPara,STDx=STDx,
			penalty.factor=rep(1, sum(is.element(colnames(x),names(TopListLassoCV)[1:HowMany]))))
	
	finalPred<-LassoOnTopKSelected$finalPred
	selGene<-LassoOnTopKSelected$selVar
	TopList<-TopListLassoCV
	
	#leave one out for final selected once.
	LoocvPredY<-rep(NA,length(y))
	for (j in 1:length(y) ) {
		fit1 <- glmnet(x[-j,is.element(colnames(x),names(TopListLassoCV)[1:HowMany])],y[-j],alpha = MixPara,standardize = STDx)
		pred.fit1 <- predict(fit1,matrix(x[j,is.element(colnames(x),names(TopListLassoCV)[1:HowMany])],ncol=HowMany),s=(1:100)*0.01)
		PRMSS<-colMeans((pred.fit1 - y[j])^2)
		OptIndex<-which.min(PRMSS)/100
		LoocvPredY[j]<- predict(fit1,matrix(x[j,is.element(colnames(x),names(TopListLassoCV)[1:HowMany])],ncol=HowMany),s=OptIndex)
	}
	
	return(list(finalPred=finalPred,selGene=selGene,nSelected=nSelected,corrtop=corrtop,HowMany=HowMany,TopList=TopList,LoocvPredY=LoocvPredY,y=y))
}

#plot functions
#setMethod("plot", signature(x="FinalRes", y="missing"),
plotFinalRes <-	function(x,  y, ...) {
			#if (class(x)!="FinalRes") stop("Invalid class object")
			nSelected<-x$nSelected
			corrtop<-x$corrtop
			dotsCall <- substitute(list(...))
			ll <- eval(dotsCall)
			if(!hasArg("ylab")) ll$ylab <- "Correlation"
			if(!hasArg("xlab")) ll$xlab <- "Top K"
			ll$type="l"
			ll$x<-nSelected[-1]
			ll$y<-corrtop
			par(mfrow=c(1,2),mar=c(4,4,5,2))
			do.call(plot,args=ll)
			points(nSelected[-1],corrtop,pch=19,col="blue")
			
			matplot(cbind(x$y,x$finalPred,x$LoocvPredY),type="l",col=1:3,ylab="Expression value",xaxt="n",main=paste("Observed Vs Predicted\n Selected Top ",
							length(x$selGene)," Variables",sep=""))
			axis(1,at=1:length(x$y),labels=names(x$y),las=2,cex.axis=0.5,cex.names=0.65)
			legend("topright",c("Observed","predicted","Pred LOOCV"),col=1:3,lty=1:3)
			
		}
#)

#----------------------------------------------------------------------------------
#setMethod("plot", signature(x="CVLassoElnet", y="missing"),
plotCVLassoElnet <-	function(x,  y, ptype=1, AtLeast=0.10, ...) {
			#if (class(x)!="CVLassoElnet") stop("Invalid class object")
			
			if (ptype==1) {
				dotsCall <- substitute(list(...))
				ll <- eval(dotsCall)
				ll$main<-"Density of Correlation"
				ll$xlab <- expression(rho)
				if(!hasArg("ylab")) ll$ylab <- "density"
				if(!hasArg("cex.lab")) ll$cex.lab <- 0.8
				if(!hasArg("cex.main")) ll$cex.main <- 1
				if(!hasArg("col")) ll$col <- "blue"
				dcor<-density(x$CORRcv,from=-1,to=1,na.rm=T)
				ll$x<-dcor
				par(mfrow=c(1,2))
				
				do.call(plot,args=ll)
#plot(dcor,main="Density of Correlation",xlab=expression(rho),col="blue")
#polygon(dcor, col="gray", border="blue") 
				
				ll$main<-"Density of MSE"
				ll$xlab <- "MSE"
				dcor<-density(x$MSEcv,na.rm=T)
				ll$x<-dcor
				do.call(plot,args=ll)
				# polygon(dcor, col="gray", border="blue")
			}
			
			if (ptype==2) {
				par(mfrow=c(1,1),mar=c(8,4,4,2))
				MsGe<-x$FeaSelFreq
				names(MsGe)<-x$feaNames 
				maxG<-max(MsGe)
				cutoff<-x$ncv*AtLeast
				if ((x$ncv*AtLeast)>=maxG)  {
					warning("percentage of maximum selection frequency ",round(max(MsGe)/x$ncv,3)*100,"%" )
					cutoff<-maxG-2
				}
				HFreqGenes<-MsGe[MsGe>cutoff]
				
				
				dotsCall <- substitute(list(...))
				lll <- eval(dotsCall)
				sHFreqGenes<-sort(HFreqGenes,decreasing=T)
				lll$height<-sHFreqGenes
				if(!hasArg("xlab" ))   lll$xlab<-""
				if(!hasArg("ylab" ))   lll$ylab<-"Frequency"
				if(!hasArg("main"))    lll$main<-"Mostly Selected Genes"
				lll$col<-rainbow(maxG)
				lll$names<-names(sHFreqGenes)
				if(!hasArg("cex.lab")) lll$cex.lab <- 1
				if(!hasArg("las"))   lll$las<-2
				lll$cex.names<-0.65
				do.call(barplot,args=lll) 
			}
			
		}#)



#----------------------------------------------------------------------------------
#setMethod("plot", signature(x="InOutCV", y="missing"),
plotInOutCV <-	function(x,  y, ptype=1, ...) {
			#if (class(x)!="InOutCV") stop("Invalid class object")
					
			if (ptype==1) {
				
				Freq<-x$OuterFeaSelFreq
				names(Freq)<-x$feaNamesOutbag
				sFreq<-sort(Freq,decreasing=T)[1:71]
				
				par(mfrow=c(1,1),mar=c(8,4,4,2))
				dotsCall <- substitute(list(...))
				ll <- eval(dotsCall)
				if(!hasArg("main")) ll$main<-"Feature selection frequency on out of bag data"
				if(!hasArg("ylab")) ll$ylab <- "frequency"
				if(!hasArg("cex.main")) ll$cex.main <- 1
				if(!hasArg("col")) ll$col <- rainbow(length(sFreq))
				if(!hasArg("cex.lab")) ll$cex.lab <- 1
				if(!hasArg("las"))   ll$las<-2
				if(!hasArg("cex.names")) ll$cex.names<-0.65
				ll$height<-sFreq
				ll$names<-names(sFreq)
				do.call(barplot,args=ll)   
				
			}
			
			if (ptype==2) {
				
				par(mfrow=c(1,2))
				dotsCall <- substitute(list(...))
				lll <- eval(dotsCall)
				
				
				lll$x<-x$CORRcv
				lll$main<-"Correlations (Out Of Bag data)"
				lll$col=3
				
				do.call(boxplot,args=lll) 
				stripchart(x$MCorrCvInner, vertical = TRUE, method = "jitter",
						pch = 21, col = "maroon", bg = "bisque",
						add = TRUE) 
				
				lll$x<-x$MSEcv
				lll$main<-"MSE (Out Of Bag data)"
				lll$col=3
				do.call(boxplot,args=lll) 
				stripchart(x$MMseCvInner, vertical = TRUE, method = "jitter",
						pch = 21, col = "maroon", bg = "bisque",
						add = TRUE) 
				
				
				
			}
			
		}#)
