library(survival)
library(survminer)
library(timeROC)
setwd("")   

bioROC=function(inputFile=null, rocFile=null){

	rt=read.table(inputFile, header=T, sep="\t", check.names=F)

	ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,
	               marker=rt$riskScore,cause=1,
	               weighting='aalen',
	               times=c(1,3,5),ROC=TRUE)

	pdf(file=rocFile, width=3.2, height=3.2)

	par(oma=c(0,0,0,0)) 

	par(mar=c(4,4,0,0) + 0.1)
	plot(ROC_rt,time=1,col="#E64B35",title=FALSE,lwd=2)
	plot(ROC_rt,time=3,col="#4DBBD5",add=TRUE,title=FALSE,lwd=2)
	plot(ROC_rt,time=5,col="#3C5488",add=TRUE,title=FALSE,lwd=2)

	legend('bottomright',
	        c(paste0('1-year AUC value: ',sprintf("%.03f",ROC_rt$AUC[1])),
	          paste0('3-year AUC value: ',sprintf("%.03f",ROC_rt$AUC[2])),
	          paste0('5-year AUC value: ',sprintf("%.03f",ROC_rt$AUC[3]))),
	        col=c("#E64B35","#4DBBD5","#3C5488"),lwd=2,bty = 'n')
	dev.off()
}

bioROC(inputFile="riskTrain.txt", rocFile="ROC.train.pdf")
bioROC(inputFile="riskTestGEO1.txt", rocFile="ROC.test1.pdf")
bioROC(inputFile="riskTestGEO2.txt", rocFile="ROC.test2.pdf")