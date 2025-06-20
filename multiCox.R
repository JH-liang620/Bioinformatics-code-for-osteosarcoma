setwd("")
library(survival)
library(survminer)

rt=read.table("lassoSigExp.txt",header=T,sep="\t",check.names=F,row.names=1) 

multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)

pdf(file="forest.pdf", onefile=FALSE,
    width = 5.5,           
    height = 2.8,          
)
ggforest(multiCox,
         main = "Hazard ratio",
         cpositions = c(0.02, 0.1, 0.3), 
         fontsize = 1, 
         refLabel = "reference", 
         noDigits = 2)
dev.off()

riskScore=predict(multiCox,type="risk",newdata=rt) 
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
medianTrainRisk=median(riskScore)
risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))
write.table(cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore,risk)),
            file="riskTrain.txt",
            sep="\t",
            quote=F,
            row.names=F)

rtTest=read.table("expTimeGEO1.txt",header=T,sep="\t",check.names=F,row.names=1)   
rtTest[,"futime"]=rtTest[,"futime"]/365
riskScoreTest=predict(multiCox,type="risk",newdata=rtTest)    
riskTest=as.vector(ifelse(riskScoreTest>medianTrainRisk,"high","low"))
write.table(cbind(id=rownames(cbind(rtTest[,outCol],riskScoreTest,riskTest)),cbind(rtTest[,outCol],riskScore=riskScoreTest,risk=riskTest)),
            file="riskTestGEO1.txt",
            sep="\t",
            quote=F,
            row.names=F)

rtTest=read.table("expTimeGEO2.txt",header=T,sep="\t",check.names=F,row.names=1)    
rtTest[,"futime"]=rtTest[,"futime"]/365
riskScoreTest=predict(multiCox,type="risk",newdata=rtTest)    
riskTest=as.vector(ifelse(riskScoreTest>medianTrainRisk,"high","low"))
write.table(cbind(id=rownames(cbind(rtTest[,outCol],riskScoreTest,riskTest)),cbind(rtTest[,outCol],riskScore=riskScoreTest,risk=riskTest)),
            file="riskTestGEO2.txt",
            sep="\t",
            quote=F,
            row.names=F)