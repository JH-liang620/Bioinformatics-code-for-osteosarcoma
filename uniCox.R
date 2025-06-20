library(limma)
library(survival)

expFile="cenExp.txt"  
cliFile="time.txt"    
setwd("")

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)           
data=data[rowMeans(data)>0.5,] 
data1=t(data)
rownames(data1)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data1))

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)  
cli$futime=cli$futime/365

sameSample=intersect(row.names(data1), row.names(cli))
data1=data1[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(cli, data1)

outTab=data.frame()
sigGenes=c()
for(i in colnames(rt[,3:ncol(rt)])){
	cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
	coxSummary = summary(cox)
	coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	if(coxP<0.05){
		sigGenes=c(sigGenes,i)
		outTab=rbind(outTab,
				         cbind(id=i,
				         HR=coxSummary$conf.int[,"exp(coef)"],
				         HR.95L=coxSummary$conf.int[,"lower .95"],
				         HR.95H=coxSummary$conf.int[,"upper .95"],
				         pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
				        )
	}
}

write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)

sigGeneExp=data[sigGenes,]
sigGeneExp=rbind(id=colnames(sigGeneExp), sigGeneExp)
write.table(sigGeneExp, file="uniSigGeneExp.txt", sep="\t", quote=F, col.names=F)

sigExpTime=rt[,c("futime", "fustat", sigGenes)]
sigExpTime=rbind(id=colnames(sigExpTime), sigExpTime)
write.table(sigExpTime, file="uniSigExpTime.txt", sep="\t", quote=F, col.names=F)
