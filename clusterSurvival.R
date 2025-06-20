library(survival)
library(survminer)

cenclusterFile="cencluster.txt" 
cliFile="time.txt"  
setwd("")  

cencluster=read.table(cenclusterFile, header=T, sep="\t", check.names=F, row.names=1)
rownames(cencluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(cencluster))
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365

sameSample=intersect(row.names(cencluster), row.names(cli))
rt=cbind(cli[sameSample,,drop=F], cencluster[sameSample,,drop=F])

length=length(levels(factor(rt$cencluster)))
diff=survdiff(Surv(futime, fustat) ~ cencluster, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ cencluster, data = rt)

bioCol=c("#E64B35FF","#4DBBD5FF","#3C5488FF")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=F,
		           pval=pValue,
		           pval.size=5.5,
		           legend.title="Cluster",
		           legend.labs=levels(factor(rt[,"cencluster"])),
		           legend = c(0.85, 0.85),
		           font.legend=14,
		           xlab="Time(years)",
		           break.time.by = 1,
		           palette = bioCol,
		           risk.table=F,
		           cumevents=F,
		           risk.table.height=.25)
pdf(file="survival.pdf",onefile = FALSE,width=3,height=3)
print(surPlot)
dev.off()
