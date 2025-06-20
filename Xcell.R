library(reshape2)
library(ggpubr)
library(limma)
expFile="mRNA_tpm.txt"
riskFile="risk.txt" 
xCellFile="xCellsig.txt"
setwd("") 

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

xCell=read.table(xCellFile, header=T, sep="\t", check.names=F, row.names=1)

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

xCell=t(xCell)
sameSample=intersect(row.names(xCell), row.names(risk))
xCell=xCell[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
scorerisk=cbind(xCell, risk)

data=melt(scorerisk, id.vars=c("risk"))
colnames(data)=c("risk", "Immune", "Fraction")

bioCol=c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000")
bioCol=bioCol[1:length(levels(factor(data[,"risk"])))]
p=ggboxplot(data, x="Immune", y="Fraction", color="risk", 
            ylab="xCell score",
            xlab="",
            legend.title="Risk",
            palette=bioCol)
p=p+rotate_x_text(50)

pdf(file="boxplot.pdf", width=5.5, height=4)
p+stat_compare_means(aes(group=risk),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")
dev.off()
