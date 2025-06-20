library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

expFile="mRNA_logfpkm.txt"           
riskFile="riskTrain.txt"    
gmtFile="c5.go.v7.4.symbols.gmt"    
setwd("")     

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

Risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
data=data[,row.names(Risk)]

dataL=data[,row.names(Risk[Risk[,"risk"]=="low",])]
dataH=data[,row.names(Risk[Risk[,"risk"]=="high",])]
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
logFC=sort(logFC,decreasing=T)
genes=names(logFC)

gmt=read.gmt(gmtFile)

kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
write.table(kkTab,file="GSEA.result.txt",sep="\t",quote=F,row.names = F)
	
termNum=5  
kkUp=kkTab[kkTab$NES>0,]
if(nrow(kkUp)>=termNum){
	showTerm=row.names(kkUp)[1:termNum]
	gseaplot=gseaplot2(kk, showTerm, base_size=12, title="Enriched in high risk group")
	pdf(file="GSEA.highRisk G.pdf", width=6, height=4)
	print(gseaplot)
	dev.off()
}

termNum=5   
kkDown=kkTab[kkTab$NES<0,]
if(nrow(kkDown)>=termNum){
	showTerm=row.names(kkDown)[1:termNum]
	gseaplot=gseaplot2(kk, showTerm, base_size=12, title="Enriched in low risk group")
	pdf(file="GSEA.lowRisk G.pdf", width=6, height=4)
	print(gseaplot)
	dev.off()
}