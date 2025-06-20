library(limma)
library(ggplot2)

expFile="uniSigGeneExp.txt"   
clusterFile="cencluster.txt"  
setwd("")  

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=t(data)

data.pca=prcomp(data, scale. = TRUE)
pcaPredict=predict(data.pca)
write.table(pcaPredict, file="newTab.xls", quote=F, sep="\t")

cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
Cluster=as.vector(cluster[,1])

bioCol=c("#E64B35FF","#4DBBD5FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF")
prgCluCol=bioCol[1:length(levels(factor(Cluster)))]

PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Cluster=Cluster)
PCA.mean=aggregate(PCA[,1:2], list(Cluster=PCA$Cluster), mean)
pdf(file="PCA.pdf", width=3, height=3)
ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Cluster)) +
	scale_colour_manual(name="Cluster", values =prgCluCol)+
    theme_bw()+
    theme(legend.position=c(0.1,0.1), plot.margin=unit(rep(0.5,4),'lines'))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
