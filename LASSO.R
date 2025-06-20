setwd("")
library("glmnet")
library("survival")

rt=read.table("uniSigExpTime.txt",header=T,sep="\t",row.names=1,check.names=F) 
rt$futime[rt$futime<=0]=1

x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))

set.seed(999)
fit <- glmnet(x, y, family = "cox", maxit = 100000)
pdf("lambda.pdf", width = 4, height = 4.5,)
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

cvfit <- cv.glmnet(x, y, family="cox", maxit = 100000, nfolds=10)
pdf("cvfit.pdf", width = 4, height = 4.5,)
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
lassoGene=c("futime","fustat",lassoGene)
lassoSigExp=rt[,lassoGene]
lassoSigExp=cbind(id=row.names(lassoSigExp),lassoSigExp)
write.table(lassoSigExp,file="lassoSigExp.txt",sep="\t",row.names=F,quote=F)