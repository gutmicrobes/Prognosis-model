library("glmnet")
library("survival")
rt=read.table("TCGAtrainexpTime.txt",header=T,sep="\t",row.names=1)

x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))
fit <- glmnet(x, y, family = "cox", maxit = 1000)
pdf("lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
pdf("cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
write.table(lassoGene,file="lassoGene.txt",sep="\t",quote=F,row.names=F,col.names=F)

opar <- par(no.readonly = TRUE)
par(mfrow=c(1,2))
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
plot(fit, xvar = "lambda", label = TRUE)
par(opar)