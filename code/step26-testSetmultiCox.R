library(survival)
rt=read.table("multiInput.txt",header=T,sep="\t",check.names=F,row.names=1)
rt$futime=rt$futime/365
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
outTab=gsub("`","",outTab)
write.table(outTab,file="immuneGeneMultiCox.xls",sep="\t",row.names=F,quote=F)

riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)

write.table(cbind(id=rownames(cbind(rt[,outCol],riskScore)),cbind(rt[,outCol],riskScore)),
    file="GeneMultiCox.txt",
    sep="\t",
    quote=F,
    row.names=F)

newrt=read.table("TCGAtest.txt",header=T,sep="\t",check.names=F,row.names=1)
newrt$futime=newrt$futime/365
riskScore=predict(multiCox,type="risk",newdata=newrt)


write.table(cbind(id=rownames(cbind(newrt[,outCol],riskScore)),cbind(newrt[,outCol],riskScore)),
            file="TCGAtestCox.txt",
            sep="\t",
            quote=F,
            row.names=F)
