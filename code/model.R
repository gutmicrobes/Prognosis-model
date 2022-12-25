#section 1 LASSO
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

#section 2 Multivariate cox analysis
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

#section 3 distribution and survival analysis of training dataset
library(survival)
library("survminer")
library(tidyverse)
rt=read.table("GeneMultiCox.txt",header=T,sep="\t")
res.cut <- surv_cutpoint(rt,time="futime",
                         event="fustat",variables=c("riskScore"),minprop = 0.2)
summary(res.cut)
plot(res.cut,"riskScore",palette="lancet")
cutoff=as.numeric(res.cut$cutpoint[1])
rt$risk=ifelse(rt[,"riskScore"]<=cutoff, "low", "high")
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
pdf(file="survival1.pdf",onefile = FALSE, 
       width = 5.5,
       height =5)
ggsurvplot(fit, 
           data=rt,
           conf.int=FALSE,
           pval=paste0("p=",pValue),
           pval.size=6,
           risk.table=TRUE,
           legend.labs=c("High risk", "Low risk"),
           legend.title="Risk",
           xlab="Time(years)",
           break.time.by = 1,
           risk.table.title="",
           palette=c("red", "blue"),
           risk.table.height=.25)
dev.off()
write.table(rt,
            file="GeneMultiCoxrisk.txt",
            sep="\t",
            quote=F,
            row.names=F)
summary(fit)

#section 4 plot heatmap and survstat map
library(pheatmap)
rt=read.table("GTExTCGA9GeneandGroup.txt",sep="\t",header=T,row.names=1,check.names=F)
rt1=rt[c(1:(ncol(rt)-1))]
rt1=log2(t(rt1)+0.001)
annotation=data.frame(sample=rt[,ncol(rt)])
rownames(annotation)=rownames(rt)
pdf(file="heatmap1.pdf",width = 10,height = 4)
pheatmap(rt1, 
         annotation=annotation, 
         cluster_cols = FALSE,
         fontsize_row=11,
         show_colnames = F,
         fontsize_col=3,
         color = colorRampPalette(c("blue", "white", "red"))(50) )
dev.off()
library(pheatmap)
rt=read.table("GeneMultiCoxrisk.txt",sep="\t",header=T,row.names=1,check.names=F)
rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
color=as.vector(rt$fustat)
color[color==1]="red"
color[color==0]="green"
pdf(file="survStat.pdf",width = 10,height = 4)
plot(rt$futime,
     pch=19,
     xlab="Patients ranked by risk socre",
     ylab="Patients' survival time (years)",
     col=color)
legend("topleft", c("Dead patients", "Alive patients"),bty="n",pch=19,col=c("red","green"),cex=1.2)
abline(v=lowLength,lty=2)
dev.off()

#section 5 ROC curve of training dataset
rm(list=ls())
library(timeROC)
library(survival)
rt=read.table("GeneMultiCoxrisk.txt",header=T,sep="\t",check.names=F,row.names=1)
result <- with(rt, timeROC(T = rt$futime,
                           delta = rt$fustat,
                           marker = rt$riskScore,
                           cause = 1,
                           weighting = "marginal",
                           times = c(1,3,5),
                           iid = TRUE))
dat = data.frame(fpr = as.numeric(result$FP),tpr = as.numeric(result$TP),
                 time = rep(as.factor(c(1,3,5)),each = nrow(result$TP)))
library(ggplot2)
gg<-ggplot()+
  geom_line(data = dat,aes(x=fpr,y=tpr,color=time),size=1)+
  scale_color_manual(name=NULL,values=c("red","green","blue"),
                     labels=paste0("AUC of ",c(1,3,5),"-year survival: ",
                                   format(round(result$AUC,2),nsmall=2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color="grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.background = element_rect(linetype=1,size=0.2,colour="black"),
        legend.position=c(0.765,0.125))+
  scale_x_continuous(expand=c(0.005,0.005))+
  scale_y_continuous(expand=c(0.005,0.005))+
  labs(x="False positive rate",
       y="True positive rate")+
  coord_fixed()
library(ggThemeAssist)
ggg<-gg + theme(axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.text = element_text(size = 20))
ggg + theme(legend.position = c(0.7, 0.125))

#section 6 Multivariate independent prognostic analysis of training dataset
library(survival)
rt=read.table("indepInput.txt",header=T,sep="\t",check.names=F,row.names=1)
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCoxSum=summary(multiCox)
outTab=data.frame()
outTab=cbind(
             HR=multiCoxSum$conf.int[,"exp(coef)"],
             HR.95L=multiCoxSum$conf.int[,"lower .95"],
             HR.95H=multiCoxSum$conf.int[,"upper .95"],
             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
pValue=outTab[,"pvalue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)
