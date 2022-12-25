#section 1 calculate risk score of testing dataset
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
newrt=read.table("TCGAtest.txt",header=T,sep="\t",check.names=F,row.names=1)
newrt$futime=newrt$futime/365
riskScore=predict(multiCox,type="risk",newdata=newrt)
write.table(cbind(id=rownames(cbind(newrt[,outCol],riskScore)),cbind(newrt[,outCol],riskScore)),
            file="TCGAtestCox.txt",
            sep="\t",
            quote=F,
            row.names=F)

#section 2 distribution and survival analysis of testing dataset
library(survival)
library("survminer")
library(tidyverse)
rt=read.table("TCGAtestCox.txt",header=T,sep="\t")
res.cut <- surv_cutpoint(rt,time="futime",
                         event="fustat",variables=c("riskScore"))
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
            file="TCGAtestCoxrisk.txt",
            sep="\t",
            quote=F,
            row.names=F)
summary(fit)

#section 3 ROC curve of testing dataset
library(timeROC)
library(survival)
rt=read.table("TCGAtestCoxrisk.txt",header=T,sep="\t",check.names=F,row.names=1)
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

#section 4 calculate risk score of entire TCGA dataset
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
newrt=read.table("TCGAexptime.txt",header=T,sep="\t",check.names=F,row.names=1)
newrt$futime=newrt$futime/365
riskScore=predict(multiCox,type="risk",newdata=newrt)
write.table(cbind(id=rownames(cbind(newrt[,outCol],riskScore)),cbind(newrt[,outCol],riskScore)),
            file="TCGAallCox.txt",
            sep="\t",
            quote=F,
            row.names=F)

#section 5 distribution and survival analysis of entire TCGA dataset
library(survival)
library("survminer")
library(tidyverse)
rt=read.table("TCGAallCox.txt",header=T,sep="\t")
res.cut <- surv_cutpoint(rt,time="futime",
                         event="fustat",variables=c("riskScore"))
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
            file="TCGAallCoxrisk.txt",
            sep="\t",
            quote=F,
            row.names=F)
summary(fit)

#section 6 ROC curve of entire TCGA dataset
library(timeROC)
library(survival)
rt=read.table("TCGAallCoxrisk.txt",header=T,sep="\t",check.names=F,row.names=1)
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
