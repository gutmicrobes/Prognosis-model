library(survival)
library("survminer")
library(tidyverse)

rt=read.table("TCGAtestCox.txt",header=T,sep="\t")

res.cut <- surv_cutpoint(rt,time="futime",
                         event="fustat",variables=c("riskScore"))
summary(res.cut)
plot(res.cut,"riskScore",palette="npg")

cutoff=as.numeric(res.cut$cutpoint[1])
rt$risk=ifelse(rt[,"riskScore"]<=cutoff, "low", "high")


diff=survdiff(Surv(futime, fustat) ~risk,data = rt)

pValue=1-pchisq(diff$chisq,df=1)

pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)


fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)


pdf(file="survival.pdf",onefile = FALSE, 
       width = 5.5,
       height =5)
ggsurvplot(fit, 
           data=rt,
           conf.int=FALSE,
           pval=paste0("p=",pValue),
           pval.size=4,
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
