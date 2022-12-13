rm(list=ls())
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
ggplot()+
  geom_line(data = dat,aes(x=fpr,y=tpr,color=time),size=1)+
  scale_color_manual(name=NULL,values=c("#92C5DE","#F4A582","#66C2A5"),
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