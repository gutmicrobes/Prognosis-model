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
