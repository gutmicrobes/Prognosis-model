inputFile="TCGA.txt"

library(limma)
library(pheatmap)

rt=read.table(inputFile,sep="\t",header=T,check.names=F)

rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

data=data[rowMeans(data)>0.2,]

tcgatumor=data[,42:514]
write.table(tcgatumor, file="TCGAtumor.txt", sep="\t", quote=F, col.names=T)
