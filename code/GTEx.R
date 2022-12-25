#section 1 extract GTEx data
library(data.table)
library(DESeq2)
pd <- fread("GTEX_phenotype.gz")
pd <- as.data.frame.matrix(pd)
table(pd[,3])
pd_N <- pd[pd[,3] == "Colon",]
pd_n <- as.character(pd_N$Sample)
exp <- fread("gtex_RSEM_gene_fpkm.gz")
exp <- as.data.frame.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:ncol(exp)]
exp_n <- exp[,colnames(exp) %in% pd_n]
exp_nc <- (2^exp_n - 0.001)
exp_nc[exp_nc<0] <- 0
write.table(exp_nc,file="GeneFPKMExp.txt",sep="\t",quote=F)

#section 2 delete genes with low expression
inputFile="GeneFPKMExpSymbol.txt"
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
write.table(data,file="GTEx.txt",sep="\t",quote=F)
