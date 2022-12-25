#section 1 Delete genes with low expression
inputFile="TCGAmRNAsymbol.txt"
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
write.table(data,file="TCGA.txt",sep="\t",quote=F)

#section 2 delete TCGA normal sample
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

#section 3 Merge TCGA and GTEx
library(limma)
library(sva)
files=c("GTEx.txt", "TCGAgene.txt")
geneList=list()
for(i in 1:length(files)){
    inputFile=files[i]
    rt=read.table(inputFile, header=T, sep="\t",check.names=F)
    header=unlist(strsplit(inputFile, "\\.|\\-"))
    geneList[[header[1]]]=as.vector(rt[,1])
}
intersectGenes=Reduce(intersect, geneList)
allTab=data.frame()
for(i in 1:length(files)){
    inputFile=files[i]
    rt=read.table(inputFile, header=T, sep="\t", check.names=F)
    rt=as.matrix(rt)
    rownames(rt)=rt[,1]
    exp=rt[,2:ncol(rt)]
    dimnames=list(rownames(exp),colnames(exp))
    data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
    rt=avereps(data)
    if(i==1){
    	allTab=rt[intersectGenes,]
    }else{
    	allTab=cbind(allTab, rt[intersectGenes,])
    }
}
outTab=rbind(geneNames=colnames(allTab), allTab)
gtex=outTab[,1:308]
write.table(gtex, file="GTEx1.txt", sep="\t", quote=F, col.names=F)
tcga=outTab[,309:699]
write.table(tcga, file="TCGA1.txt", sep="\t", quote=F, col.names=F)
gtextcga=outTab[,1:699]
write.table(gtextcga, file="GTExTCGA1.txt", sep="\t", quote=F, col.names=F)

#section 4 differentially expressed genes
library("limma")
inputFile="GTExTCGA1.txt"
fdrFilter=0.05
logFCfilter=1
conNum=308
treatNum=391
outTab=data.frame()
grade=c(rep(1,conNum),rep(2,treatNum))
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.2,]
for(i in row.names(data)){
  geneName=unlist(strsplit(i,"\\|",))[1]
  geneName=gsub("\\/", "_", geneName)
  rt=rbind(expression=data[i,],grade=grade)
  rt=as.matrix(t(rt))
  wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  pvalue=wilcoxTest$p.value
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
	if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
		  outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
	 }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)
write.table(outTab,file="all.txt",sep="\t",row.names=F,quote=F)
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="diff.xls",sep="\t",row.names=F,quote=F)
heatmap=rbind(ID=colnames(data[as.vector(outDiff[,1]),]),data[as.vector(outDiff[,1]),])
write.table(heatmap,file="diffGeneExp.txt",sep="\t",col.names=F,quote=F)

#section 5 differentially expressed immune genes
fdrFilter=0.05
logFCfilter=1
conNum=308
treatNum=391
rt=read.table("all.txt",sep="\t",header=T,check.names=F,row.names=1)
diffExp=read.table("diffGeneExp.txt",sep="\t",header=T,check.names=F,row.names=1)
gene=read.table("gene.txt",sep="\t",header=F)
immuneDiffAll=rt[intersect(gene[,1],row.names(rt)),]
immuneDiffGene=intersect(gene[,1],row.names(diffExp))
hmExp=diffExp[immuneDiffGene,]
immuneDiffResult=immuneDiffAll[immuneDiffGene,]
immuneDiffResult=rbind(ID=colnames(immuneDiffResult),immuneDiffResult)
write.table(immuneDiffResult,file="immuneDiff.xls",sep="\t",col.names=F,quote=F)
immuneGeneExp=rbind(ID=colnames(hmExp),hmExp)
write.table(immuneGeneExp,file="immuneGeneExp.txt",sep="\t",col.names=F,quote=F)

#section 6 differentially immune genes expression
diffExp=read.table("TCGA1.txt",sep="\t",header=T,check.names=F,row.names=1)
gene=read.table("immunediffgenelist.txt",sep="\t",header=F)
immuneDiffGene=intersect(gene[,1],row.names(diffExp))
hmExp=diffExp[immuneDiffGene,]
immuneGeneExp=rbind(ID=colnames(hmExp),hmExp)
write.table(immuneGeneExp,file="immuneGeneExp.txt",sep="\t",col.names=F,quote=F)

#section 7 Divided into training group and verification group
data=read.table("expTime.txt",sep="\t",header=T,check.names=F,row.names=1)
set.seed(1)
train <- sample(nrow(data), nrow(data)*0.7)
test <- c(1:nrow(data))[-train]
train_data <- data[train,]
test_data <- data[test,]
traindata=rbind(ID=colnames(train_data),train_data)
write.table(traindata,file="TCGAtrainexpTime.txt",sep="\t",col.names=F,quote=F)
testdata=rbind(ID=colnames(test_data),test_data)
write.table(testdata,file="TCGAtestexpTime.txt",sep="\t",col.names=F,quote=F)