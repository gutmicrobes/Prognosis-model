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

pdf(file="vol.pdf",height=5,width=5)
xMax=max(abs(as.numeric(as.vector(immuneDiffAll$logFC))))
yMax=max(-log10(immuneDiffAll$fdr))+1
plot(as.numeric(as.vector(immuneDiffAll$logFC)), -log10(immuneDiffAll$fdr), xlab="logFC",ylab="-log10(fdr)",
     main="Volcano", ylim=c(0,yMax),xlim=c(-xMax,xMax),yaxs="i",pch=20, cex=0.8)
diffSub=subset(immuneDiffAll, fdr<fdrFilter & as.numeric(as.vector(logFC))>logFCfilter)
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="#1500FF",cex=0.8)
diffSub=subset(immuneDiffAll, fdr<fdrFilter & as.numeric(as.vector(logFC))<(-logFCfilter))
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="#FF0102",cex=0.8)
abline(v=0,lty=2,lwd=3)
dev.off()

library(pheatmap)
geneNum=50
diffSig=immuneDiffResult[order(as.numeric(as.vector(immuneDiffResult$logFC))),]
diffGeneName=as.vector(rownames(diffSig))
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(geneNum*2)){
  hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
  hmGene=diffGeneName
}
hmExp=hmExp[hmGene,]

hmExp=log2(hmExp+0.001)
Type=c(rep("N",conNum),rep("T",treatNum))
names(Type)=colnames(diffExp)
Type=as.data.frame(Type)
pdf(file="heatmap1.pdf",height=12,width=15)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         show_rownames = F,
         fontsize = 12,
         fontsize_row=3,
         fontsize_col=10)
dev.off()