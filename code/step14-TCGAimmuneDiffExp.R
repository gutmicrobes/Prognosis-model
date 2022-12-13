
diffExp=read.table("TCGA1.txt",sep="\t",header=T,check.names=F,row.names=1)

gene=read.table("immunediffgenelist.txt",sep="\t",header=F)

immuneDiffGene=intersect(gene[,1],row.names(diffExp))

hmExp=diffExp[immuneDiffGene,]

immuneGeneExp=rbind(ID=colnames(hmExp),hmExp)
write.table(immuneGeneExp,file="immuneGeneExp.txt",sep="\t",col.names=F,quote=F)
