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
