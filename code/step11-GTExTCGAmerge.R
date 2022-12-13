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