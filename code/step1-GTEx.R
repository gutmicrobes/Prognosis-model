Sys.setenv(LANGUAGE = "en")
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
