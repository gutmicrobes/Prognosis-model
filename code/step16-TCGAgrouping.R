
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