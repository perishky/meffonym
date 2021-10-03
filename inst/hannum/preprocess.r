x <- read.csv("hannum-model.csv",stringsAsFactors=F)
coefs <- x[,c("Marker","Coefficient")]
colnames(coefs) <- c("cpg","coef")
write.csv(coefs, file="coefs.csv",row.names=F)
