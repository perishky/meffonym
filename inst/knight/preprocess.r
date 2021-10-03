x <- read.csv("13059_2016_1068_MOESM3_ESM.csv", stringsAsFactors=F)

coefs <- x[,c("CpGmarker","CoefficientTraining")]
colnames(coefs) <- c("cpg","coef")

write.csv(coefs, file="coefs.csv", row.names=F)


