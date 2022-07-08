
x <- read.csv("AdditionalFile3.csv", stringsAsFactors=F)
coefs <- x[,c("CpGmarker","CoefficientTraining")]
colnames(coefs) <- c("pred.var","coef")

write.csv(coefs, file="coefs.csv", row.names=F)
