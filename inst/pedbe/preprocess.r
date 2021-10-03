## https://doi.org/10.1073/pnas.1820843116
## https://raw.githubusercontent.com/kobor-lab/Public-Scripts/master/datcoefInteresting94.csv

coefs <- read.csv("datcoefInteresting94.csv")
coefs <- coefs[,c("ID","Coef")]
colnames(coefs) <- c("cpg","coef")
write.csv(coefs, file="coefs.csv", row.names=F)


