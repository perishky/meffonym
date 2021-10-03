## Zhang
## 10.1186/s13073-019-0667-1
## https://github.com/qzhang314/DNAm-based-age-predictor
## https://raw.githubusercontent.com/qzhang314/DNAm-based-age-predictor/master/en.coef
## https://raw.githubusercontent.com/qzhang314/DNAm-based-age-predictor/master/blup.coef

coef.en <- read.table("en.coef",stringsAsFactor=F,header=T)
colnames(coef.en) <- c("cpg","coef")
write.csv(coef.en, file="coef-en.csv", row.names=F)

coef.blup <- read.table("blup.coef",stringsAsFactor=F,header=T)
colnames(coef.blup) <- c("cpg","coef")
write.csv(coef.blup, file="coef-blup.csv", row.names=F)

## standardize methylation
## meth <- apply(meth, 1, scale)


