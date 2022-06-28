## Levine https://dx.doi.org/10.18632/aging.101414
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5940111/bin/aging-10-101414-s002.csv
levine <- read.csv("aging-10-101414-s002.csv",stringsAsFactors=F)
levine <- levine[,c("CpG","Weight")]
colnames(levine) <- c("pred.var","coef")

write.csv(levine, file="coefs.csv", row.names=F)
