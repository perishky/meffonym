## https://doi.org/10.18632/aging.101508
## Supplementary Dataset 2
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6075434/bin/aging-10-101508-s005.csv

coefs <- read.csv("aging-10-101508-s005.csv")
coefs <- coefs[,c("ID", "Coef")]
colnames(coefs) <- c("pred.var","coef")

write.csv(coefs, file="coefs.csv", row.names=F)

