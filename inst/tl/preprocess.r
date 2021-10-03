library(readxl)
## Lu's DNAmTL (140 CpG sites) https://dx.doi.org/10.18632/aging.102173
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6738410/bin/aging-11-102173-s003.xlsx
coefs <- read_xlsx("aging-11-102173-s003.xlsx",skip=5)
coefs <- coefs[,c("Variable","Coefficient")]
colnames(coefs) <- c("cpg","coef")

write.csv(coefs, file="coefs.csv", row.names=F)
