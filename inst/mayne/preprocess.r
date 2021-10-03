## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6040051
## 10.2217/epi-2016-0103
## Supplementary Table 5
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6040051/bin/epi-09-279-s6.xlsx
library(readxl)
x <- read_xlsx("epi-09-279-s6.xlsx") 
x <- x[-1,]
x <- x[-1,]
coefs <- data.frame(cpg=x[[1]], coef=as.numeric(as.character(x[[2]])))
write.csv(coefs, file="coefs.csv", row.names=F)
