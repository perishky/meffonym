
## Yang clock (epiTOC; 385 CpG sites) 
## https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1064-3/MediaObjects/13059_2016_1064_MOESM2_ESM.xls
library(readxl)
yang <- read_xls("13059_2016_1064_MOESM2_ESM.xls")
yang <- yang[[1]][-1]
write.csv(data.frame(cpg=yang, coef=1), file="coefs.csv", row.names=F)
