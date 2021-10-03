## MiAge
## 10.1080/15592294.2017.1389361
## http://www.columbia.edu/~sw2206/softwares/mitotic_age_R_code.zip
sites <- read.csv("code/Additional_File1.csv", stringsAsFactors=F)$CpG_site_ID
load("code/site_specific_parameters.Rdata")

coefs <- data.frame(
    cpg=sites,
    b=methyl.age[[1]],
    c=methyl.age[[2]],
    d=methyl.age[[3]],
    stringsAsFactors=F)

write.csv(coefs, file="coefs.csv", row.names=F)
