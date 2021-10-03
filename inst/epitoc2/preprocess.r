## epiTOC2
## 10.1186/s13073-020-00752-3
## download code and model from here 10.5281/zenodo.2632938

## code:
## download.file("https://zenodo.org/record/2632938/files/epiTOC2.R?download=1", destfile="epitoc2.r")

download.file("https://zenodo.org/record/2632938/files/dataETOC2.Rd?download=1", destfile="dataETOC2.Rd")
load("dataETOC2.Rd")
x <- as.data.frame(dataETOC2.l$epiTOC2)
x$cpg <- rownames(x)
write.csv(x, file="coefs.csv", row.names=F)






