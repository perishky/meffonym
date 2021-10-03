## Lee et al.
## https://doi.org/10.18632/aging.102049
## Supplementary File 1
## https://www.aging-us.com/article/102049/supplementary/SD2/0/aging-v11i12-102049-supplementary-material-SD2.csv

x <- read.csv("aging-v11i12-102049-supplementary-material-SD2.csv")

models <- list(cpc=x[,c("CpGs","Coefficient_CPC")],
               rpc=x[,c("CpGs","Coefficient_RPC")],
               rrpc=x[,c("CpGs","Coefficient_refined_RPC")])

for (id in names(models)) {
    colnames(models[[id]]) <- c("cpg","coef")
    models[[id]] <- models[[id]][which(abs(models[[id]]$coef) > 0),]
}

sapply(models, nrow)
## cpc  rpc rrpc 
## 547  559  396 

for (id in names(models)) {
    write.csv(models[[id]], file=paste0("coefs-", id, ".csv"), row.names=F)
