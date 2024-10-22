## Cheng diabetes cox lasso (145 CpG sites) 
#https://static-content.springer.com/esm/art%3A10.1038%2Fs43587-023-00391-4/MediaObjects/43587_2023_391_MOESM3_ESM.xlsx
cheng <- read.csv("43587_2023_391_MOESM3_ESM.csv", skip = 2) # supplementary table 5
cheng <- cheng[,c("Name","Beta")]
colnames(cheng) <- c("pred.var","coef")

write.csv(cheng, file="coefs.csv", row.names=F)
