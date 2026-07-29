## Bhardwaj M, Sun YQ, Frick C, Schöttker B, Røe OD, Holleczek B, Mai
##   XM, Brenner H. Blood-based DNA methylation marker model for
##   short-term and long-term lung cancer risk prediction. BMC Med. 2026
##   Jun 6;24(1):344. doi: 10.1186/s12916-026-04973-y. PMID: 42251373;
##   PMCID: PMC13242670.
## Table 2

model = read.csv("pmc13242670-table2.csv")
model = model[,c("sites","coef")]
colnames(model) = c("pred.var","coef")

write.csv(model,file="bbdmm.csv",row.names=F)
