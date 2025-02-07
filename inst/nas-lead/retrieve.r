suptab = "https://pmc.ncbi.nlm.nih.gov/articles/instance/7170756/bin/NIHMS1540834-supplement-SupTab3-4.xlsx"
if (!file.exists(basename(suptab)))
    download.file(suptab,basename(suptab))

library(readxl)

patella = as.data.frame(read_xlsx(basename(suptab),sheet="Table S3 Patella_EN_coeff"))
tibia = as.data.frame(read_xlsx(basename(suptab),sheet="Table S4 Tibia_EN_coeff"))

cols = c("CpG Name", "Elastic Net Coeff")
patella = patella[,cols]
tibia = tibia[,cols]

colnames(patella) = colnames(tibia) = c("pred.var","coef")

patella = patella[abs(patella$coef) > 0,]

write.csv(patella,row.names=F,file="patella.csv")
write.csv(tibia,row.names=F,file="tibia.csv")
