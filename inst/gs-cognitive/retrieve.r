
library(data.table)

model = as.data.frame(fread("https://gitlab.com/danielmccartney/ewas_of_cognitive_function/-/raw/master/EWAS_Betas/g_processed_Mean_Beta_PIP.txt"))

model = model[model$PIP > 0.05,]
model = model[,c("CpG","Mean_Beta")]
colnames(model) = c("pred.var","coef")
write.csv(model, file="coefs.csv", row.names=F)