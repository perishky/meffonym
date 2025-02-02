
library(data.table)
model = fread("https://github.com/JacobBergstedt/CMVest/raw/refs/heads/main/model_coefficients.tsv")
model = as.data.frame(model)
colnames(model) = c("pred.var", "coef")

write.csv(
    model,
    file="coefs.csv",
    row.names=F)

