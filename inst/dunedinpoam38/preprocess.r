## devtools::install_github("danbelsky/DunedinPoAm38")
library("DunedinPoAm38")
model <- as.list(DunedinPoAm38::mPOA_Models)
coefs <- model$model_weights$DunedinPoAm_38
intercept <- model$model_intercept$DunedinPoAm_38
model <- data.frame(pred.var=c("intercept",names(coefs)),
                    coef=c(intercept, coefs))
write.csv(model, "coefs.csv", row.names=F)
