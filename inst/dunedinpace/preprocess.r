## remotes::install_github("danbelsky/DunedinPACE")
library("DunedinPACE")
model <- as.list(DunedinPACE::mPACE_Models)
pace.standard <- data.frame(
    name=model$gold_standard_probes$DunedinPACE,
    mean=model$gold_standard_means$DunedinPACE,
    stringsAsFactors=F)
write.csv(pace.standard, "standard.csv", row.names=F)
## note: the function DunedinPACE::PoAmProjector() refers
## to gold standard for the PoAm model, however this is an
## error because mPOA_Models does not exist in the package.
## The PoAm package DunedinPoAm38 does not refer to nor
## provide a gold standard either.

coefs <- model$model_weights$DunedinPACE
intercept <- model$model_intercept$DunedinPACE
model <- data.frame(
    pred.var=c("intercept",names(coefs)),
    coef=c(intercept, coefs))
write.csv(model, "coefs.csv", row.names=F)


