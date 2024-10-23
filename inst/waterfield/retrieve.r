
#Scott Waterfield, Paul Yousefi, Matt Suderman
#DNA methylation models of protein abundance across the lifecourse
#medRxiv 2024.06.13.24308877; doi: https://doi.org/10.1101/2024.06.13.24308877

#https://www.medrxiv.org/content/medrxiv/early/2024/06/13/2024.06.13.24308877/DC2/embed/media-2.zip?download=true

sets <- list(
    "age24-450k"="Age24ModelCoef450K.csv",
    age24="Age24ModelCoef.csv",
    age9="Age9ModelCoef.csv",
    age48="MiddleAgeModelCoef.csv")

dir.create("models")

models <- do.call(rbind, lapply(sets, function(filename) {
    coefs <- read.csv(filename)
    coefs$Variable <- sub("\\(intercept\\)", "intercept", coefs$Variable,ignore.case=T)
    coefs[,c("Variable","Coefficient","Protein")]
}))

for (protein in unique(models$Protein)) {
    model <- models[models$Protein==protein,c("Variable","Coefficient")]
    colnames(model) <- c("pred.var","coef")
    write.csv(
        model,
        file=file.path("models",paste(protein, "coefs.csv", sep="-")),
        row.names=F)
}

proteins <- unique(models$Protein)

write.csv(
    data.frame(
        name=paste("waterfield",sub("_","-",proteins),sep="-"),
        tissue="blood",
        target=sub("_.*", "", proteins),
        publication="10.1101/2024.06.13.24308877",
        source="https://www.medrxiv.org/content/medrxiv/early/2024/06/13/2024.06.13.24308877/DC2/embed/media-2.zip?download=true",
        filename=file.path("waterfield","models",paste(proteins,"coefs.csv",sep="-"))),
    file="models.csv",
    row.names=F)




