require(Matrix)

options(download.file.method = "wget")
url <- "https://zenodo.org/record/60498/files/predictGA-v1.0.0.zip"
filename <- basename(url)
download.file(url, filename)
filenames <- unzip(filename, list=T)
unzip(filename, files=filenames$Name[grep("sysdata.rda", filenames$Name)], junkpaths=T)
load("sysdata.rda")

## should be equivalent to 'coef(fit, s=lambda)' 
extract.model <- function(fit, s) {
    idx <- which(fit$lambda == s)
    beta <- fit$beta[,idx]
    beta <- beta[which(abs(beta) > 0)]
    intercept <- fit$a0[[idx]]
    c("(Intercept)"=intercept, beta)
}

model.1se <- extract.model(UL.mod.cv$glmnet.fit, s=UL.mod.cv$lambda.1se)
model.min <- extract.model(UL.mod.cv$glmnet.fit, s=UL.mod.cv$lambda.min)

write.csv(data.frame(cpg=names(model.1se), coef=model.1se), row.names=F, file="model-1se.csv")
write.csv(data.frame(cpg=names(model.min), coef=model.min), row.names=F, file="model-min.csv")

unlink(filename)
unlink("sysdata.rda")
