grimage.link <- "https://github.com/bio-learn/biolearn/raw/refs/heads/master/biolearn/data/GrimAgeV1.csv"
grimagev2.link <- "https://github.com/bio-learn/biolearn/raw/refs/heads/master/biolearn/data/GrimAgeV2.csv"

grimage <- read.csv(grimage.link)
grimagev2 <- read.csv(grimagev2.link)

calculate.grimage <- function(meth,female,age,coefs) {
    dat = cbind(t(meth), Female=female, Age=age)

    transform = coefs[coefs$Y.pred == "transform",]
    transform = setNames(transform$beta, transform$var)
    
    cox = coefs[coefs$Y.pred == "COX",]
    cox = setNames(cox$beta, cox$var)

    subscorenames = setdiff(coefs$Y.pred, c("transform","COX"))

    subscores = sapply(subscorenames, function(name) {
        coefs = coefs[coefs$Y.pred == name & (coefs$var=="Intercept" | coefs$var %in% colnames(dat)),]
         coefs = setNames(coefs$beta, coefs$var)
         intercept = coefs['Intercept']
         coefs = coefs[names(coefs) != "Intercept"]
         intercept + dat[,names(coefs)] %*% matrix(coefs,ncol=1)
    })
    subscores = cbind(subscores, dat[,c("Female","Age")])
    scores = subscores[,names(cox)] %*% matrix(cox, ncol=1)
    std_scores = (scores-transform["m_cox"])/transform["sd_cox"]
    std_scores*transform["sd_age"] + transform["m_age"]
}

sites = unique(c(grimage$var, grimagev2$var))
sites = sites[grep("^cg",sites)]

meth0 = matrix(0,ncol=length(sites),nrow=length(sites))
rownames(meth0) = sites
meth = meth0
diag(meth) = 1

models <- lapply(list("grimage"=grimage,"grimagev2"=grimagev2), function(grimage) {
    intercept=calculate.grimage(meth0,female=0,age=0,grimage)[1]
    Female=calculate.grimage(meth0,female=1,age=0,grimage)[1]-intercept
    Age=calculate.grimage(meth0,female=0,age=1,grimage)[1]-intercept
    cg=calculate.grimage(meth,female=0,age=0,grimage)-intercept
    names(cg)=rownames(meth)
    c(intercept=intercept,female=Female,age=Age,cg)
})

subs <- lapply(list("grimage"=grimage,"grimagev2"=grimagev2), function(coefs) {
    sapply(setdiff(coefs$Y.pred, c("transform","COX")), function(subname) {
        coefs <- coefs[coefs$Y.pred == subname,]
        coefs <- setNames(coefs$beta, coefs$var)
        if ("Intercept" %in% names(coefs))
            names(coefs)[names(coefs) == "Intercept"] <- "intercept"
        coefs
    }, simplify=F)
})

for (name in names(subs)) {
    names(subs[[name]]) <- paste(name, names(subs[[name]]), sep="-")
    models <- c(models, subs[[name]])
}

for (name in names(models)) 
    write.csv(
        data.frame(pred.var=names(models[[name]]), coef=models[[name]])
        file=paste(name, "coefs.csv", sep="-"),
        row.names=F)

