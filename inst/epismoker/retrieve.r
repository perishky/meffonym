

load(url("https://github.com/sailalithabollepalli/EpiSmokEr/raw/refs/heads/master/data/MLM_coefficients.RData"))

models = sapply(c(current="CS",former="FS",never="NS"), function(group) {
    coefs = get(paste0(group,"_final_coefs"))
    data.frame(
        pred.var = c("intercept",names(coefs)[-1]),
        coef=unname(coefs))
},simplify=F)

for (i in 1:length(models))
    write.csv(models[[i]], file=paste0(names(models)[i],"-coefs.csv"),row.names=F)
