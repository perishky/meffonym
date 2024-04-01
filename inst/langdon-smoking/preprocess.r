library(readxl)

url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8597260/bin/13148_2021_1191_MOESM1_ESM.xlsx"
download.file(url, destfile=basename(url))

x <- read_excel(basename(url), skip=1)

models <- c(
    "ever-vs-never-candidate",
    "ever-vs-never-agnostic",
    "current-vs-former-candidate",
    "current-vs-former-agnostic")
start <- sort(
    c(which(x$Classes == "Ever vs never"),
      which(x$Classes == "Current vs former"),
      which(x[["Model name"]] == "Agnostic LASSO")))
x$model <- c(
    rep(models[1],start[2]-start[1]),
    rep(models[2],start[3]-start[2]),
    rep(models[3],start[4]-start[3]),
    rep(models[4],nrow(x)-start[4]+1))

for (model in models) {
    coefs <- x[x$model == model,c("CpG","Beta")]
    colnames(coefs) <- c("pred.var","coef")
    write.csv(coefs, file=paste(model, "coefs.csv", sep="-"),row.names=F)
}
