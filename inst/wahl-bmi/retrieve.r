## Wahl et al. Supplementary Table 3
## European = discovery Europeans
## Indian = discovery Indian Asians
## Combined = discovery Combined

url = "https://pmc.ncbi.nlm.nih.gov/articles/instance/5570525/bin/NIHMS70428-supplement-Supplementary_Tables.xlsx"

filename = basename(url)
download.file(url, filename)

library(readxl)

stats = read_xlsx(filename, "ST-3", skip=2)
stats = as.data.frame(stats)
stats = stats[stats$Replicate == "Yes" & stats[["Sentinel marker"]]=="Yes",]

nrow(stats)
## [1] 187

models = c("european","indian","combined")
for (i in 1:length(models)) {
    effect.col = colnames(stats)[grep("Effect",colnames(stats))[i]]
    model = stats[,c("CpG",effect.col)]
    colnames(model) = c("pred.var","coef")
    model$coef = as.double(sub("([^ ]+) .*", "\\1", model$coef))
    write.csv(model, file=paste0(models[i],".csv"), row.names=F)
}


