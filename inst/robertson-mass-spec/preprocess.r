library(data.table)

## download models
url = "https://zenodo.org/records/16924748/files/training_effect_sizes.csv.gz?download=1"
download.file(url, destfile="training_effect_sizes.csv.gz")

## load models (ignore anything with  < 0.05) 
stats = fread("training_effect_sizes.csv.gz",data.table=F)
stats = stats[stats$pip > 0.05,]
stats = stats[,c("CpGs","effect_size", "protein")]
colnames(stats) = c("pred.var","coef", "target")

## split by target (protein)
models = split(stats, stats$target)
models = sapply(models, function(i) i[,c("pred.var","coef")], simplify=F)
names(models) = tolower(names(models))

## save each model to a csv file
dir.create("models")
for (name in names(models))
  fwrite(models[[name]],file.path("models", paste0(name, ".csv")), row.names=F)

## assemble metadata 
metadata = data.frame(
  name = paste("robertson", names(models), sep="-"),
  tissue = "blood",
  target = names(models),
  publication = "10.1186/s13059-025-03892-0",
  source = "https://zenodo.org/records/16924748/files/training_effect_sizes.csv.gz?download=1",
  filename = file.path("robertson-mass-spec", "models", paste0(names(models), ".csv")))
write.csv(metadata, "metadata.csv", row.names=F, quote = F)
