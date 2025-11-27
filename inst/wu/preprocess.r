library(readxl)
url <-"https://www.aging-us.com/article/102399/supplementary/SD3/0/aging-v11i22-102399-supplementary-material-SD3.xlsx"
download.file(url, destfile=basename(url))
if (!file.exists(basename(url))) download.file(url, destfile=basename(url))

wu <- read_excel(basename(url))

wu <- wu[,c("CpGs", "Active.coefficients")]
colnames(wu) <- c("pred.var","coef")

write.csv(wu, file="coefs.csv", row.names=F)

## assemble meta data for adding to models.csv
models <- data.frame(name = "wu",
					tissue = "blood",
					target = "age",
					publication = "10.1038/s43587-023-00557-0",
					source = url,
					filename = file.path("wu", "coefs.csv"))

## save the meta data to append to models.csv
write.csv(models, "wu-model.csv", row.names=F, quote = F)

