## Gadd episcores for GDF15 and nt-probnp 
#https://static-content.springer.com/esm/art%3A10.1038%2Fs43587-023-00391-4/MediaObjects/43587_2023_391_MOESM3_ESM.xlsx
library(readxl)
url <-"https://www.medrxiv.org/content/medrxiv/early/2023/10/19/2023.10.18.23297200/DC3/embed/media-3.xlsx"
download.file(url, destfile=basename(url))
if (!file.exists(basename(url))) download.file(url, destfile=basename(url))

scores <- read_excel(basename(url), 
                sheet = "Supplementary Table 5",
                skip = 3)

scores <- scores[,c("CpG", "Coefficient", "Predictor")]
colnames(scores) <- c("pred.var","coef", "target")

# write file for each protein target
score.list <- split(scores, scores$target)
score.list <- sapply(score.list, 
			function(i) i[,-which(colnames(i)=="target")], simplify=F)
names(score.list) <- tolower(names(score.list))

sapply(names(score.list), function(i){
	write.csv(score.list[[i]], 
		paste0(i, ".csv"), row.names=F) 
}, simplify=F)

# assemble meta data for adding to models.csv
models <- data.frame(name = names(score.list),
					tissue = "blood",
					target = unique(scores$target),
					publication = "10.1101/2023.10.18.23297200",
					source = "https://www.medrxiv.org/content/medrxiv/early/2023/10/19/2023.10.18.23297200/DC3/embed/media-3.xlsx",
					filename = file.path("gadd-biomarkers", paste0(names(score.list), ".csv")))

write.csv(models, "gadd-biomarkers-models.csv", row.names=F, quote = F)