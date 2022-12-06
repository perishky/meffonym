
## read raw scores from zenodo v10.0 record
scores.raw <- read.csv("https://zenodo.org/record/5869576/files/Predictors_Shiny_by_Groups.csv")
scores.raw <- scores.raw[, -grep("Mean_Beta_Value", colnames(scores.raw))]
colnames(scores.raw) <- c("pred.var", "coef", "target")

## keep only the protein episcores
predictors <- unique(scores.raw$target)
proteins <- predictors[grep("ADAMTS", predictors):length(predictors)]
scores <- subset(scores.raw, target %in% proteins)

# drop spaces from protein target names
library(stringr)
scores$target <- str_replace(scores$target, " ", "_")

# write file for each protein target
score.list <- split(scores, scores$target)
score.list <- sapply(score.list, 
			function(i) i[,-which(colnames(i)=="target")], simplify=F)

sapply(names(score.list), function(i){
	write.csv(score.list[[i]], 
		paste0(i, ".csv"), row.names=F) 
}, simplify=F)

# assemble meta data for adding to models.csv
models <- data.frame(name = unique(scores$target),
					tissue = "blood",
					target = unique(scores$target),
					publication = "10.7554/ELIFE.71802",
					source = "https://zenodo.org/record/5869576/files/Predictors_Shiny_by_Groups.csv",
					filename = file.path("episcores", paste0(unique(scores$target), ".csv")))

write.csv(models, "episcores-models.csv", row.names=F, quote = F)

