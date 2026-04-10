## download and save the Ying models from Albert's biolearn repo:
## 	https://github.com/albert-ying/biolearn

# prep urls and model names
url <- "https://raw.githubusercontent.com/albert-ying/biolearn/b332d123e77ffb21d443a9984332afdb88adbb17/biolearn/data"
model.names <- c("YingAdaptAge", "YingCausAge", "YingDamAge")
model.urls <- sapply(model.names, function(i) file.path(url, paste0(i, ".csv")), simplify=F)

## download model csv files
scores.raw <- sapply(model.urls, function(i) {
	read.csv(i, stringsAsFactors=F)
}, simplify=F)

## coerce to meffonym format
scores <- sapply(scores.raw, function(i) {
	colnames(i) <- c("pred.var", "coef")
	i
	}, simplify=F)

## write out a csv for each model
sapply(names(scores), function(i){
	write.csv(scores[[i]], 
		paste0(i, ".csv"), row.names=F) 
}, simplify=F)

## assemble meta data for adding to models.csv
models <- data.frame(name = model.names,
					tissue = "blood",
					target = gsub("^Ying", "", model.names),
					publication = "10.1038/s43587-023-00557-0",
					source = unlist(model.urls),
					filename = file.path("ying", paste0(names(scores), ".csv")))

## save the meta data to append to models.csv
write.csv(models, "ying-models.csv", row.names=F, quote = F)

