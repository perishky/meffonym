## Hillary https://www.medrxiv.org/content/10.1101/2023.11.02.23298000v1
## https://raw.githubusercontent.com/robertfhillary/dnam-crp/main/Projections_Test_Cohorts/Predictors_by_groups.csv
url <- "https://raw.githubusercontent.com/robertfhillary/dnam-crp/main/Projections_Test_Cohorts/Predictors_by_groups.csv"
hillary.raw <- read.csv(url,
				stringsAsFactors=F)
elnet <- subset(hillary.raw, 
			Predictor == "Elnet", 
			select = c(CpG, Beta))

colnames(elnet) <- c("pred.var","coef")
write.csv(elnet, file="coefs.csv", row.names=F)
