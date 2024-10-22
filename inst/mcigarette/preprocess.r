## Chybowska mCigarette smoking pack years score
library(readxl)
url <-"https://www.medrxiv.org/content/medrxiv/early/2024/05/21/2024.05.21.24307663/DC1/embed/media-1.xlsx"
download.file(url, destfile=basename(url))
if (!file.exists(basename(url))) download.file(url, destfile=basename(url))

mcigarette <- read_excel(basename(url), 
                sheet = "6",
                skip = 2)

mcigarette <- mcigarette[,c("CpG site", "Coefficient")]
colnames(mcigarette) <- c("pred.var","coef")

write.csv(mcigarette, file="coefs.csv", row.names=F)