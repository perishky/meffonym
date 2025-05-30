library(meffonym)

## save methylation dataset
if (!file.exists("GSE145254-meth.csv"))
    source("save-GSE145254.r")

# load the example dataset generated by save-GSE145254.r
library(data.table)
data <- list(
    dnam=as.matrix(fread("GSE145254-meth.csv"),rownames=1),
    samples=as.data.frame(fread("GSE145254-samples.csv")))

str(data)
# List of 2
#  $ dnam   : num [1:866092, 1:23] 0.652 0.92 0.916 0.96 0.965 ...
#   ..- attr(*, "dimnames")=List of 2
#   .. ..$ : chr [1:866092] "cg00000029" "cg00000103" "cg00000109" "cg00000155" ...
#   .. ..$ : chr [1:23] "GSM4310213" "GSM4310214" "GSM4310215" "GSM4310216" ...
#  $ samples:'data.frame':        23 obs. of  7 variables:
#   ..$ Basename : chr [1:23] "201501980023_R01C01" "201501980036_R01C01" "201501970005_R01C01" "201501980023_R02C01" ...
#   ..$ id       : chr [1:23] "GSM4310213" "GSM4310214" "GSM4310215" "GSM4310216" ...
#   ..$ Female   : int [1:23] 1 1 0 1 1 1 1 0 0 1 ...
#   ..$ status   : chr [1:23] "healthy control" "psychiatric case" "healthy control" "healthy control" ...
#   ..$ Age      : int [1:23] 80 60 53 71 69 61 74 52 75 67 ...
#   ..$ ethnicity: chr [1:23] "Black" "Black" "White" "White" ...
#   ..$ Tissue   : chr [1:23] "Whole Blood" "Whole Blood" "Whole Blood" "Whole Blood" ...

# check samples match 
identical(data$samples$id, colnames(data$dnam))

## check a verifiable score
ret <- meffonym.score(data$dnam, "hannum")
cor(ret$score, data$samples$Age)
# [1] 0.8069383

## print the sites not used in the score
missing.inputs <- ret$vars[!ret$vars %in% ret$sites]
missing.inputs
##[1] "intercept"  "cg24079702" "cg14361627" "cg07927379" "cg18473521"
##[6] "cg09651136" "cg21139312"

ret$coefs[names(ret$coefs) %in% missing.inputs]
##cg24079702 cg14361627 cg07927379 cg18473521 cg09651136 cg21139312 
##      2.48      10.70      -1.42       8.85     -15.80      17.10 

## check meffonym.accel
accel <- meffonym.accel(data$dnam, "hannum", age = data$samples$Age)

## check dunedinpace score
packages <- c("preprocessCore") 
lapply(packages, require, character.only=T)
# [1] TRUE

ret <- meffonym.score(data$dnam, "dunedinpace")
dunedin <- meffonym.score(data$dnam, "dunedinpace", calibrate=TRUE)
dunedin <- meffonym.dunedinpace.estimate(data$dnam)

