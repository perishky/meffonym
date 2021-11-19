## Compare age estimates against Horvath's age calculator

library(meffonym)

## Load an example dataset from Horvath's clock website.
url <- "https://horvath.genetics.ucla.edu/html/dnamage/"
dnam.url <- file.path(url, "dat0BloodIllumina450K.zip")
sample.url <- file.path(url, "datSampleBloodIllumina450K.csv")

filename <- basename(dnam.url)
download.file(dnam.url, destfile=filename)
dnam <- read.csv(unz(filename, sub("zip","csv",filename)),row.names=1)

data <- list(dnam=as.matrix(dnam),
             samples=read.csv(sample.url))

## Calculate DNAmAge using Steve Horvath's scripts
source("horvath-example.r") ## creates horvath.example (2 minutes)

## Normalize the data to Horvath's 'gold standard' using BMIQ. (2 minutes)
data$norm <- meffonym.bmiq.calibration(data$dnam, meffonym.horvath.standard())

## Estimate age.
norm.ret <- meffonym.score(data$norm, "horvath")

## Our normalized estimates should be exactly the same
## as those obtained using the Horvath script
stopifnot(all(abs(norm.ret$score - horvath.example$DNAmAge) < 1e-12))

## Load age estimates obtained using the Horvath's clock website.
website <- read.csv("dat0BloodIllumina450K.output.csv",stringsAsFactors=F)

## The Horvath website appears to use a different random seed
## than the published script 
## when generating random numbers for the BMIQ normalization.
quantile(norm.ret$score - horvath.example$DNAmAge)
quantile(norm.ret$score - website$DNAmAge)
## > quantile(norm.ret$score - horvath.example$DNAmAge)
##            0%           25%           50%           75%          100% 
## -1.847411e-13 -4.618528e-14 -2.131628e-14  4.618528e-14  1.847411e-13 
## > quantile(norm.ret$score - website$DNAmAge)
##          0%         25%         50%         75%        100% 
## -0.00301415  0.01991511  0.04975749  0.05892616  0.14197597 

stopifnot(all(abs(norm.ret$score - horvath.example$DNAmAge) < 1e-12))

calculator.clocks <- read.csv(text="
meffonym,calculator,accel
horvath,DNAmAge,AgeAccelerationResidual
hannum,DNAmAgeHannum,DNAmAgeHannumAdjAge
phenoage,DNAmPhenoAge,AgeAccelPheno
skin,DNAmAgeSkinBloodClock,DNAmAgeSkinBloodClockAdjAge
dnamtl,DNAmTL,DNAmTLAdjAge", stringsAsFactors=F)

scores <- sapply(
    calculator.clocks$meffonym,
    function(model) {
        cat(date(), model, "\n")
        meffonym.score(data$dnam, model)$score
    })

scores[,"horvath"] <- meffonym.score(data$norm, "horvath")$score

r <- cor(scores)
r <- data.frame(
    a=rep(colnames(r),ncol(r)),
    b=rep(colnames(r),each=ncol(r)),
    r=as.vector(r),
    stringsAsFactors=F)
r <- r[r$a < r$b,]
r <- r[order(r$r),]

r <- cor(data$samples$Age, scores, use="p")
r[,order(r)]
##     dnamtl   phenoage    horvath     hannum       skin 
## -0.5163266  0.7749623  0.7879994  0.8308086  0.9180756 

r <- cor(website[,calculator.clocks$calculator],
         scores[,calculator.clocks$meffonym],
         use="p")
data.frame(name=colnames(r),r=diag(r))
##       name         r
## 1  horvath 0.9999925
## 2   hannum 1.0000000
## 3 phenoage 0.9999960
## 4     skin 0.9999988
## 5   dnamtl 1.0000000

stopifnot(all(diag(r) > 0.9999))

scores.adj <- apply(scores,2,function(ests) {
    residuals(lm(ests ~ data$samples$Age, na.action=na.exclude))
})

r <- cor(website[,calculator.clocks$accel],
         scores.adj[,calculator.clocks$meffonym],
         use="p")
data.frame(name=colnames(r),r=diag(r))
##       name         r
## 1  horvath 0.9999910
## 2   hannum 1.0000000
## 3 phenoage 0.9999904
## 4     skin 0.9999942
## 5   dnamtl 1.0000000

stopifnot(all(diag(r) > 0.9999))

diff <- (website[,calculator.clocks$calculator]
         - scores[,calculator.clocks$meffonym])

colnames(diff) <- calculator.clocks$meffonym

apply(abs(diff), 2, quantile)
##          horvath       hannum   phenoage         skin       dnamtl
## 0%   0.002352094 0.000000e+00 0.00012748 7.105427e-15 8.881784e-16
## 25%  0.019915105 7.105427e-15 0.00012748 1.421085e-14 1.776357e-15
## 50%  0.049757491 1.421085e-14 0.00012748 2.842171e-14 3.552714e-15
## 75%  0.058926158 1.421085e-14 0.00012748 4.263256e-14 4.440892e-15
## 100% 0.141975969 2.131628e-14 0.09781089 4.303148e-02 5.329071e-15

