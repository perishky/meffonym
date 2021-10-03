## Test all functions, particularly DNA methylation age estimates

library(meffonym)

## Load an example dataset from Horvath's clock website.
url <- "https://dnamage.genetics.ucla.edu/sites/all/files/tutorials"
dnam.url <- file.path(url, "MethylationDataExample55.csv")
sample.url <- file.path(url, "SampleAnnotationExample55.csv")
data <- list(dnam=as.matrix(read.csv(dnam.url, row.names=1)),
             samples=read.csv(sample.url))

## Load age estimates obtained using the Horvath's clock website.
website <- read.csv("MethylationDataExample55.output.csv")

## Calculate DNAmAge using Steve Horvath's scripts
source("horvath-example.r") ## creates horvath.example (2 minutes)

## Normalize the data to Horvath's 'gold standard' using BMIQ. (2 minutes)
data$norm <- meffonym.bmiq.calibration(data$dnam, meffonym.horvath.standard())

## Estimate age.
norm.ret <- meffonym.score(data$norm, "horvath")

## Our normalized estimates should be exactly the same
## as those obtained using the Horvath script
stopifnot(all(abs(norm.ret$score - horvath.example$DNAmAge) < 1e-6))

## The Horvath website appears to use a different random seed
## than the published script 
## when generating random numbers for the BMIQ normalization.
quantile(norm.ret$score - horvath.example$DNAmAge)
quantile(norm.ret$score - website$DNAmAge)
## > quantile(norm.ret$score - horvath.example$DNAmAge)
##            0%           25%           50%           75%          100% 
## -6.394885e-14 -2.131628e-14  7.993606e-15  7.549517e-14  1.421085e-13 
## > quantile(norm.ret$score - website$DNAmAge)
##           0%          25%          50%          75%         100% 
## -0.042524709 -0.018469904 -0.007308543  0.002707157  0.007430722 


## Calculate DNA methylation age acceleration.
acc <- cbind(norm=residuals(lm(norm.ret$score ~ data$samples$Age)),
             example=residuals(lm(horvath.example$DNAmAge ~ data$samples$Age)),
             website=website$AgeAccelerationResidual)

## Our normalized acceleration estimates should be exactly
## the same as those obtained from the Horvath script.
stopifnot(all(abs(acc[,"norm"] - acc[,"example"]) < 1e-6))

cor(acc)
##              norm   example   website
## norm    1.0000000 1.0000000 0.9999948
## example 1.0000000 1.0000000 0.9999948
## website 0.9999948 0.9999948 1.0000000
apply(acc, 2, quantile)
##            norm    example    website
## 0%   -4.0585129 -4.0585129 -4.0779644
## 25%  -1.8685057 -1.8685057 -1.8801605
## 50%  -1.2822465 -1.2822465 -1.2693940
## 75%   0.1959615  0.1959615  0.2043791
## 100% 11.0172491 11.0172491 11.0293354

## We now test the remaining functions.
rowtstats <- function(x, group) {
    stats <- function(x) {
        mean <- rowMeans(x)
        n <- ncol(x)
        var <- rowSums((x-mean)^2)/(n-1)
        list(mean=mean, n=n, var=var)
    }
    s0 <- stats(x[,which(group==min(group))])
    s1 <- stats(x[,which(group==max(group))])
    se <- sqrt(s1$var/s1$n + s0$var/s0$n)
    diff <- s1$mean - s0$mean
    data.frame(diff=diff, se=se, t=diff/se)
}
tstats <- rowtstats(data$norm, data$samples[,"diseaseStatus"])
idx <- order(abs(tstats$t), decreasing=T)[1:50]

meffonym.add.model("cc", rownames(data$norm)[idx], tstats$diff[idx], "example")

stopifnot("cc" %in% meffonym.models())

cc.ret <- meffonym.score(data$norm, "cc",transform=NULL)
....................

stats <- coef(summary(lm(cc.ret$score ~ data$samples[,"diseaseStatus"])))
stats
## > stats 
##                                   Estimate  Std. Error   t value     Pr(>|t|) 
## (Intercept)                     0.07775667 0.002401027 32.384751 1.450308e-14 
## data$samples[, "diseaseStatus"] 0.02785393 0.003395565  8.203031 1.021697e-06 

stopifnot(stats[2,"Pr(>|t|)"] < 10e-5)

scores <- sapply(sort(meffonym.models()), function(model) {
    meffonym.score(data$norm, model)$score
})
cor(scores)
## > cor(scores) 
##               hannum horvath          cc bohlin.1se bohlin.min 
## hannum     1.0000000   0.9223671 -0.24629895   -0.76579610   -0.81488944 
## horvath    0.9223671   1.0000000 -0.26196565   -0.57994926   -0.61587929 
## cc            -0.2462989  -0.2619657  1.00000000   -0.08397825   -0.02775553 
## bohlin.1se -0.7657961  -0.5799493 -0.08397825    1.00000000    0.94593012 
## bohlin.min -0.8148894  -0.6158793 -0.02775553    0.94593012    1.00000000 
## knight      0.1657923   0.2402457 -0.47841075    0.10338741    0.20935622 
##                knight 
## hannum     0.1657923 
## horvath    0.2402457 
## cc            -0.4784107 
## bohlin.1se  0.1033874 
## bohlin.min  0.2093562 
## knight      1.0000000 

