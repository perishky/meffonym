## Test all functions, particularly DNA methylation age estimates

library(meffonym)

## Load an example dataset from Horvath's clock website.
url <- "https://dnamage.genetics.ucla.edu/sites/all/files/tutorials"
dnam.url <- file.path(url, "MethylationDataExample55.csv")
sample.url <- file.path(url, "SampleAnnotationExample55.csv")
data <- list(dnam=as.matrix(read.csv(dnam.url, row.names=1)),
             samples=read.csv(sample.url))

## Load age estimates obtained using the Horvath's clock website.
website.ret <- read.csv(system.file("MethylationDataExample55.output.csv", package="meffonym"))

## Normalize the data to Horvath's 'gold standard' using BMIQ.
data$norm <- meffonym.bmiq.calibration(data$dnam, meffonym.horvath.standard()) ## 2 minutes

## Estimate age.
raw.ret <- meffonym.score(data$dnam, "age.horvath")
norm.ret <- meffonym.score(data$norm, "age.horvath")

## Our normalized estimates should be exactly the same as those obtained from the clock website.
stopifnot(all(abs(norm.ret$score - website.ret$DNAmAge) < 1e-6))

## Normalization is clearly necessary to exactly match clock website output.
est <- cbind(raw=raw.ret$score, norm=norm.ret$score, website=website.ret$DNAmAge)
cor(est)
apply(est, 2, quantile)
## > cor(est) 
##               raw      norm   website 
## raw     1.0000000 0.9986507 0.9986507 
## norm    0.9986507 1.0000000 1.0000000 
## website 0.9986507 1.0000000 1.0000000 
## > apply(est, 2, quantile) 
##            raw       norm    website 
## 0%    1.375200  0.9844838  0.9844838 
## 25%   9.008508  7.6804399  7.6804399 
## 50%  25.219689 26.0938764 26.0938764 
## 75%  35.676135 39.5862988 39.5862988 
## 100% 55.990342 62.1743800 62.1743800 

## Calculate DNA methylation age acceleration.
acc <- cbind(raw=residuals(lm(raw.ret$score ~ data$samples$Age)),
             norm=residuals(lm(norm.ret$score ~ data$samples$Age)),
             website=website.ret$AgeAccelerationResidual)

## Our normalized acceleration estimates should be exactly
## the same as those obtained from the clock website.
stopifnot(all(abs(acc[,"norm"] - acc[,"website"]) < 1e-6))


cor(acc)
apply(acc, 2, quantile)
## > cor(acc) 
##               raw      norm   website 
## raw     1.0000000 0.9673395 0.9673395 
## norm    0.9673395 1.0000000 1.0000000 
## website 0.9673395 1.0000000 1.0000000 
## > apply(acc, 2, quantile) 
##             raw       norm    website 
## 0%   -4.4520747 -4.0779644 -4.0779644 
## 25%  -1.8575085 -1.8801605 -1.8801605 
## 50%  -0.4587683 -1.2693940 -1.2693940 
## 75%   0.3506519  0.2043791  0.2043791 
## 100%  8.9151850 11.0293354 11.0293354 

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

cc.ret <- meffonym.score(data$norm, "cc")

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
##               age.hannum age.horvath          cc ga.bohlin.1se ga.bohlin.min 
## age.hannum     1.0000000   0.9223671 -0.24629895   -0.76579610   -0.81488944 
## age.horvath    0.9223671   1.0000000 -0.26196565   -0.57994926   -0.61587929 
## cc            -0.2462989  -0.2619657  1.00000000   -0.08397825   -0.02775553 
## ga.bohlin.1se -0.7657961  -0.5799493 -0.08397825    1.00000000    0.94593012 
## ga.bohlin.min -0.8148894  -0.6158793 -0.02775553    0.94593012    1.00000000 
## ga.knight      0.1657923   0.2402457 -0.47841075    0.10338741    0.20935622 
##                ga.knight 
## age.hannum     0.1657923 
## age.horvath    0.2402457 
## cc            -0.4784107 
## ga.bohlin.1se  0.1033874 
## ga.bohlin.min  0.2093562 
## ga.knight      1.0000000 

