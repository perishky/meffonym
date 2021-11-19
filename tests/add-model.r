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

## function for applying a t-test to each row of a matrix
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

## simulate disease status variable
set.seed(20211117)
data$samples$disease <- sample(0:1,size=nrow(data$samples),replace=T)

## compare DNAm between disease groups
tstats <- rowtstats(data$dnam, data$samples[,"disease"])

## identify the CpG sites with the top 50 differences
idx <- order(abs(tstats$t), decreasing=T)[1:50]

## create a DNAm model for these top CpG sites
meffonym.add.model(
    name="cc",
    variables=rownames(data$dnam)[idx],
    coefficients=tstats$diff[idx],
    description="example")

## verify that the new model has been added
stopifnot("cc" %in% meffonym.models())

## apply the new model to the methylation data
cc.ret <- meffonym.score(data$dnam, "cc",transform=NULL)

## verify that the resulting model score differentiates by disease status
stats <- coef(summary(lm(cc.ret$score ~ data$samples[,"disease"])))
stats
##                            Estimate  Std. Error   t value     Pr(>|t|)
## (Intercept)               0.9667640 0.009139452 105.77922 1.807144e-20
## data$samples[, "disease"] 0.2506084 0.011798981  21.23983 1.775961e-11

stopifnot(stats[2,"Pr(>|t|)"] < 1e-10)
