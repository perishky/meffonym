## Compare age estimates against Horvath's age calculator

## save methylation dataset
if (file.exists("GSE145254-meth.csv"))
    source("save-GSE145254.r")

library(datatable)
data <- list(
    dnam=as.matrix(fread("GSE145254-meth.csv"),rownames=1),
    samples=as.data.frame(fread("GSE145254-samples.csv")))

dat1 <- data.frame(ProbeID=rownames(data$dnam),data$dnam)
source("horvath-example.r")
## in: dat1
## out: horvath.example

calculator.clocks <- read.csv(text="
meffonym,calculator,accel
horvath,DNAmAge,AgeAccelerationResidual
hannum,DNAmAgeHannum,AgeAccelerationResidualHannum
phenoage,DNAmPhenoAge,AgeAccelPheno
skin,DNAmAgeSkinBloodClock,
grimage,DNAmGrimAgeBasedOnRealAge,DNAmGrimAgeBasedOnRealAgeAdjAge
grimagev2,DNAmGrimAge2BasedOnRealAge,
dnamtl,DNAmTL,DNAmTLAdjAge", stringsAsFactors=F)

## https://dnamage.clockfoundation.org/
website <- read.csv("GSE145254-dnamage-calculator.csv")

library(meffonym)
clock.data <- rbind(data$dnam,female=data$samples$Female,age=data$samples$Age)
scores <- sapply(
    calculator.clocks$meffonym,
    function(model) {
        cat(date(), model, "\n")
        meffonym.score(clock.data, model)$score
    })

data$norm <- meffonym.bmiq.calibration(data$dnam,meffonym.horvath.standard()) ## 2 minutes

scores[,"horvath"] <- meffonym.score(data$norm, "horvath")$score

cor(cbind(meffonym=scores[,"horvath"],
          script=horvath.example[,"DNAmAge"],
          website=website[,"DNAmAge"],
          age=data$samples$Age))
##           meffonym    script   website       age
## meffonym 1.0000000 0.9979165 0.9986968 0.8437749
## script   0.9979165 1.0000000 0.9997274 0.8404264
## website  0.9986968 0.9997274 1.0000000 0.8394921
## age      0.8437749 0.8404264 0.8394921 1.0000000

quantile(horvath.example[,"DNAmAge"]-scores[,"horvath"])
quantile(horvath.example[,"DNAmAge"]-website[,"DNAmAge"])
## > quantile(horvath.example[,"DNAmAge"]-scores[,"horvath"])
##          0%         25%         50%         75%        100% 
## -1.47335212 -0.27982687 -0.13511150  0.05374547  0.87773815 
## > quantile(horvath.example[,"DNAmAge"]-website[,"DNAmAge"])
##       0%      25%      50%      75%     100% 
## 5.401105 5.893641 5.961252 6.016499 6.179862 

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
##     dnamtl  grimagev2    grimage   phenoage     hannum    horvath       skin 
## -0.4699341  0.6096642  0.7235472  0.7535686  0.8069383  0.8437749  0.9104319 

r <- cor(website[,calculator.clocks$calculator],
         scores[,calculator.clocks$meffonym],
         use="p")
data.frame(name=colnames(r),r=diag(r))
##        name         r
## 1   horvath 0.9986968
## 2    hannum 1.0000000
## 3  phenoage 1.0000000
## 4      skin 1.0000000
## 5   grimage 1.0000000
## 6 grimagev2 1.0000000
## 7    dnamtl 1.0000000

stopifnot(all(diag(r) > 0.99))

scores.adj <- apply(scores,2,function(ests) {
    residuals(lm(ests ~ data$samples$Age, na.action=na.exclude))
})

accel.exists <- calculator.clocks$accel != ""
r <- cor(website[,calculator.clocks$accel[accel.exists]],
         scores.adj[,calculator.clocks$meffonym[accel.exists]],
         use="p")
data.frame(name=colnames(r),r=diag(r))
##       name        r
## 1  horvath 0.995639
## 2   hannum 1.000000
## 3 phenoage 1.000000
## 4  grimage 1.000000
## 5   dnamtl 1.000000

stopifnot(all(diag(r) > 0.99))

diff <- (website[,calculator.clocks$calculator]
         - scores[,calculator.clocks$meffonym])

colnames(diff) <- calculator.clocks$meffonym

apply(abs(diff), 2, quantile)
##       horvath       hannum   phenoage         skin    grimage    grimagev2
## 0%   5.231281 0.000000e+00 0.00012748 0.000000e+00 0.01222936 1.421085e-14
## 25%  5.925395 7.105427e-15 0.00012748 1.421085e-14 0.01323864 9.237056e-14
## 50%  6.078209 7.105427e-15 0.00012748 2.842171e-14 0.01389040 1.421085e-13
## 75%  6.226368 1.421085e-14 0.00012748 4.263256e-14 0.01423449 1.882938e-13
## 100% 6.874457 2.131628e-14 0.00012748 5.684342e-14 0.01496891 2.557954e-13
##            dnamtl
## 0%   0.000000e+00
## 25%  8.881784e-16
## 50%  1.776357e-15
## 75%  3.552714e-15
## 100% 5.329071e-15










