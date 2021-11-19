library(meffonym)

## load example dataset from the estimage website
url <- "https://estimage.iac.rm.cnr.it/EstimAge/static/examples/"
dnam.url <- file.path(url, "GSE72776_reduced.data.csv.gz")
sample.url <- file.path(url, "GSE72776_reduced.metadata.csv")

data <- list(dnam=read.csv(dnam.url,row.names=1),
             samples=read.csv(sample.url,row.names=1))
data$dnam <- as.matrix(data$dnam)

## estimage.iac.rm.cnr.it with imputation by 'knn', no cell count adjustment
estimage <- read.csv("Table_Age.2021_11_19_10_35_40.csv",stringsAsFactors=F,check.names=F)

estimage.clocks <- read.csv(comment.char="#",text="
meffonym,estimage
#CorticalClock (y)#Shireby et al., 2020
#EPM_0.8 (y)# Snir et al., 2019
#FAradas16 (y)#Freire-Aradas et al., 2016
hannum,Hannum13 (y)
horvath,Horvath13 (y)
skin,Horvath18 (y)
#MEAT (y)#Voisin et al., 2020
pedbe,PedBE (y)
phenoage,PhenoAge (y)
zhang,Zhang19.enpred (y)
zhang.blup,Zhang19.blupred (y)
#ZPiekarska15 (y)#Zbieć-Piekarska et al., 2015
bohlin.min,Bohlin16 (w)
knight,Knight16 (w)
lee.rpc,Lee19.RPC (w)
lee.cpc,Lee19.CPC (w)
lee.rrpc,Lee19.refRPC (w)
mayne,Mayne17 (w)
epitoc,epiTOC (cc)
miage,MiAge (cc)
dnamtl,DNAmTL (kb)",stringsAsFactors=F)

setdiff(meffonym.models(), estimage.clocks$meffonym)
## [1] "bohlin.1se"    "brenner"       "dunedinpoam38" "epitoc2"

scores <- sapply(
    estimage.clocks$meffonym,
    function(model) {
        cat(date(), model, " ")
        ret <- meffonym.score(data$dnam, model)
        cat(" used ", length(ret$sites), "/", length(ret$vars), "sites\n")
        ret$score
    })

data$norm <- meffonym.bmiq.calibration(data$dnam, meffonym.horvath.standard())
scores[,"horvath"] <- meffonym.score(data$norm, "horvath")$score

scores[,"bohlin.min"] <- scores[,"bohlin.min"]/52
scores[,"dnamtl"] <- 20 + 21*scores[,"dnamtl"]
scores[,"knight"] <- meffonym.score(data$norm, "knight")$score
scores[,"knight"] <- 20 + 21*scores[,"knight"]
scores[,"zhang"] <- meffonym.score(data$dnam, "zhang", scale=F)$score
scores[,"zhang"] <- -112.373 + 2.663*scores[,"zhang"]
scores[,"zhang.blup"] <- meffonym.score(data$dnam, "zhang.blup", scale=F)$score
scores[,"zhang.blup"] <- -139.866 + 2.634*scores[,"zhang.blup"]

r <- diag(cor(scores, estimage[,estimage.clocks$estimage]))

stopifnot(all(r[!grepl("lee",colnames(scores))] > 0.99))

data.frame(
    name=colnames(scores),
    r=r)
##          name          r
## 1      hannum  1.0000000
## 2     horvath  1.0000000
## 3        skin  1.0000000
## 4       pedbe  1.0000000
## 5    phenoage  1.0000000
## 6       zhang  0.9998752
## 7  zhang.blup  0.9990261
## 8  bohlin.min  1.0000000
## 9      knight  1.0000000
## 10    lee.rpc -0.1414252
## 11    lee.cpc -0.1117716
## 12   lee.rrpc -0.2649742
## 13      mayne  1.0000000
## 14     epitoc  1.0000000
## 15      miage  0.9948304
## 16     dnamtl  1.0000000

t(apply(
    scores - estimage[,estimage.clocks$estimage],
    2,
    quantile, probs=c(0,0.5,1)))
##                                0%           50%         100%
## Hannum13 (y)        -5.684342e-14  1.065814e-14 5.684342e-14
## Horvath13 (y)       -1.136868e-13  0.000000e+00 8.526513e-14
## Horvath18 (y)       -7.105427e-14  7.105427e-15 7.105427e-14
## PedBE (y)           -4.263256e-14  2.664535e-15 4.085621e-14
## PhenoAge (y)        -9.947598e-14  1.421085e-14 7.105427e-14
## Zhang19.enpred (y)  -5.518230e-01  7.650258e-03 2.499857e-01
## Zhang19.blupred (y) -7.117101e-01  6.975951e-02 1.484786e+00
## Bohlin16 (w)        -6.217249e-15  9.769963e-15 2.664535e-14
## Knight16 (w)        -3.865352e-12  0.000000e+00 4.774847e-12
## Lee19.RPC (w)        7.656738e+00  1.515305e+01 2.143021e+01
## Lee19.CPC (w)        7.248691e+00  1.243348e+01 1.724373e+01
## Lee19.refRPC (w)     6.568341e+00  1.075734e+01 1.318655e+01
## Mayne17 (w)         -6.394885e-14  3.552714e-15 6.394885e-14
## epiTOC (cc)         -4.163336e-16  0.000000e+00 3.608225e-16
## MiAge (cc)          -1.976493e+02 -1.196336e+02 2.872462e+02
## DNAmTL (kb)         -4.831691e-13 -9.947598e-14 4.831691e-13

#library(devtools)
#install_github("isglobal-brge/methylclock")
library(methylclock)

mc <- DNAmGA(data$dnam, cell.count=F)

stopifnot(
    all(c(
        cor(meffonym.score(data$dnam, "knight")$score, mc$Knight, use="p"),
        cor(meffonym.score(data$dnam, "bohlin.1se")$score, mc$Bohlin, use="p"),
        cor(scores[,"mayne"], mc$Mayne, use="p"),
        cor(scores[,"lee.cpc"], mc$Lee.CPC, use="p"),
        cor(scores[,"lee.rpc"], mc$Lee.RPC, use="p"),
        cor(scores[,"lee.rrpc"], mc$Lee.refRPC, use="p"))
        >= 1)
    )
