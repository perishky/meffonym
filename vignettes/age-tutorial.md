# DNA methylation age and acceleration



## Dataset
We will analyze data used to derive the Hannum et al.
DNA methylation age predictor.

> Hannum G, Guinney J, Zhao L, Zhang L et al. Genome-wide methylation
> profiles reveal quantitative views of human aging rates. Mol Cell 2013
> Jan 24;49(2):359-67. PMID: 23177740

Retrieve the data from the Gene Expression Omnibus (GEO) website
(http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279).


```r
source(system.file("read-gse-matrix-file.r", package="meffonym"))

file <-  "GSE40279_series_matrix.txt.gz"
gse <- "GSE40279"
url <- file.path("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn", gse,
                 "matrix", file)

if (!file.exists(file))
    download.file(url, destfile=file) ## 1.2Gb file

if (!exists("hannum"))
    hannum <- read.gse.matrix.file(file) ## ~15 minutes
```

Extract phenotypes/exposures.

```r
extract.characteristic <- function(ch) {
    sub("[^:]+: (.*)", "\\1", ch)
}
hannum$samples$age <- as.numeric(extract.characteristic(hannum$samples$characteristics_ch1))
hannum$samples$gender <- extract.characteristic(hannum$samples$characteristics_ch1.3)
hannum$samples$ethnicity <- extract.characteristic(hannum$samples$characteristics_ch1.4)
```

## Calculate scores for all available models


```r
library(meffonym)

scores <- cbind(age=hannum$samples$age,
                sapply(meffonym.models(), function(model) {
                    meffonym.score(hannum$dnam, model)$score
                }))
```

Calculate correlations between scores.

```r
cor(scores)
```

```
##                         age  age.hannum age.horvath
## age             1.000000000  0.94607431  0.91834357
## age.hannum      0.946074305  1.00000000  0.94750672
## age.horvath     0.918343573  0.94750672  1.00000000
## ga.bohlin.1se  -0.186930069 -0.17344962 -0.18016336
## ga.bohlin.min  -0.237881118 -0.22397250 -0.23879848
## ga.knight      -0.004742566  0.01194531  0.02103212
##                ga.bohlin.1se ga.bohlin.min    ga.knight
## age               -0.1869301    -0.2378811 -0.004742566
## age.hannum        -0.1734496    -0.2239725  0.011945310
## age.horvath       -0.1801634    -0.2387985  0.021032121
## ga.bohlin.1se      1.0000000     0.9440823  0.434868330
## ga.bohlin.min      0.9440823     1.0000000  0.420777076
## ga.knight          0.4348683     0.4207771  1.000000000
```

Calculate differences between biological and estimated age.

```r
apply(scores[,-1], 2, function(estimate) quantile(estimate-scores[,"age"]))
```

```
##      age.hannum age.horvath ga.bohlin.1se ga.bohlin.min
## 0%   -17.953761  -30.427276      176.8630      164.2975
## 25%    1.315673   -5.777963      215.3069      209.8056
## 50%    4.906997   -2.317217      225.8057      220.6148
## 75%    7.870733    1.037909      238.0843      232.8203
## 100%  27.105114   31.887875      281.9626      284.1008
##       ga.knight
## 0%   -55.013400
## 25%  -29.301776
## 50%  -18.833782
## 75%   -8.492379
## 100%  27.836865
```

## Calculate age acceleration

Age accelaration is obtained by adjusting age estimates for actual age.

```r
age.scores <- scores[,c("age.horvath","age.hannum")]
    
acc <- apply(age.scores, 2, function(estimate) {
    residuals(lm(estimate ~ scores[,"age"]))
})
```

Test age acceleration differences between the sexes.
Both Hannum and Horvath scores show
greater acceleration in males than females.

```r
t(sapply(colnames(acc), function(name) {
    acc <- acc[,name]
    coef(summary(lm(acc ~ gender, data=hannum$samples)))["genderM",]
}))
```

```
##             Estimate Std. Error  t value     Pr(>|t|)
## age.horvath 1.546429  0.3908939 3.956134 8.449921e-05
## age.hannum  2.215048  0.3206417 6.908173 1.166350e-11
```

## Obtaining estimates identical to the online Horvath clock

The following estimates are equivalent to that of the Horvath clock.
The calculation can be time consuming, so we will obtain estimates
for only the first 10.

```r
horvath.est <- meffonym.horvath.age(hannum$dnam[,1:10]) ## 1-2 minutes
```

These estimates are obtained in two steps: 
normalizing the methylation data to 'gold standard' used by Horvath when
developing the model, and then applying the model in the normalized data.

```r
standard <- meffonym.horvath.standard()
## step 1: normalization (1-2 minutes)
dnam.norm <- meffonym.bmiq.calibration(hannum$dnam[,1:10], standard) 
## step 2: model application
equiv.est <- meffonym.score(dnam.norm, "age.horvath")
```

They are the same.

```r
cor(equiv.est$score, horvath.est$score)
```

```
## [1] 1
```

```r
equiv.est$score - horvath.est$score
```

```
##  [1] 0 0 0 0 0 0 0 0 0 0
```

Normalization in this case has a fairly minor effect
on the age estimates.

```r
cor(age.scores[1:10,"age.horvath"], horvath.est$score)
```

```
## [1] 0.9996503
```

```r
sort(age.scores[1:10, "age.horvath"] - horvath.est$score)
```

```
##  [1] -1.5863928 -1.4992797 -1.4157137 -1.3949272 -0.9008319 -0.8854542
##  [7] -0.8767604 -0.7215569 -0.7183931 -0.4808198
```

## Calculate scores for a new model


```r
model <- meffonym.get.model("age.horvath")
idx <- order(abs(model$coefficients), decreasing=T)[1:10]
model$coefficients <- model$coefficients[idx]
meffonym.add.model("age.horvath.10",
                     variables=names(model$coefficients),
                     coefficients=model$coefficients,
                     description="Horvath model with only the top 10 effect sizes.")
```

Estimate age with that new model.

```r
ten.est <- meffonym.score(hannum$dnam, "age.horvath.10")
```

The result is obviously similar but 10 CpG sites is obviously
not the same as 353.

```r
cor(ten.est$score, scores[,"age.horvath"])
```

```
## [1] 0.7324617
```

## Document generated from R Markdown


```r
library(knitr)
library(markdown)

options(markdown.HTML.options=union('toc', getOption("markdown.HTML.options")))

knit("age-tutorial.rmd", "age-tutorial.md")
markdownToHTML("age-tutorial.md","age-tutorial.html", stylesheet="style.css")
```
