# DNA methylation indices of exposure and phenotype (meffonym)

## Installation 
```
library(devtools)
install_github("perishky/meffonym")
```

## Typical use

In the example below, we
show how to apply the Hannum clock (Hannum et al. 2013).

```
library(meffonym)

meth <- ... ## methylation matrix, rows=CpG site, columns=samples
age <- ... ## age for each sample

ret <- meffonym.score(meth, "hannum")

cor(ret$score, age)
## [1] 0.9342
```

## Current list of models

|name          |tissue            |target              |publication                        |source                                                                                                                |filename                |
|:-------------|:-----------------|:-------------------|:----------------------------------|:---------------------------------------------------------------------------------------------------------------------|:-----------------------|
|knight        |cord              |gestational age     |10.1186/s13059-016-1068-z          |https://static-content.springer.com/esm/art:10.1186/s13059-016-1068-z/MediaObjects/13059_2016_1068_MOESM3_ESM.csv     |knight/coefs.csv        |
|horvath       |multi             |age                 |10.1186/10.1186/gb-2013-14-10-r115 |https://horvath.genetics.ucla.edu/html/dnamage/AdditionalFile3.csv                                                    |horvath/coefs.csv       |
|hannum        |blood             |age                 |10.1016/j.molcel.2012.10.016       |https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3780611/bin/NIHMS418935-supplement-02.xlsx                               |hannum/coefs.csv        |
|bohlin.min    |cord              |gestational age     |10.1186/s13059-016-1063-4          |http://dx.doi.org/10.5281/zenodo.60498                                                                                |bohlin/model-min.csv    |
|bohlin.1se    |cord              |gestational age     |10.1186/s13059-016-1063-4          |http://dx.doi.org/10.5281/zenodo.60498                                                                                |bohlin/model-1se.csv    |
|epitoc        |blood             |stem cell divisions |10.1186/s13059-016-1064-3          |https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1064-3/MediaObjects/13059_2016_1064_MOESM2_ESM.xls |epitoc/coefs.csv        |
|dunedinpoam38 |blood             |pace of aging       |10.7554/eLife.54870                |https://github.com/danbelsky/DunedinPoAm38                                                                            |dunedinpoam38/coefs.csv |
|brenner       |blood             |mortality           |10.1038/ncomms14617                |Supplementary Figure 1                                                                                                |brenner/coefs.csv       |
|epitoc2       |blood             |stem cell divisions |10.1186/s13073-020-00752-3         |https://zenodo.org/record/2632938/files/dataETOC2.Rd?download=1                                                       |epitoc2/coefs.csv       |
|lee.cpc       |placenta          |gestational age     |10.18632/aging.102049              |https://www.aging-us.com/article/102049/supplementary/SD2/0/aging-v11i12-102049-supplementary-material-SD2.csv        |lee/coefs-cpc.csv       |
|lee.rpc       |placenta          |gestational age     |10.18632/aging.102049              |https://www.aging-us.com/article/102049/supplementary/SD2/0/aging-v11i12-102049-supplementary-material-SD2.csv        |lee/coefs-rpc.csv       |
|lee.rrpc      |placenta          |gestational age     |10.18632/aging.102049              |https://www.aging-us.com/article/102049/supplementary/SD2/0/aging-v11i12-102049-supplementary-material-SD2.csv        |lee/coefs-rrpc.csv      |
|phenoage      |blood             |age                 |10.18632/aging.101414              |https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5940111/bin/aging-10-101414-s002.csv                                     |levine/coefs.csv        |
|mayne         |placenta          |gestational age     |10.2217/epi-2016-0103              |https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6040051/bin/epi-09-279-s6.xlsx                                           |mayne/coefs.csv         |
|pedbe         |buccal epithelial |age                 |10.1073/pnas.1820843116            |https://raw.githubusercontent.com/kobor-lab/Public-Scripts/master/datcoefInteresting94.csv                            |pedbe/coefs.csv         |
|skin          |fibroblasts       |age                 |10.18632/aging.101508              |https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6075434/bin/aging-10-101508-s005.csv                                     |skin/coefs.csv          |
|dnamtl        |blood             |telomere length     |10.18632/aging.102173              |https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6738410/bin/aging-11-102173-s003.xlsx                                    |tl/coefs.csv            |
|zhang         |blood             |age                 |10.1186/s13073-019-0667-1          |https://github.com/qzhang314/DNAm-based-age-predictor                                                                 |zhang/coef-en.csv       |
|zhang.blup    |blood             |age                 |10.1186/s13073-019-0667-1          |https://github.com/qzhang314/DNAm-based-age-predictor                                                                 |zhang/coef-blup.csv     |
|miage         |blood             |stem cell divisions |10.1080/15592294.2017.1389361      |http://www.columbia.edu/~sw2206/softwares/mitotic_age_R_code.zip                                                      |miage/coefs.csv         |

