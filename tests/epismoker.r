filename = "epismoker-results.csv"
if (!file.exists(filename)) {
    ## remotes::install_github("sailalithabollepalli/EpiSmokEr")
    library(EpiSmokEr)
    data("dummyBetaData")
    samplesheet = read.csv(system.file("extdata", "samplesheet_GSE42861.csv", package="EpiSmokEr"), header=TRUE, sep=",")
    eresults = epismoker(dataset=dummyBetaData, samplesheet = samplesheet, method = "SSt")
    write.csv(eresults,file=filename,row.names=F)
} else
    eresults = read.csv(filename)

load(url("https://github.com/sailalithabollepalli/EpiSmokEr/raw/refs/heads/master/data/dummyBetaData.rda"))

samplesheet = read.csv("https://github.com/sailalithabollepalli/EpiSmokEr/raw/refs/heads/master/inst/extdata/samplesheet_GSE42861.csv")

library(meffonym)
mresults = meffonym.epismoker(dummyBetaData,samplesheet$Gender == "f")

stopifnot(
    max(eresults[,c("probs_CS","probs_FS","probs_NS")]
        -mresults[,c("current","former","never")])
    < 1e-12)
          
