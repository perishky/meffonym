printFlush <- print
library(impute)

url <- "https://horvath.genetics.ucla.edu/html/dnamage"

probeAnnotation21kdatMethUsed <- read.csv(file.path(url,"probeAnnotation21kdatMethUsed.csv"))
                                          
probeAnnotation27k <- read.csv(file.path(url,"datMiniAnnotation27k.csv"))

datClock <- read.csv(file.path(url,"AdditionalFile3.csv"))

source(file.path(url, "NORMALIZATION.R"))

dat1 <- read.csv(file.path(url, "MethylationDataExample55.csv"))
nSamples <- ncol(dat1)-1
normalizeData <- TRUE

trafo= function(x,adult.age=20) {
    x=(x+1)/(1+adult.age);
    ifelse(x<=1,log(x),x-1);
}

anti.trafo= function(x,adult.age=20) {
    ifelse(x<0,
           (1+adult.age)*exp(x)-1,
           (1+adult.age)*x+adult.age)
}

XchromosomalCpGs <- as.character(probeAnnotation27k$Name[probeAnnotation27k$Chr=="X"])
selectXchromosome <- is.element(dat1[,1], XchromosomalCpGs )
selectXchromosome[is.na(selectXchromosome)] <- FALSE
meanXchromosome <- rep(NA, dim(dat1)[[2]]-1)
meanXchromosome= as.numeric(apply( as.matrix(dat1[selectXchromosome,-1]),2,mean,na.rm=TRUE)) 

dat1 <- dat1[match(probeAnnotation21kdatMethUsed$Name,
                   dat1$ProbeID),]

set.seed(1)
source(file.path(url, "StepwiseAnalysis.txt"))

horvath.example <- datout

rm(anti.trafo,
   betaEst2,
   blc2,
   BMIQ,
   BMIQcalibration,
   CalibrateUnitInterval,
   CheckBMIQ,
   Comment,
   dat1,
   datClock,
   datMethClock,
   datMethClock0,
   datMethUsedNormalized,
   datout,
   fastImputation,
   lab1,
   maxMethBySample,
   meanMethBySample,
   meanXchromosome,
   minMethBySample,
   noMissingPerSample,
   normalizeData,
   nSamples,
   predictedAge,
   predictedGender,
   printFlush,
   probeAnnotation21kdatMethUsed,
   probeAnnotation27k,
   restSamples,
   selectCpGsClock,
   selectXchromosome,
   trafo,
   url,
   website.ret,
   XchromosomalCpGs)
