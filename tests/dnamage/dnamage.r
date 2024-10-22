## USE:
## source("dnamage/dnamage.r", chdir=T)
## ret <- dnamage(beta)

## http://labs.genetics.ucla.edu/horvath/dnamage/TUTORIAL1.pdf
## http://labs.genetics.ucla.edu/horvath/dnamage/probeAnnotation21kdatMethUsed.csv
## http://labs.genetics.ucla.edu/horvath/dnamage/datMiniAnnotation27k.csv
## http://labs.genetics.ucla.edu/horvath/dnamage/AdditionalFile3.csv
## http://labs.genetics.ucla.edu/horvath/dnamage/StepwiseAnalysis.txt
## http://labs.genetics.ucla.edu/horvath/dnamage/NORMALIZATION.R

require(RPMM)
require(sqldf) 
## http://cran.r-project.org/src/contrib/Archive/sqldf/sqldf_0.4-6.4.tar.gz
## http://cran.r-project.org/src/contrib/Archive/RSQLite.extfuns/RSQLite.extfuns_0.0.1.tar.gz
require(impute) ## from bioconductor
require(WGCNA)

source("NORMALIZATION.R")

dnamage.probeAnnotation21kdatMethUsed=read.csv("probeAnnotation21kdatMethUsed.csv")
dnamage.probeAnnotation27k=read.csv("datMiniAnnotation27k.csv")
dnamage.datClock=read.csv("AdditionalFile3.csv")

dnamage <- function(x,normalizeData=TRUE) {
    dat0 <- data.frame(id=rownames(x), x)
    nSamples=dim(dat0)[[2]]-1
    nProbes= dim(dat0)[[1]]
    ## the following command may not be needed. But it is sometimes useful when you use read.csv.sql
    dat0[,1]= gsub(x=dat0 [,1],pattern="\"",replacement="")
    cat(paste( "The methylation data set contains", nSamples, "samples (e.g. arrays) and ",
              nProbes, " probes.\n"))
    if (nSamples==0) {
        stop("ERROR: There must be a data input error since there seem to be no samples.\n Make sure that you input a comma delimited file (.csv file)\n that can be read using the R command read.csv.sql(). Samples correspond to columns in that file .")
    }
    if (nProbes==0) {
        stop("ERROR: There must be a data input error since there seem to be zero probes.\n Make sure that you input a comma delimited file (.csv file)\n that can be read using the R command read.csv.sql CpGs correspond to rows.")
    }
    if ( nSamples > nProbes ) {
        warning("MAJOR WARNING: It worries me a lot that there are more samples than CpG probes.\n Make sure that probes correspond to rows and samples to columns.\n I wonder whether you want to first transpose the data and then resubmit them? In any event, I will proceed with the analysis.")
    }
    if ( is.numeric(dat0[,1]) ) {
        stop(paste("Error: The first column does not seem to contain probe identifiers (cg numbers from Illumina) since these entries are numeric values. Make sure that the first column of the file contains probe identifiers such as cg00000292. Instead it contains ", dat0[1:3,1]))
    }
    if ( !is.character(dat0[,1]) ) {
        warning(paste("Major Warning: The first column does not seem to contain probe identifiers (cg numbers from Illumina) since these entries are numeric values. Make sure that the first column of the file contains CpG probe identifiers such as cg00000292. Instead it contains ", dat0[1:3,1]))
    }

    nonNumericColumn=rep(FALSE, dim(dat0)[[2]]-1)
    for (i in 2:dim(dat0)[[2]] ){ nonNumericColumn[i-1]=! is.numeric(dat0[,i]) }
    if ( sum(nonNumericColumn) >0 ) {
        warning(paste("MAJOR WARNING: Possible input error. The following samples contain nonnumeric beta values: ", colnames(dat0)[-1][ nonNumericColumn], "\n Hint: Maybe you use the wrong symbols for missing data. Make sure to code missing values as NA in the Excel file. To proceed, I will force the entries into numeric values but make sure this makes sense."))
    }
    XchromosomalCpGs=as.character(dnamage.probeAnnotation27k$Name[dnamage.probeAnnotation27k$Chr=="X"])
    selectXchromosome=is.element(dat0[,1], XchromosomalCpGs )
    selectXchromosome[is.na(selectXchromosome)]=FALSE
    meanXchromosome=rep(NA, dim(dat0)[[2]]-1)
    if ( sum(selectXchromosome) >=500 ) {
        meanXchromosome= as.numeric(apply( as.matrix(dat0[selectXchromosome,-1]),2,mean,na.rm=TRUE))
    }
    if ( sum(is.na(meanXchromosome)) >0 ) {
        warning("Comment: There are lots of missing values for X chromosomal probes for some of the samples. This is not a problem when it comes to estimating age but I cannot predict the gender of these samples.")
        match1=match(dnamage.probeAnnotation21kdatMethUsed$Name , dat0[,1])
        if ( sum( is.na(match1))>0 ) {
            missingProbes= dnamage.probeAnnotation21kdatMethUsed$Name[!is.element( dnamage.probeAnnotation21kdatMethUsed$Name , dat0[,1])]
            ##MATT used to stop on this error, changed to a warning
            warning(paste("Input error: You forgot to include the following ", length(missingProbes), " CpG probes (or probe names):\n "))#MATT##, paste( missingProbes, sep="",collapse=", ")))
        }
    }
    
    match1=match(dnamage.probeAnnotation21kdatMethUsed$Name , dat0[,1])
    if ( sum( is.na(match1))>0 ) ##MATT used to stop on this error, changed to a warning
        warning(paste(sum( is.na(match1)), "CpG probes cannot be matched"))
    dat1= dat0[match1,]
    asnumeric1=function(x) {as.numeric(as.character(x))}
    dat1[,-1]=apply(as.matrix(dat1[,-1]),2,asnumeric1)
    
    set.seed(1)
    ## Do you want to normalize the data (recommended)?
    ##MATT##normalizeData=TRUE

    ##MATT added this so that all age CpGs are at least included (check for all CpGs above is a warning now)
    selectCpGsClock=is.element(dat0[,1], as.character(dnamage.datClock$CpGmarker[-1]))
    if ( sum( selectCpGsClock) < dim(dnamage.datClock)[[1]]-1 ) {stop("The CpGs listed in column 1 of the input data did not contain the CpGs needed for calculating DNAm age. Make sure to input cg numbers such as cg00075967.")}
    if ( sum( selectCpGsClock) > dim(dnamage.datClock)[[1]]-1 ) {stop("ERROR: The CpGs listed in column 1 of the input data contain duplicate CpGs. Each row should report only one unique CpG marker (cg number).")}


    ##MATT##source("dnamage/StepwiseAnalysis.txt")
    dnamage.stepwise(dat1, meanXchromosome,
                     normalizeData=normalizeData)
}

dnamage.stepwise <- function(dat1, meanXchromosome,
                             normalizeData) {

    nSamples=dim(dat1)[[2]]-1
    nProbes= dim(dat1)[[1]]


    trafo <- function(x,adult.age=20) { 
        x=(x+1)/(1+adult.age);
        y=ifelse(x<=1, log( x),x-1);
        y
    }
    anti.trafo <- function(x,adult.age=20) {
        ifelse(x<0,
               (1+adult.age)*exp(x)-1,
               (1+adult.age)*x+adult.age)
    }

    ## Steve Horvath: Estimating DNAm age.
    ## This file assumes a data frame exists called dat1 whose rows correspond to CpGs
    ## and whose first column reports the CpG identifier
    ## and whose remaining columns corresponds to samples (e.g. Illumina arrays).
    
    fastImputation=FALSE
    
    ##STEP 1: DEFINE QUALITY METRICS
    
    meanMethBySample =as.numeric(apply(as.matrix(dat1[,-1]),2,mean,na.rm=TRUE))
    minMethBySample   =as.numeric(apply(as.matrix(dat1[,-1]),2,min,na.rm=TRUE))
    maxMethBySample  =as.numeric(apply(as.matrix(dat1[,-1]),2,max,na.rm=TRUE))
    
    datMethUsed= t(dat1[,-1])
    colnames(datMethUsed)=as.character(dat1[,1])
    
    noMissingPerSample=apply(as.matrix(is.na(datMethUsed)),1,sum)
    print(table(noMissingPerSample))
    
    ##STEP 2: Imputing 
    if (! fastImputation & nSamples>1 & max(noMissingPerSample,na.rm=TRUE)<3000 ){
        
        ## run the following code if there is at least one missing
        if ( max(noMissingPerSample,na.rm=TRUE)>0 ){
            dimnames1=dimnames(datMethUsed)
            datMethUsed= data.frame(t(impute.knn(t(datMethUsed))$data))
            dimnames(datMethUsed)=dimnames1
        } # end of if
    } # end of if (! fastImputation )
    
    if ( max(noMissingPerSample,na.rm=TRUE)>=3000 ) fastImputation=TRUE
    
    if ( fastImputation | nSamples==1 ){
        noMissingPerSample=apply(as.matrix(is.na(datMethUsed)),1,sum)
        table(noMissingPerSample)
        if ( max(noMissingPerSample,na.rm=TRUE)>0 & max(noMissingPerSample,na.rm=TRUE) >= 3000 ) {normalizeData=FALSE}

        ## run the following code if there is at least one missing
        if ( max(noMissingPerSample,na.rm=TRUE)>0 & max(noMissingPerSample,na.rm=TRUE) < 3000 ){
            dimnames1=dimnames(datMethUsed)
            for (i in which(noMissingPerSample>0) ){
                selectMissing1=is.na(datMethUsed[i,])
                datMethUsed[i,selectMissing1] = as.numeric(dnamage.probeAnnotation21kdatMethUsed$goldstandard2[selectMissing1])
            } # end of for loop
            dimnames(datMethUsed)=dimnames1
        } # end of if
    } # end of if (! fastImputation )


    ## STEP 3: Data normalization (each sample requires about 8 seconds). It would be straightforward to parallelize this operation.
    
    if (normalizeData ){
        datMethUsedNormalized=BMIQcalibration(datM=datMethUsed,goldstandard.beta= dnamage.probeAnnotation21kdatMethUsed$goldstandard2,plots=FALSE)
    }
    if (!normalizeData ){ datMethUsedNormalized=datMethUsed }
    rm(datMethUsed); gc()
    
    ##STEP 4: Predict age and create a data frame for the output (referred to as datout)
    selectCpGsClock=is.element(dimnames(datMethUsedNormalized)[[2]], as.character(dnamage.datClock$CpGmarker[-1]))
    if ( sum( selectCpGsClock) < dim(dnamage.datClock)[[1]]-1 ) {stop("The CpGs listed in column 1 of the input data did not contain the CpGs needed for calculating DNAm age. Make sure to input cg numbers such as cg00075967.")}
    if ( sum( selectCpGsClock) > dim(dnamage.datClock)[[1]]-1 ) {stop("ERROR: The CpGs listed in column 1 of the input data contain duplicate CpGs. Each row should report only one unique CpG marker (cg number).")}
    if (nSamples>1 ) {
        datMethClock0=data.frame(datMethUsedNormalized[,selectCpGsClock])
        datMethClock= data.frame(datMethClock0[ as.character(dnamage.datClock$CpGmarker[-1])])
        dim(datMethClock)
        predictedAge=as.numeric(anti.trafo(dnamage.datClock$CoefficientTraining[1]+as.matrix(datMethClock)%*% as.numeric(dnamage.datClock$CoefficientTraining[-1])))
    } # end of if
    
    
    if (nSamples==1 ) {
        datMethUsedNormalized2=data.frame(rbind(datMethUsedNormalized,datMethUsedNormalized))
        datMethClock0=data.frame(datMethUsedNormalized2[,selectCpGsClock])
        datMethClock= data.frame(datMethClock0[ as.character(dnamage.datClock$CpGmarker[-1])])
        dim(datMethClock)
        predictedAge=as.numeric(anti.trafo(dnamage.datClock$CoefficientTraining[1]+as.matrix(datMethClock)%*% as.numeric(dnamage.datClock$CoefficientTraining[-1])))
        predictedAge=predictedAge[1]
    } # end of if
    
    
    
    ## Let's add comments to the age prediction
    Comment=ifelse ( predictedAge <0, "Negative DNAm age.", ifelse ( predictedAge >100, "Old DNAm age.", rep("",length(predictedAge))))
    
    Comment[is.na(predictedAge)]="Age prediction was not possible. "
    
    
    if ( sum( selectCpGsClock) < dim(dnamage.datClock)[[1]]-1 ) {
        Comment=rep("ERROR: The CpGs listed in column 1 of the input data did not contain the CpGs needed for calculating DNAm age. Make sure to input cg numbers such as cg00075967.",length(predictedAge) )}
    
    
    if ( sum( selectCpGsClock) > dim(dnamage.datClock)[[1]]-1 ) {
        Comment=rep("ERROR: The CpGs listed in column 1 of the input data contain duplicate CpGs. Each row should report only one unique CpG marker (cg number).",length(predictedAge) )}
    
    
    restSamples=-minMethBySample>0.05 | maxMethBySample>1.05;
    restSamples[is.na(restSamples)]=FALSE
    lab1="MAJOR WARNING: Probably you did not input beta values since either minMethBySample<-0.05 or maxMethBySample>1.05.";Comment[restSamples]= paste(Comment[restSamples],lab1)
    
    restSamples= noMissingPerSample >0 & noMissingPerSample <=100;lab1="WARNING: Some beta values were missing, see noMissingPerSample."; Comment[restSamples]= paste(Comment[restSamples],lab1)
    restSamples= noMissingPerSample >3000;lab1="MAJOR WARNING: More than 3k missing values!!"; Comment[restSamples]= paste(Comment[restSamples],lab1)
    
    restSamples= noMissingPerSample >100 & noMissingPerSample <=3000 ;lab1="MAJOR WARNING: noMissingPerSample>100"
    Comment[restSamples]= paste(Comment[restSamples],lab1)
    restSamples=meanMethBySample>.35;
    restSamples[is.na(restSamples)]=FALSE
    lab1="Warning: meanMethBySample is >0.35";Comment[restSamples]= paste(Comment[restSamples],lab1)
    restSamples=meanMethBySample<.25;
    restSamples[is.na(restSamples)]=FALSE; lab1="Warning: meanMethBySample is <0.25"
    Comment[restSamples]= paste(Comment[restSamples],lab1)
    datout=data.frame(SampleID=colnames(dat1)[-1], DNAmAge=predictedAge, Comment, noMissingPerSample,meanMethBySample, minMethBySample, maxMethBySample)
    
    if ( !is.null( meanXchromosome) ){  
        
        if ( length( meanXchromosome)==dim(datout)[[1]] ){
            predictedGender=ifelse(meanXchromosome>.4,"female",
                ifelse(meanXchromosome<.38,"male","Unsure"))
            datout=data.frame(datout,predictedGender=predictedGender,meanXchromosome=meanXchromosome)
            
        } # end of if 
        
    } # end of if

    datout
}
