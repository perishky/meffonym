#### epiTOC2.R

#### Author: Andrew E Teschendorff (a.teschendorff@ucl.ac.uk)
#### Date: 8th Apr.2019
#### Copyright 2019 Andrew Teschendorff
#### Copyright permission: epiTOC2 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version-3 as published by the Free Software  Foundation. epiTOC2 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details (http://ww.gnu.org/licenses/).


#### DESCRIPTION
#### An R-function to estimate the cumulative number of stem-cell divisions in a sample using the epiTOC2 model. Only required input argument is a DNA methylation data matrix, normalized with BMIQ or another type-2 probe correction method. If the ages of the samples are provided and assuming all samples are from the same tissue-type, function also returns the estimated intrinsic rate of stem-cell division for each sample as well as a median estimate for the tissue. The function returns mitotic age estimates for two different epiTOC2 models, which differ in the ground-state methylation values for the epiTOC2 CpGs. In the full model, we use the estimated methylation ground state values, whereas in the simplified model we assume that these are all zero. Reason for returning the values estimated using the simplified model is because it could happen that some methylation beta-values in the data matrix are lower than the ground-state values, which in principle is not allowed. If this is the case, then more reliable estimates are provided by using the simplified model. The functions alerts the user to this, if it detects lots of beta-values less than ground-state values. Finally, the function also returns the pcgtAge-score and HypoClock score for each column (sample) in the input data matrix using the previous epiTOC model and solo-WCGWs, respectively.

#### REQUIRED OBJECTS:
#### dataETOC2.Rd: this object file will be loaded and must reside in the working directory. It is a list with the 3 elements. The first element is a matrix with rows labeling the 163 epiTOC2 CpGs, and columns labeling the estimated de-novo and ground-state methylation parameters. The 2nd element is a vector of CpG identifiers for the 385 epiTOC CpGs. The 3rd element is a vector of CpG identifiers for the 678 solo-WCGWs which are constitutively methylated in fetal tissue. 

#### INPUT:
#### data.m: DNAm data beta-valued matrix with rownames labeling CpGs and columns labeling samples. All samples should be from the same tissue-type. 
#### ages.v: Optional argument representing the chronological ages or surrogates thereof of the samples. Vector must be of same length as the number of samples.



#### OUPTUT:
#### A list containing the following entries
#### tnsc: the estimated cumulative number of stem-cell divisions per stem-cell per year and per sample using the full epiTOC2 model.
#### tnsc2: the estimated cumulative number of stem-cell divisions per stem-cell per year and per sample using an approximation of epiTOC2 which assumes all epiTOC2 CpGs have beta-values exactly 0 in the fetal stage.
#### irS: this is returned only if the ages are provided, and gives the estimated average lifetime intrinsic rate of stem-cell division per sample, as derived from epiTOC2
#### irS2: as irS, but for the approximation.
#### irT: the median estimate over all irS values, yielding a median estimate for the intrinsic rate of stem-cell division for the tissue.
#### irT2: as irT, but for the approximation.
#### pcgtAge: this is the mitotic-score obtained using our previous epiTOC model.
#### hypoSC: the HypoClock score over the 678 solo-WCGWs

epiTOC2 <- function(data.m,ages.v=NULL){
    load("dataETOC2.Rd"); ## this loads the CpG information
    cpgETOC.v <- dataETOC2.l[[2]];
    estETOC2.m <- dataETOC2.l[[1]];
    soloCpG.v <- dataETOC2.l[[3]];
    ### do epiTOC
    common.v <- intersect(rownames(data.m),cpgETOC.v);
    print(paste("Number of represented epiTOC CpGs (max=385)=",length(common.v),sep=""));
    map.idx <- match(common.v,rownames(data.m));
    pcgtAge.v <- colMeans(data.m[map.idx,],na.rm=TRUE);
    ### do epiTOC2
    map.idx <- match(rownames(estETOC2.m),rownames(data.m));
    rep.idx <- which(is.na(map.idx)==FALSE);
    print(paste("Number of represented epiTOC2 CpGs (max=163)=",length(rep.idx),sep=""))
    tmp.m <- data.m[map.idx[rep.idx],];
    TNSC.v <- 2*colMeans(diag(1/(estETOC2.m[rep.idx,1]*(1-estETOC2.m[rep.idx,2]))) %*% (tmp.m - estETOC2.m[rep.idx,2]),na.rm=TRUE);
    TNSC2.v <- 2*colMeans(diag(1/estETOC2.m[rep.idx,1]) %*% tmp.m,na.rm=TRUE);
    ### do HypoClock
    common.v <- intersect(rownames(data.m),soloCpG.v);
    print(paste("Number of represented solo-WCGWs (max=678)=",length(common.v),sep=""));
    map.idx <- match(common.v,rownames(data.m));
    hypoSC.v <- colMeans(data.m[map.idx,],na.rm=TRUE);

    estIR.v <- NULL; estIR2.v <- NULL;
    estIR <- NULL;  estIR2 <- NULL;
    if(!is.null(ages.v)){
      estIR.v <- TNSC.v/ages.v;
      estIR <- median(estIR.v,na.rm=TRUE);
      estIR2.v <- TNSC2.v/ages.v;
      estIR2 <- median(estIR2.v,na.rm=TRUE);
    }

    
    return(list(tnsc=TNSC.v,tnsc2=TNSC2.v,irS=estIR.v,irS2=estIR2.v,irT=estIR,irT2=estIR2,pcgtAge=pcgtAge.v,hypoSC=hypoSC.v));
}

