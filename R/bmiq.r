#' BMIQ normalization to a standard
#'
#' Methylation levels are normalized so that they are comparable to a given dataset. It is based on
#' BMIQ:
#'
#' A beta-mixture quantile normalization method for correcting probe design bias in Illumina Infinium 450 k DNA methylation data.
#' Teschendorff AE, Marabita F, Lechner M, Bartlett T, Tegner J, Gomez-Cabrero D, Beck S.
#' Bioinformatics. 2013 Jan 15;29(2):189-96. doi: 10.1093/bioinformatics/bts680. Epub 2012 Nov 21.
#' PMID: 23175756
#' 
#' The code was obtained from \url{http://labs.genetics.ucla.edu/horvath/dnamage/NORMALIZATION.R}.
#'
#' @param x DNA methylation matrix (rows=CpG sites, columns=samples).
#' @param standard A named vector of methylation levels corresponding to the row names of \code{x}.
#' @param ... See \url{http://code.google.com/p/bmiq/} for more details.
#' @return Matrix identical to \code{x} but with methylation levels normalized
#' by BMIQ to the \code{standard} and rows corresponding to \code{standard}.
#' Missing values in \code{x} are imputed by k-nearest neighbor.
#' 
#' @export 
meffonym.bmiq.calibration <- function(x,standard,...) {
    missing.sites <- setdiff(names(standard), rownames(x))
    if (length(missing.sites) > 0) {
        if (length(missing.sites) == length(standard))
            stop("All CpG sites for DNAm clock are missing.")
        x <- rbind(x,
                   t(sapply(standard[missing.sites], rep, ncol(x))))
    }
    x <- x[names(standard),]
    
    if (any(is.na(x))) {
        dn <- dimnames(x)
        x <- impute.knn(x)$data
        dimnames(x) <- dn
    }
    t(BMIQcalibration(t(x),
                      goldstandard.beta=standard,
                      ...))
}

