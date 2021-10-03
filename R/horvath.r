#' Horvath 21k methylation standard
#'
#' Obtained from
#' \url{http://labs.genetics.ucla.edu/horvath/dnamage/probeAnnotation21kdatMethUsed.csv}.
#' 
#' @export
meffonym.horvath.standard <- function() {
    horvath.standard
}

load.horvath.standard <- function(pkgname) {
    path <- system.file("horvath", package=pkgname)
    ## http://labs.genetics.ucla.edu/horvath/dnamage/probeAnnotation21kdatMethUsed.csv
    filename <- file.path(path, "probeAnnotation21kdatMethUsed.csv")
    if (!file.exists(filename))
        return(NULL)
    
    standard <- read.csv(filename, stringsAsFactors=F)
    sites <- standard$Name
    standard <- standard$goldstandard2
    names(standard) <- sites
    standard
}

#' Horvath DNA methylation age estimate
#' 
#' @param x DNA methylation matrix (rows=CpG sites, columns=samples).
#' @return Scores derived for each sample in \code{x}.
#'
#' @export
meffonym.horvath.age <- function(x) {
    meffonym.score(x, "horvath", calibrate=T)$score
}


