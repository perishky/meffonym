#' Horvath 21k methylation standard
#'
#' Obtained from
#' \url{http://labs.genetics.ucla.edu/horvath/dnamage/probeAnnotation21kdatMethUsed.csv}.
#' 
#' @examples
#' x <- ... ## methylation matrix
#' standard <- meffonym.horvath.standard()
#' x.norm <- meffonym.bmiq.calibration(x, standard)
#' ret <- meffonym.score(x, "age.horvath")
#' 
#' @export
meffonym.horvath.standard <- function() {
    horvath.standard
}

load.horvath.standard <- function(pkgname) {
    path <- system.file("horvath", package=pkgname)
    ## http://labs.genetics.ucla.edu/horvath/dnamage/probeAnnotation21kdatMethUsed.csv
    standard <- read.csv(file.path(path, "probeAnnotation21kdatMethUsed.csv"), stringsAsFactors=F)
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
#' @examples
#' ## Equivalent to the following:
#' x <- ... ## methylation matrix
#' standard <- meffonym.horvath.standard()
#' x.norm <- meffonym.bmiq.calibration(x, standard)
#' ret <- meffonym.score(x, "age.horvath")
#' 
#' @export
meffonym.horvath.age <- function(x) {
    standard <- meffonym.horvath.standard()
    sites <- intersect(names(standard), rownames(x))
    if (length(sites) < 2)
        stop("Little or no overlap between methylation data and Horvath methylation standard.")    
    x.norm <- meffonym.bmiq.calibration(x, standard[sites])
    meffonym.score(x.norm, "age.horvath")
}


