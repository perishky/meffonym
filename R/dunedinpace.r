#' DunedinPACE methylation standard
#'
#' Obtained from
#' \url{https://github.com/danbelsky/DunedinPACE}.
#' 
#' @export
meffonym.dunedinpace.standard <- function() {
    dunedinpace.standard
}

load.dunedinpace.standard <- function(pkgname) {
    path <- system.file("dunedinpace", package=pkgname)
    filename <- file.path(path, "standard.csv")
    if (!file.exists(filename))
        return(NULL)
    
    standard <- read.csv(filename, stringsAsFactors=F)
    sites <- standard$name
    standard <- standard$mean
    names(standard) <- sites
    standard
}

#' DunedinPACE pace of aging estimate
#' 
#' @param x DNA methylation matrix (rows=CpG sites, columns=samples).
#' @return Scores derived for each sample in \code{x}.
#'
#' @export
meffonym.dunedinpace.estimate <- function(x) {
    meffonym.score(x, "dunedinpace", calibrate=T)$score
}
