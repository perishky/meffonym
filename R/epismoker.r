#' Epismoker linear predictor for a specified smoking status
#' 
#' @param x DNA methylation matrix (rows=CpG sites, columns=samples).
#' @param group Either 'current', 'former' or 'never'
#' @param female Vector of booleans indicating whether each individual is female.
#' @return Epismoker smoking group linear predictor derived for each sample in \code{x}.
#' 
#' @export
meffonym.epismoker.score <- function(x,group,female,calibrate=F) {
    stopifnot(group %in% c("current","former","never"))
    stopifnot(is.vector(female) && length(female) == ncol(x))
    meffonym.score(
        rbind(x,sexM=sign(!female)),
        paste0("epismoker-",group),
        calibrate=calibrate)$score
}

#' Epismoker smoking status probabilities
#'
#' @param x DNA methylation matrix (rows=CpG sites, columns=samples).
#' @param female Vector of booleans indicating whether each individual is female.
#' @return Matrix of epismoker smoking status probabilities (rows=samples, columns=smoking status) 
#' 
#' @export 
meffonym.epismoker <- function(x,female,calibrate=F) {
    stopifnot(is.vector(female) && length(female) == ncol(x))
    scores = cbind(
        current=meffonym.epismoker.score(x,"current",female),
        former=meffonym.epismoker.score(x,"former",female),
        never=meffonym.epismoker.score(x,"never",female))
    scores = exp(scores)
    scores/rowSums(scores,na.rm=T)
}
