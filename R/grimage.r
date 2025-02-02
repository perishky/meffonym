#' Grimage
#'
#' Obtained from
#' \url{https://github.com/bio-learn/biolearn}.
#'
#' @export
meffonym.grimage <- function(x,female,age,version=1,calibrate=F) {
    stopifnot(is.vector(female) && length(female) == ncol(x))
    stopifnot(is.vector(age) && length(age) == ncol(x))
    modelname <- ifelse(version==1, "grimage", "grimagev2")
    meffonym.score(
        rbind(x,age=age,female=female),
        modelname,
        calibrate=calibrate)$score
}



