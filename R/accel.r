#' Calculate age acceleration 
#'
#' @param x,model,calibrate,transform,adjust See \code{\link{meffonym.score}}.
#' @param age Chronological age in years.
#' @return Age-adjusted model scores.
#' 
#' @export
meffonym.accel <- function(
    x,
    model,
    age,
    calibrate=F,
    transform=ifelse(model %in% c("horvath","pedbe"), anti.trafo, NULL),
    adjust=NULL) {
    
    stopifnot(length(age) == ncol(x))
    
    if (is.null(adjust))
        adjust <- data.frame(age=age)
    else {
        stopifnot(is.matrix(adjust) || is.data.frame(adjust))
        stopifnot(nrow(adjust) == length(age))
        adjust <- cbind(adjust, age=age)
    }
    meffonym.score(x,model,adjust)$score
}

