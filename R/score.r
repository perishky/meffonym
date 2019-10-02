#' Calculate DNA methylation scores
#'
#' @param x DNA methylation matrix (rows=CpG sites, columns=samples).
#' @param model Model name from available list (\code{\link{meffonym.models}()}).
#' @return Scores derived for each sample in \code{x}.
#' Missing values in \code{x} are replaced with the mean value of the row.
#'
#' @export
meffonym.score <- function(x, model) {
    stopifnot(is.matrix(x))

    x <- impute.matrix(x,1) ## replace missing values with mean values
    
    ret <- meffonym.get.model(model)
    ret$name <- model

    sites <- intersect(rownames(x), names(ret$coefficients))

    score <- rep(NA, ncol(x))
    if (length(sites) == 1)
        score <- ret$intercept + ret$coefficients[sites] * x[sites,]
    if (length(sites) > 1)
        score <- ret$intercept + as.vector(rbind(ret$coefficients[sites]) %*% x[sites,,drop=F])

    if (model == "age.horvath") {
        anti.trafo <- function(x,adult.age=20) {
            ifelse(x<0,
                   (1+adult.age)*exp(x)-1,
                   (1+adult.age)*x+adult.age)
        }
        score <- anti.trafo(score)
    }

    ret$sites <- sites
    ret$score <- score
    ret
}
