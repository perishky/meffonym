#' Calculate DNA methylation scores
#'
#' @param x A matrix of all score input variables, such as DNA methylation values 
#' (rows=Input variables like CpG sites, columns=samples).
#' Missing values in \code{x} will be replaced with the mean value of the row.
#' @param model Model name from available list
#' (\code{\link{meffonym.models}()}).
#' @param calibrate Prior to calculating scores, if TRUE then
#' calibrate DNA methylation, otherwise use DNA methylation unmodified
#' (Default: FALSE).
#' If `model=="dunedinpace"`, then calibration involves quantile
#' normalization to the DunedinPACE gold standard;
#' otherwise, calibration involves applying the BMIQ algorithm
#' to the Horvath gold standard as in Horvath's DNAm age calculator.
#' @param scale Prior to calculating scores, standardize methylation levels,
#' i.e. subtract mean and divide by standard deviation for each CpG site.
#' Default is \code{TRUE} for 'zhang' and 'zhang.blup',  otherwise \code{FALSE}.
#' @param transform After calculating scores, transform scores
#' using this function. Default is \code{\link{anti.trafo}} for
#' the "horvath", "skin" and "pedbe" models, otherwise \code{NULL}.
#' @param adjust Data frame with variables for adjusting the final score (Default: \code{NULL}).
#' @return List containing the list of score variables (e.g. CpG sites) actually used
#' (\code{sites}) in case
#' some score variables included in the model are missing from \code{x},
#' the resulting score (\code{score}) for each sample in the dataset.
#' If \code{adjust} is not \code{NULL}, then the score prior to adjustment
#' will also be included in the return list (\code{raw}).
#'
#' @export
meffonym.score <- function(
    x,
    model,
    calibrate=F,
    scale=model %in% c("zhang","zhang.blup"),
    transform=NULL,
    adjust=NULL) {

    ## methylation data should be a matrix
    stopifnot(is.matrix(x))

    if (is.null(transform)
        &&  model %in% c("horvath","skin","pedbe"))
        transform <- anti.trafo
  
    if (!is.null(adjust)) {
        ## ensure number of samples in
        ## adjustment variables matches 
        ## same as in methylation data
        stopifnot(is.matrix(adjust) || is.data.frame(adjust))
        stopifnot(nrow(adjust) == ncol(x))
    }
    
    ret <- meffonym.get.model(model)
    ret$name <- model

    ## retrieve CpG sites available for model
    sites <- intersect(rownames(x), ret$vars)
    if (length(sites) == 0) 
        stop("x does not contain data for score input variables
         in the model")

    if (calibrate) {
        x.orig <- x
        if (model == "dunedinpace") {
            std <- meffonym.dunedinpace.standard()
            x <- meffonym.quantile.normalization(x,std)
        } else {
            ## calibrate the methylation data by Horvath standard
            std <- meffonym.horvath.standard()
            x <- meffonym.bmiq.calibration(x, std)
        }
        skipped.vars <- setdiff(rownames(x.orig), rownames(x))
        if (length(skipped.vars) > 0) 
            x <- rbind(x, x.orig[skipped.vars,,drop=F])
    }

    ## restrict methylation data to CpG sites in model
    x <- x[sites,,drop=F]

    if (any(is.na(x))) {
        ## impute means for missing methylation values
        x <- impute.mean(x,1,na.rm=T)
        sites <- rownames(x)
        if (length(sites) < 1)
            stop("x does not contain data for score input variables
             in the model")
    }

    num.missing <- length(setdiff(names(ret$coefs), sites))
    if (num.missing)
        ## let user know how many sites being used for the model
        warning(paste("Dataset missing",
                      num.missing,
                      "score input variables for model", model))    

    if (scale) {
        cols <- colnames(x)
        x <- t(apply(x,1,scale))
        colnames(x) <- cols
    }
    
    if (model == "miage") {
        score <- miage(
            x=x,
            b=ret$coefs$b,
            c=ret$coefs$c,
            d=ret$coefs$d)
    }
    else if (model == "epitoc") {    
        score <- colMeans(x)
    }
    else if (model == "epitoc2") {
        beta0 <- ret$coefs[sites,"beta0"]
        delta <- ret$coefs[sites,"delta"]
        coefs <- diag(1/(delta*(1-beta0)))
        score <- 2*colMeans(coefs %*% (x-beta0))
    }
    else {
        ## default: linear model
        score <- ret$intercept
        if (length(sites) == 1)
            score <- score + ret$coefs[sites] * as.vector(x)
        else
            score <- score + as.vector(rbind(ret$coefs[sites]) %*% x)
    }

    if (!is.null(transform))
        ## transform score using input transform() function
        score <- transform(score)
    
    ret$sites <- sites
    ret$score <- score

    if (!is.null(adjust)) {
        ## adjust score with supplied adjustment variables
        if (is.matrix(adjust))
            adjust <- as.data.frame(adjust)
        fit <- lm(score ~ ., adjust)
        ret$raw <- score
        ret$score <- residuals(fit)
    }

    ret
}

