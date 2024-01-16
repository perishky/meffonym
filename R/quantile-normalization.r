#' Quantile normalization to a standard
#'
#' Methylation levels are quantile normalized
#' so that they are comparable to a given dataset.
#'
#' @param x DNA methylation matrix (rows=CpG sites, columns=samples).
#' @param standard A named vector of methylation levels corresponding to the row names of \code{x}.
#' @return Matrix identical to \code{x} but with methylation levels
#' quantile normalized to the \code{standard}
#' and rows corresponding to \code{standard}.
#' Missing values in \code{x} are imputed by k-nearest neighbor.
#' 
#' @export 
meffonym.quantile.normalization <- function(x,standard,...) {
    missing.sites <- setdiff(names(standard), rownames(x))
    if (length(missing.sites) > 0) {
        if (length(missing.sites) == length(standard))
            stop("All CpG sites for DNAm clock are missing.")
        warning("Some missing sites:", length(missing.sites))
        x <- rbind(x,
                   t(sapply(standard[missing.sites], rep, ncol(x))))
    }
    x <- x[names(standard),]
    
    if (any(is.na(x))) 
        x <- impute.matrix(x)

    x.norm <- preprocessCore::normalize.quantiles.use.target(x, target=standard)
    dimnames(x.norm) <- dimnames(x)
    x.norm
}

