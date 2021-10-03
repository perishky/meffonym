# to impute missing values with row means
# x <- impute.mean(x,1)
# to impute missing values with column means
# x <- impute.mean(x,2)
impute.mean <- function(
    x,
    margin=1,
    fun=function(x) mean(x, na.rm=T),
    na.rm=F) {
    
    if (margin == 2) x <- t(x)
    
    idx <- which(is.na(x) | !is.finite(x), arr.ind=T)
    if (length(idx) > 0) {
        na.idx <- unique(idx[,"row"])
        v <- apply(x[na.idx,,drop=F],1,fun) ## v = summary for each row
        v[which(is.na(v))] <- fun(v)      ## if v[i] is NA, v[i] = fun(v)
        x[idx] <- v[match(idx[,"row"],na.idx)] ##
        stopifnot(all(!is.na(x)))
    }

    if (margin == 2) x <- t(x)

    if (na.rm==T) {
        if (margin == 1) { ## drop rows with all missing values
            idx <- which(!is.na(x[,1]))
            x <- x[idx,,drop=F]
        }
        else { ## drop columns with all missing values
            idx <- which(!is.na(x[1,]))
            x <- x[,idx,drop=F]
        }
    }
    
    x
}

## impute matrix values using knn
impute.matrix <- function(x) {
    dn <- dimnames(x)
    x <- impute.knn(x)$data
    dimnames(x) <- dn
    x
}
