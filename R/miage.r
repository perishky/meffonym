#' Estimate mitotic age from DNAm using the MiAge algorithm
#' 
#' Define the following:
#' * x[i,j] be the methylation level of CpG site i in sample j
#' * E(X[n,i]) be the expected methylation level of CpG site i after n generations,
#' * a[i] be the rate of de novo methylation of CpG site i
#' * b[i] be the fidelity of methylation maintenance of CpG site i
#' * c[i] = a[i]/(1-b[i]) 
#' * d[i] = E(X[1,i]) - c[i]
#' 
#' This function estimates the number of generations n[j] for each sample j
#' from the DNA methylation levels x[,j] of 268 CpG sites.
#' The estimate is obtained by selecting the value n[j] that minimizes the sum of
#' (E(X[n,i])-x[i,j])^2 across CpG sites i=1..268
#' where E(X[n,i]) is estimated as c[i] + b[i]^(n[j]-1)*d[i]
#' for given b[i], c[i] and d[i].
#'
#' @param x DNA methylation matrix (rows=CpG sites, columns=samples).
#' Missing values in \code{x} will be replaced with the mean value of the row.
#' @param b As defined above (Default: published values).
#' @param c As defined above (Default: published values).
#' @param d As defined above (Default: published values).
#' @param minage Lowest age to consider (Default: 10).
#' @param maxage Highest age to consider (Default: 10000).
#' @param inits Initial ages to consider at each iteration.
#' (Default: 500 and lowerage + (1:4)*(maxage-minage)/5).
#' @return Estimated mitotic age for each sample.
#' 
miage <- function(
    x,
    b=NULL,
    c=NULL,
    d=NULL,
    minage=10,
    maxage=10000,
    inits=c(500, minage + (1:4)*(maxage-minage)/5)) {
    
    objective.f <- function(nj,b,c,d,xj) {
        return(sum((c+b^(nj-1)*d-xj)^2,na.rm=T))
    }
    
    derivative.f <- function(nj,b,c,d,xj) {
        return(2*sum((c+b^(nj-1)*d-xj)*b^(nj-1)*log(b)*d,na.rm=T))
    }

    optimize.f <- function(init, xj) {
        tryCatch({
            ret <- stats::optim(
                par=init, 
                fn=objective.f,
                gr=derivative.f,
                b=b,
                c=c,
                d=d,
                xj=xj,
                method="L-BFGS-B",
                lower=minage,
                upper=maxage,
                control=list(factr=1))
            list(ss=ret$value, n=ret$par)
        }, error=function(e) {
            list(ss=NA, n=NA)
        })
    }
    
    sapply(1:ncol(x), function(j) {
        ret <- lapply(inits, optimize.f, xj=x[,j])
        idx <- which.min(sapply(ret, function(ret) ret$ss))
        if (!is.na(idx))
            ret[[idx]]$n
        else
            NA
    })
}

