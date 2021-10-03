read.gse.matrix.file <- function(filename) {
    dat <- readLines(filename)
    nseries <- sum(grepl("^!Series_", dat))
    nsamples <- sum(grepl("^!Sample_", dat))
    ndata <- length(dat) - match("!series_matrix_table_begin", dat) - 2
    con <- file(filename, "r")
    header <- read.table(con, sep="\t", header=F, nrows=nseries, stringsAsFactors=F)
    samples <- read.table(con, sep="\t", header=F, nrows=nsamples, stringsAsFactors=F)
    samples <- t(samples)
    colnames(samples) <- samples[1,]
    colnames(samples) <- sub("!Sample_", "", colnames(samples))
    samples <- data.frame(samples[-1,], stringsAsFactors=F)

    rm(dat)
    gc()

    readLines(con,1)
    dnam <- read.table(con, sep="\t", header=TRUE, quote="\"", dec=".", fill=TRUE,
                       na.strings = c("NA", "null", "NULL", "Null"), comment.char = "")
    close(con)
    if (ndata == 0)
        dnam <- dnam[-(1:nrow(dnam)),]
    
    if (nrow(dnam) > 0) 
        rownames(dnam) <- dnam[,1]
    dnam <- as.matrix(dnam[,-1])
    
    rownames(samples) <- colnames(dnam)
    colnames(dnam) <- samples$geo_accession

    list(dnam=dnam, samples=samples)
}

