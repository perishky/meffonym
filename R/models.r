#' List names of available models
#'
#' @export
meffonym.models <- function() ls(models.global)

#' Add or update set of available models
#'
#' @param name Name of the model.
#' @param variables Names of CpG sites included in the model.
#' @param coefficients Coefficients of the CpG sites in the model.
#' @param intercept Model intercept (Default: 0).
#' @param description Text description of the model, e.g. publication.
#'
#' @export
meffonym.add.model <- function(name, variables, coefficients, description, intercept=0) {
    model <- create.model(variables, coefficients, intercept=intercept, description=description)
    assign(name, model, envir=models.global)
    invisible(TRUE)
}

#' Retrieve model details
#'
#' @param name Name of the model.
#' @return Model details formatted as a named list.
#' 
#' @export
meffonym.get.model <- function(name) {
    stopifnot(name %in% meffonym.models())
    get(name, envir=models.global)
}

create.model <- function(vars, coefs, intercept, description) {
    model <- list()
    if (is.vector(coefs)) {
        stopifnot(length(vars) == length(coefs))
        names(coefs) <- vars
        model$intercept=intercept
        vars <- c('intercept', vars)
    } else {
        stopifnot(is.matrix(coefs) || is.data.frame(coefs))
        rownames(coefs) <- vars
    }
    model$vars <- vars
    model$coefs <- coefs
    model$description <- description
    model
}

load.models <- function(pkgname) {
    filename <- system.file("models.csv", package=pkgname)
    if (!file.exists(filename))
        return
    
    path <- dirname(filename)
    models <- read.csv(filename, stringsAsFactors=F)

    for (i in 1:nrow(models)) {
        filename <- file.path(path, models$filename[i])
        cat("loading", filename, " ...\n")
        coefs <- read.csv(filename, stringsAsFactors=F)
        if (ncol(coefs) == 2) {
            cpgs <- coefs$cpg
            coefs <- coefs$coef
            if (!any(grepl("intercept", cpgs, ignore.case=T))) {
                cpgs <- c("intercept", cpgs)
                coefs <- c(0,coefs)
            }
            intercept <- coefs[1]
            cpgs <- cpgs[-1]
            coefs <- coefs[-1]
        }
        else {
            cpgs <- coefs$cpg
            coefs$cpg <- NULL
            intercept <- 0
        }
        meffonym.add.model(
            models$name[i],
            cpgs,
            coefs,
            intercept=intercept,
            description=paste(paste(colnames(models), models[i,], sep=":"), collapse=","))
    }
}

