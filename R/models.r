#' List names of available models
#'
#' @export
meffonym.models <- function() ls(models.global)

#' Add or update set of available models
#'
#' @param name Name of the model.
#' @param variables Names of CpG sites included in the model.
#' @param coefficients Coefficients of the CpG sites in the model.
#' @param description Text description of the model, e.g. publication.
#'
#' @export
meffonym.add.model <- function(name, variables, coefficients, description) {
    model <- create.model(variables, coefficients, description)
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

create.model <- function(variables, coefficients, description) {
    stopifnot(length(variables) == length(coefficients))
    names(coefficients) <- variables
    list(intercept=coefficients[1], coefficients=coefficients[-1], description=description)
}

load.models <- function(pkgname) {
    ## https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1068-z/MediaObjects/13059_2016_1068_MOESM3_ESM.csv
    path <- system.file("knight", package=pkgname)                               
    model <- read.csv(file.path(path, "model.csv"), stringsAsFactors=F)
    meffonym.add.model("ga.knight", model$CpGmarker, model$CoefficientTraining,
                         readLines(file.path(path, "readme.txt")))
    
    ## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3780611/bin/NIHMS418935-supplement-02.xlsx
    path <- system.file("hannum", package=pkgname)
    model <- read.csv(file.path(path, "hannum-model.csv"), stringsAsFactors=F)
    meffonym.add.model("age.hannum", c("intercept", model$Marker), c(0, model$Coefficient),
                         readLines(file.path(path, "readme.txt")))
    
    ## see inst/bohlin/retrieve.r
    path <- system.file("bohlin", package=pkgname)
    model <- read.csv(file.path(path, "model-1se.csv"), stringsAsFactors=F)
    meffonym.add.model("ga.bohlin.1se", model$cpg, model$coefficient,
                         readLines(file.path(path, "readme.txt")))
    
    path <- system.file("bohlin", package=pkgname)
    model <- read.csv(file.path(path, "model-min.csv"), stringsAsFactors=F)
    meffonym.add.model("ga.bohlin.min", model$cpg, model$coefficient, readLines(file.path(path, "readme.txt")))

    ## http://labs.genetics.ucla.edu/horvath/dnamage/AdditionalFile3.csv
    path <- system.file("horvath", package=pkgname)
    model <- read.csv(file.path(path, "AdditionalFile3.csv"), stringsAsFactors=F)
    meffonym.add.model("age.horvath", model$CpGmarker, model$CoefficientTraining,
                         readLines(file.path(path, "readme.txt")))
}

