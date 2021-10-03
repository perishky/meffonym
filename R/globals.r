.onLoad <- function(libname, pkgname) {
    env <- parent.env(environment())

    cat("loading models ...\n")
    assign("models.global", new.env(), envir=env)
    load.models(pkgname)

    cat("load horvath standard ...\n")
    standard <- load.horvath.standard(pkgname)
    assign("horvath.standard", standard, envir=env)
}

