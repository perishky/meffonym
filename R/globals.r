.onLoad <- function(libname, pkgname) {
    env <- parent.env(environment())

    assign("models.global", new.env(), envir=env)
    load.models(pkgname)
    
    standard <- load.horvath.standard(pkgname)
    assign("horvath.standard", standard, envir=env)
}

