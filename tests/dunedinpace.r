## Compare pace of age estimates against the DunedinPACE R package

library(meffonym)
library(DunedinPACE) ## remotes::install_github("danbelsky/DunedinPACE")
library(readr)

sites <- c(
    DunedinPACE::mPACE_Models$gold_standard_probes$DunedinPACE,
    DunedinPACE::mPACE_Models$model_probes$DunedinPACE)
sites <- unique(sites)
meth <- matrix(runif(length(sites)*10), ncol=10)
rownames(meth) <- sites
colnames(meth) <- paste0("s",1:ncol(meth))

noob_betas <- meth ## PACEProjector refers to this erroneously
pkg.est <- DunedinPACE::PACEProjector(meth)
us.est <- meffonym.dunedinpace.estimate(meth)

stopifnot(all(abs(pkg.est$DunedinPACE-us.est) < 1e-10))
