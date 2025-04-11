library(meffonym)

#remotes::install_github("chiaraherzog/WID.smk")
library(WID.smk)

wid_models = meffonym.models()[grep("wid-smk", meffonym.models())]
wid_models = sapply(wid_models, meffonym.get.model, simplify=F)
sites = unique(unlist(lapply(wid_models, function(x)names(x$coefs))))

meth = matrix(runif(length(sites)*20), nrow=length(sites))
rownames(meth) = sites
colnames(meth) = paste0("s",1:ncol(meth))

wid_out = sapply(WID_SMK(meth)[1:4], function(v) v)
mef_out = sapply(names(wid_models),  function(name) meffonym.score(meth,name)$score)

stopifnot(all(diag(cor(wid_out,mef_out)) > 1-1e-12))


