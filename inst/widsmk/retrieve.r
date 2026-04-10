## Herzog C, Jones A, Evans I, Raut JR, Zikan M, Cibula D, Wong A,
##  Brenner H, Richmond RC, Widschwendter M. Cigarette Smoking and
##  E-cigarette Use Induce Shared DNA Methylation Changes Linked to
##  Carcinogenesis. Cancer Res. 2024 Jun 4;84(11):1898-1914. doi:
##  10.1158/0008-5472.CAN-23-2957. PMID: 38503267; PMCID: PMC11148547.


load(url("https://github.com/chiaraherzog/WID.smk/raw/refs/heads/updated-sites/R/sysdata.rda"))
## [1] "df"           "models_exsmk" "models_nsmk"  "models_smk"   "mu"          
## [6] "sigma"        "WID_SMK_cpgs"

dat = WID_SMK_cpgs
dat = dat[dat$chr != "chrX",]

dat$set = paste0("WID_SMK_", dat$set)

WID_SMK_cpgs <- WID.smk:::WID_SMK_cpgs |>
    dplyr::filter(chr != "chrX")

for (name in dat$set) {
    in_model = dat$set == name
    model = data.frame(
        pred.var=dat$cg[in_model],
        coef=1/sum(in_model))
    write.csv(model, file=paste0(name,".csv"), row.names=F)
}


#' The WID.smk R package provides three functions:
#' 
#' - WID_SMK(meth)
#'   Calculates the four methylation scores from the input methylation matrix from
#'   a model with the same weight for each CpG site.
#'
#' - scale(scores)
#'   Standardizes calculated smoking scores by subtracting score-specific
#'   mu and dividing by a score-specific sigma calculated in the training dataset
#' 
#' - correct.ic(scores,ic)
#'   Regresses immune cell variation (ic) from calculated scores
#'   separately by smoking status.
#'
#' Important notes:
#' 
#' * 'scale' will not change the correlation of the scores with other variables.
#'
#' * 'correct.ic' is a standard linear adjustment except for the fact that
#'   it is performed by smoking group.
#'
#' The meffonym packages scores available for use as predictors and estimators.
#' 'scale' would have no effect on these tasks, and 'correct.ic' would
#' include the main exposure being estimated.
#' Consequently, both functions are omitted from the implemented here. 
