library(readxl)

link <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8056641/bin/13148_2021_1055_MOESM4_ESM.xlsx"
filename <- basename(link)
if (!file.exists(filename))
    download.file(link, destfile=filename)

clocks <- list(
    "epic"=read_xlsx(filename, sheet="4A EPIC GA clock", skip=1),
    "450k"=read_xlsx(filename, sheet="4B 450K EPIC overlap clock", skip=1),
    "etd"=read_xlsx(filename, sheet="4C ETD-based clock",skip=1))

## intercept 293.455727634503

for (clock in names(clocks)) {
    with(clocks[[clock]], {
        write.csv(
            data.frame(pred.var=c("(Intercept)", as.character(CpG)), coef=c(0,Coefficient), stringsAsFactors=F),
            row.names=F,
            file=paste0(clock, ".csv"))
    })
}

unlink(filename)

