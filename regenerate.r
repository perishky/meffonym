library(devtools)
library(roxygen2)

setwd("..")
document("meffonym")

system("R CMD INSTALL meffonym")
