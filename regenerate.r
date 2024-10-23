library(roxygen2)

setwd("..")
roxygenise("meffonym")

system("R CMD INSTALL meffonym")
