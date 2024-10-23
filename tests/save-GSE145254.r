## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145254
library(geograbi) #https://github.com/yousefi138/geograbi

samples = geograbi.get.samples("GSE145254")
samples = cbind(samples, geograbi.extract.characteristics(samples, "characteristics_ch"))
samples = samples[,c("title","geo_accession","gender","disease status","age","race /ethnicity")]
colnames(samples) = c("Basename","id","Female","status","Age","ethnicity")
samples$Female = sign(samples$Female == "Female")
samples$Tissue = "Whole Blood"
write.csv(samples,file="GSE145254-samples.csv",row.names=F)

data = geograbi.get.data("GSE145254")
write.csv(data.frame(ProbeID=rownames(data),data), file="GSE145254-meth.csv", row.names=F)
          

