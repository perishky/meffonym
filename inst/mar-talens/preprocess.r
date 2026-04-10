library(data.table)

metacols = c("name","tissue","target","publication","source","filename")
metadata = fread("metadata.csv",data.table=F)
stopifnot(all(colnames(metadata) == metacols))

## check that all files exist
for (i in 1:nrow(metadata)) {
  if (!file.exists(metadata$filename[i])) {
    cat("File does not exist: ",metadata$filename[i], "\n")
  }
}

## harmonize column names and remove rows coefficients equal to 0
new_cols = c("pred.var","coef")
old_cols1 = c("cpg_or_intercept","coefficient")
old_cols2 = c("Coefficient","CpG")
for (i in 1:nrow(metadata)) {
  filename = metadata$filename[i]
  data = fread(filename,data.table=F)
  changed = F
  if (all(colnames(data) != new_cols)) {
    changed =T
    if (all(colnames(data) == old_cols1)) {
      colnames(data) = new_cols
    } else if (all(colnames(data) == old_cols2)) {
      data = data[,rev(old_cols2)]
      colnames(data) = new_cols
    } else {
      stop(paste("Unknown column names in file: ",filename))
    }
  }
  idx = which(data$coef == 0)
  if (length(idx) > 0) {
    cat("removing ",length(idx),"rows from file: ",filename,"\n")
    data = data[-idx,]
    changed = T
  }
  if (changed) {
    cat("writing file: ",filename,"\n")
    fwrite(data,filename)
  }
}



  
 
