# DNA methylation indices of exposure and phenotype (meffonym)

## Installation 
```
library(devtools)
install_github("perishky/meffonym")
```

## Typical use

In the example below, we
show how to apply the Hannum clock (Hannum et al. 2013).

```
library(meffonym)

meth <- ... ## methylation matrix, rows=CpG site, columns=samples
age <- ... ## age for each sample

ret <- meffonym.score(meth, "hannum")

cor(ret$score, age)
## [1] 0.9342
```

## [inst/model.csv](Current list of models)

## For developers

Models can be added temporarily to the package using `meffonym.add.models()`.  However, these models will disappear as soon as R is exited or `meffonym` is reloaded.  To add a model so that it persists, create a new folder and add the model definition to [inst](inst) and add a new row for the model in [inst/models.csv](inst/models.csv).  The model definition is encoded in a csv file with two columns, `cpg` and `coef`.  One row can optionally provide a model intercept with the `cpg` name 'intercept'. 
