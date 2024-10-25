# DNA methylation indices of exposure and phenotype (meffonym)

## Installation 
```
library(remotes)
install_github("perishky/meffonym")
```

## Typical use

In the example below, we
show how to apply the Hannum clock (Hannum et al. 2013).

```
library(meffonym)

meth <- ... ## methylation matrix, rows=CpG site, columns=samples
age <- ... ## age for each sample
female <- ... ## sex for each sample (1=female,0=male)

ret <- meffonym.score(meth, "hannum")

cor(ret$score, age)
## [1] 0.9342
```

## Special cases

Some models require additional inputs or have optional steps:

* The original epigenetic clock, **Horvath DNAmAge** ('horvath' model)
  includes an optional normalization step using the BMIQ algorithm and a
  DNA methylation 'gold standard'
  (http://labs.genetics.ucla.edu/horvath/dnamage/probeAnnotation21kdatMethUsed.csv).
  This can be achieved here by setting the `calibrate` parameter to `TRUE`
  in `meffonym.score()` or by using the `meffonym.horvath.age()` function.

```
meffonym.score(meth,"horvath",calibrate=TRUE)$score
== meffonym.horvath.age(meth)
```

* Similarly, **DunedinPACE** ('dunedinpace' model) allows a similar optional normalization
  step using a standard (from https://github.com/danbelsky/DunedinPACE).
  This can be achieved here by setting the `calibrate` parameter to `TRUE`
  in `meffonym.score()` or by using the `meffonym.dunedinpace.estimate()` function.

```
meffonym.score(meth, "dunedinpace", calibrate=TRUE)$score
== meffonym.dunedinpace.estimate(meth)
```

* The **Grimage** model ('grimage' and 'grimagev2') includes sex and age in addition to DNA methylation.
  These can be specified either by adding additional rows
  called 'female' (1 if female and 0 of male) and 'age' to the data matrix passed to `meffonym.score()`.
  Alternatively, the `meffonym.grimage()` function can be used.

```
meffonym.score(rbind(meth,female=female,age=age))$score
== meffonym.grimage(meth,female,age)
```

* A few clocks apply a log-transformation to linear model estimates.
  This is default behavior by `meffonym.score()` for
  the **Horvath DNAmAge** ('horvath'), **Skin and Blood** ('skin') and **PedBE** ('pedbe') models.
  This behavior can be changed for all clocks by passing an alterative transformation
  function to the `transform` parameter.


## Current list of models

[models.csv](inst/models.csv)

## For developers

Models can be added temporarily to the package using
`meffonym.add.models()`.  However, these models will disappear as soon
as R is exited or `meffonym` is reloaded.  To add a model so that it
persists, create a new folder and add the model definition to
[inst](inst) and add a new row for the model in
[inst/models.csv](inst/models.csv).  The model definition is encoded
in a csv file with two columns, `pred.var` and `coef`.  One row can
optionally provide a model intercept with the `pred.var` name
'intercept'.

## Models to add
* Shokhirev MN, Torosin NS, Kramer DJ, Johnson AA, Cuellar TL. CheekAge: a next-generation buccal epigenetic aging clock associated with lifestyle and health. Geroscience. 2024 Jun;46(3):3429-3443.
* Smith HM, Ng HK, Moodie JE, Gadd DA, McCartney DL, Bernabeu E, Campbell A, Redmond P, Taylor A, Page D, Corley J, Harris SE, Tay D, Deary IJ, Evans KL, Robinson MR, Chambers JC, Loh M, Cox SR, Marioni RE, Hillary RF. Methylome-wide studies of six metabolic traits. medRxiv [Preprint]. 2024 May 29:2024.05.29.24308103. doi: 10.1101/2024.05.29.24308103.
* Li X, Zhang Y, Gao X, Holleczek B, Schottker B, Brenner H. Comparative validation of three DNA methylation algorithms of ageing and a frailty index in relation to mortality: results from the ESTHER cohort study. EBioMedicine. 2021;74 doi: 10.1016/j.ebiom.2021.103686.
