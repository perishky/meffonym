#!/bin/bash

## test adding a model
Rscript add-model.r

## verify Horvath clocks match calculator website
Rscript horvath.r

## compare clock outputs to estimage and methylclock output
## https://estimage.iac.rm.cnr.it/
## https://github.com/isglobal-brge/methylclock
Rscript estimage.r
