#!/bin/bash

## verify Horvath clocks match calculator website
Rscript horvath.r

## verify the DunedinPACE estimates match the R package
Rscript dunedinpace.r

## verify epismoker probabilities match originals
Rscript epismoker.r

## compare clock outputs to estimage and methylclock output
## https://estimage.iac.rm.cnr.it/
## https://github.com/isglobal-brge/methylclock
Rscript estimage.r


