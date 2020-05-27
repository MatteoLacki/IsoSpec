#!/bin/bash

cd $(git rev-parse --show-toplevel)
rm -rf IsoSpecR_*.tar.gz IsoSpecR.Rcheck
R CMD build IsoSpecR
R CMD check IsoSpecR_*.tar.gz --as-cran
