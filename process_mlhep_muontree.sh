#!/bin/bash

outputDir=/cmsuf/data/store/user/t2/users/klo/MLHEP/MuonTree/201211_mg_dyll/
for f in $(ls /cmsuf/data/store/user/t2/users/klo/MLHEP/Delphes/201211_mg_dyll/*.root) ;
do
    bf="$(basename -- ${f})"
    root -b -q 'MuonTreeProducer.C("'${f}'","'${outputDir}${bf}'")' ;
done
