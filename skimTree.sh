#!/bin/bash

outputDir=/cmsuf/data/store/user/t2/users/klo/Delphes/ALP_HToZaTo2l2g_M1/2020-06-02/MuonTreeProducer/
for f in $(ls /cmsuf/data/store/user/t2/users/klo/Delphes/ALP_HToZaTo2l2g_M1/2020-06-02/2020-06-02_ALP_HToZaTo2l2g_M1_13TeV_*.root) ;
do
    bf="$(basename -- ${f})"
    root -b -q 'MuonTreeProducer.C("'${f}'","'${outputDir}${bf}'")' ;
done
