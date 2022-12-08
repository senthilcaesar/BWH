#!/bin/bash

# Folder location variable containing stages original annotation file
DIR=/data/nsrr/working/stages-annotations/

for sub in "BOGN" "GSBB" "GSDV" "GSLH" "GSSA" "GSSW" "MSMI" "MSNF" "MSQW" "MSTH" "MSTR" "STLK" "STNF" "MAYO"
do
  for i in ${DIR}/${sub}/*.csv
  do
    echo "cohort is ${i}"
    awk -F"," '{print $3}' ${i} >> ${sub}_all_annots.txt
  done
  cat ${sub}_all_annots.txt | sort | uniq -c > ${sub}_all_annots_count.txt
done
