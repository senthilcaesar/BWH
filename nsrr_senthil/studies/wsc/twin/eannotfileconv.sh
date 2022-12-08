#!/bin/bash

DIR=/data/nsrr/working/wsc-scoring-annotations

for f in ${DIR}/twin/*allScore.txt
do
id=`cut -d'/' -f7 <<< $f | cut -d'.' -f1 | awk '{print substr($1,1, length($1)-8)}'` #get id for each file 
awk -F"\t" ' { print $2 } ' ${DIR}/twin/${id}allScore.txt | sort | uniq -c #extract staging to temporary .txt file
tr -d '\r' < ${DIR}/twin/${id}allScore.txt | awk -F"\t" \
    ' BEGIN { printf "# N1\n# N2\n# N3\n# R\n# W\n# ?\n" } \
      $2 == "START RECORDING" { t = $1 ; s = "" } \
      $2 == "STAGE - N1" && s != "" { print s , "." , t , $1 ; s = "N1" ; t = $1 } \
      $2 == "STAGE - N2" && s != "" { print s , "." , t , $1 ;s = "N2" ; t = $1 } \
      $2 == "STAGE - N3" && s != "" { print s , "." , t , $1 ;s = "N3" ; t = $1 } \
      $2 == "STAGE - R" && s != "" { print s , "." , t , $1 ;s = "R" ; t = $1 } \
      $2 == "STAGE - W" && s != "" { print s , "." , t , $1 ;s = "W" ; t = $1 } \
      $2 == "STAGE - NO STAGE" && s != "" { print s , "." , t , $1 ;s = "?" ; t = $1 } \
      $2 == "LIGHTS ON" && s != "" {  print s , "." , t , $1 } \
      $2 == "STAGE - N1" && s == "" { s = "N1" ; t = $1 } \
      $2 == "STAGE - N2" && s == "" { s = "N2" ; t = $1 } \
      $2 == "STAGE - N3" && s == "" { s = "N3" ; t = $1 } \
      $2 == "STAGE - R" && s == "" { s = "R" ; t = $1 } \
      $2 == "STAGE - W" && s == "" { s = "W" ; t = $1 } \
      $2 == "STAGE - NO STAGE" && s == "" { s = "?" ; t = $1 } '  OFS="\t" > tmp.txt
sed 's/\.[0-9][0-9]//g' < tmp.txt > tmp2.annot #create temporary .annot file from temporary .txt file
luna ${DIR}/twin/${id}.edf annot-file=tmp2.annot -o out.db -s STAGE #run luna's STAGE command on temporary .annot file
destrat out.db +STAGE -r E -v STAGE | awk ' NR != 1 { print $3 } ' | sed 's/NREM1/N1/g; s/NREM2/N2/g; s/NREM3/N3/g; s/wake/W/g; s/REM/R/g' >  ${DIR}/twin/${id}.eannot #create final .eannot file from luna output
done
