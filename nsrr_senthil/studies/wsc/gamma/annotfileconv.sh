#!/bin/bash

DIR=/data/nsrr/working/wsc-scoring-annotations

for f in ${DIR}/gamma/*sco.txt
do
 id=`cut -d'/' -f7 <<< $f | cut -d'.' -f1 | awk '{print substr($1,1, length($1)-3)}'` #get id for each file 
 echo "Working with id ${id}"
 
 tr -d ' ' < $f | awk -F"\t" \
  'function timeadd(t, dur){ 
        start_time="jan 1 1970 + ";
        end="seconds";
        cmd="date -d ";
        space=" ";
        q="\"";
        format=" +\"%H:%M:%S.%2N\"" 
        (cmd q t space start_time dur end q format) | getline end_date;
        return end_date;
   }
   BEGIN { printf "# desat | SaO2 desaturations | min[num] drop[num]\n"; \
          printf "# arousal_spontaneous | Spontaneous Arousal\n"; \
          printf "# arousal_standard | Standard Arousal\n"; \
          printf "# arousal_respiratory | Respiratory Arousal\n"; \
          printf "# arousal_lm | LM Arousal\n"; \
          printf "# hypopnea | Hypopnea | min[num]\n"; \
          printf "# apnea_obstructive | Obstructive Apnea | min[num]\n"; \
          printf "# apnea_central | Central Apnea | min[num]\n"; \
          printf "# apnea_mixed | Mixed Apnea | min[num]\n"; \
          printf "# lm | LM\n"; \
          printf "# misc | Miscellaneous | notes[txt]\n";} \
      $5 == "SaO2" { print "desat" ,  ".", ".", $7, timeadd($7, $10), $8"|"$9} \
      $5 == "SponArousal"|| $5 == "sponarousal"|| $5 == "SPONArousal" { print "arousal_spontaneous" ,  ".", ".", $7, $7, "." } \
      $5 == "Arousal" { print "arousal_standard" , ".",  "." , $7 , $7, "." } \
      $5 == "RespA"|| $5 == "resparousal"|| $5 == "RespArousal"|| $5 == "RESPArousal" { print "arousal_respiratory" , ".", "." , $7 , $7, "."} \
      $5 == "LMA" { print "arousal_lm" , ".", "." , $7 , $7, "." } \
      $5 == "Hypopnea" || $5 == "CentralHypopnea"|| $5 == "Obst.Hypopnea" { print "hypopnea" , ".", "." , $7 , timeadd($7, $10), $8} \
      $5 == "OA"|| $5 == "ObsApnea"||$5 == "Obst.Apnea"|| $5 == "ObstApnea"|| $5 == "OBSApnea"|| $5 == "Apnea" { print "apnea_obstructive" , ".", "." , $7, timeadd($7, $10), $8 } \
      $5 == "CA"|| $5 == "CentralApnea" { print "apnea_central" , ".", "." , $7, timeadd($7, $10), $8} \
      $5 == "MA"|| $5 ==  "MixedApnea" { print "apnea_mixed" , ".", "." , $7, timeadd($7, $10), $8} \
      $5 == "LM" { print "lm" , ".", "." , $7, $7, "."} \
      $5 == "SnoreA" || $5 == "Snore" ||$5 == "PLME" ||$5 == "PLM" ||$5 == "BadSaO2Epoch" ||$5 == "BadECGEpoch" {print "misc",".",".",$7,$7,$5}'  OFS="\t" > ${DIR}/gamma/${id}.annot
done

# Time add function within AWK command will return stop time
# provided start time and dur as input values

# We have identified following files with issues
# negative dur values:
# wsc-visit1-10191
# wsc-visit1-11162
# wsc-visit1-41115
# wsc-visit1-75614
# wsc-visit1-82488
# wsc-visit1-89175
# wsc-visit2-12325
# wsc-visit2-64948

# missing start time:
# wsc-visit2-77724
