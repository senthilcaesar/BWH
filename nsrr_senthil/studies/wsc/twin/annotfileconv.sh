#!/bin/bash

# Folder location variable containing wsc original scoring files
DIR=/data/nsrr/working/wsc-scoring-annotations/twin

for f in ${DIR}/*allScore.txt
do

id=`cut -d'/' -f7 <<< $f | cut -d'.' -f1 | awk '{print substr($1,1, length($1)-8)}'` #get id for each file
  echo "Working on file with id ${id}"

 
tr -d '\r' <  $f | sed 's/ - /\t/g' | awk -F"\t" \
   ' function timeadd(t, dur){
        gsub(/SEC./,"",dur); 
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
             printf "# arousal_respiratory | Respiratory Arousal\n"; \
             printf "# arousal_lm | LM Arousal\n"; \
             printf "# hypopnea | Hypopnea  | min[num]\n"; \
             printf "# apnea_obstructive | Obstructive Apnea | min[num] \n"; \
             printf "# apnea_central | Central Apnea | min[num]\n"; \
             printf "# apnea_mixed | Mixed Apnea | min[num]\n"; \
             printf "# lm | LM\n"; \
             printf "# lights_off| Lights Out\n"; \
             printf "# lights_on | Lights On\n"; \
             printf "# paused | Paused\n"; \
             printf "# startrecording | Start Recording\n"; \
             printf "# misc | Miscellaneous | notes[txt]\n"; } \
         $2 == "DESATURATION" {split($3,a,":");split($4,b," ");split($5, c, " "); print "desat",".",".",$1,timeadd($1,a[2]),b[2]"|"c[2]} \
         $2 == "AROUSAL" && $4 == "SPONTANEOUS" {split($3,a, ":"); print "arousal_spontaneous",".",".",$1,timeadd($1,a[2]),"."} \
         $2 == "AROUSAL" && $4 == "RESPIRATORY EVENT" {split($3,a,":"); print "arousal_respiratory",".",".",$1,timeadd($1,a[2]),"."} \
         $2 == "AROUSAL" && $4 == "LM" {split($3, a, ":"); print "arousal_lm",".",".",$1,timeadd($1,a[2]),"."} \
         $2 == "RESPIRATORY EVENT" && $4 == "HYPOPNEA" {split($3,a,":");split($5,b," "); print "hypopnea",".",".",$1,timeadd($1,a[2]),b[2]} \
         $2 == "RESPIRATORY EVENT" && $4 == "OBSTRUCTIVE APNEA" {split($3,a,":");split($5,b," "); print "apnea_obstructive",".",".",$1,timeadd($1,a[2]),b[2]} \
         $2 == "RESPIRATORY EVENT" && $4 == "CENTRAL APNEA" {split($3,a,":");split($5,b," "); print "apnea_central",".",".",$1,timeadd($1,a[2]),b[2]} \
         $2 == "RESPIRATORY EVENT" && $4 == "MIXED APNEA" {split($3,a,":");split($5,b," "); print "apnea_mixed",".",".",$1,timeadd($1,a[2]),b[2]} \
         $2 == "LM" {split($3,a,":"); print "lm",".",".",$1,timeadd($1,a[2]),"."} \
         $2 == "LIGHTS OUT" {print "lights_off",".",".",$1,$1,"."} \
         $2 == "LIGHTS ON" {print "lights_on",".",".",$1,$1,"."} \
         $2 == "PAUSED" {print "paused",".",".",$1,$1,"."} \
         $2 == "START RECORDING" {print "startrecording",".",".",$1,$1,"."} \
         $2 != "START RECORDING" && $2 != "RESPIRATORY EVENT" && $2!= "DESATURATION" && $2!= "STAGE" && $2!= "AROUSAL" && $2!= "LM" && $2!= "LIGHTS OUT" && $2!= "LIGHTS ON" && $2!= "PAUSED" {print "misc",".",".",$1,$1,$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "$10" "$11" "$12" "$13} ' OFS="\t" > ${DIR}/${id}tmp.txt 
 
 # Replacing spaces with underscores for data rows - targeted for misc events
 cat ${DIR}/${id}tmp.txt | grep -e '^#' > ${DIR}/${id}.annot
 cat ${DIR}/${id}tmp.txt | grep -v '^#'  | tr -s ' ' | tr ' ' '_' | sed 's/_$//g' >> ${DIR}/${id}.annot
 rm ${DIR}/${id}tmp.txt
done

# Time add function within AWK command will return stop time
# provided start time and dur as input values

# Following manual changes are applied after running the above script
# For wsc-visit3-24698 → we have changed <NA> to "."
# For wsc-visit3-13061  → we have removed last line that was empty

# Apena events desat metadata is converted to min (to standardize across cohorts) 
