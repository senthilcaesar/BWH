#!/bin/bash

arrText=($(cat sub.txt))


for i in "${!arrText[@]}"; do

        idx=$(( $i + 1 ))

        echo $idx

        #val=`luna /data/purcell/projects/saps/sl/mnc_ssc.lst $idx silent=T -s 'CONTAINS sig=C4' | grep -i PRESENT | awk '{print $NF}'`
        val=`luna /data/purcell/projects/saps/sl/mnc_ssc.lst -i ${arrText[$i]} silent=T -s 'CONTAINS sig=C4,C3' | grep -i PRESENT | awk '{sum += $6} END {print sum}'`

        if [ "$val" -eq 2 ]; then
                echo ${arrText[$i]} >> c3c4.txt
        fi
done
