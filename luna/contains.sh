#!/bin/bash

arrText=($(cat sub.txt))


for i in "${!arrText[@]}"; do

    idx=$(( $i + 1 ))

    val=`luna s.lst $idx silent=T -s 'CONTAINS sig=${airflow}' | grep -i PRESENT | wc -l`
    
    if [ "$val" -eq "2" ]; then
        echo ${arrText[$i]} >> 2.txt
    
    elif [ "$val" -eq "3" ]; then
        echo ${arrText[$i]} >> 3.txt

    else
        echo ${arrText[$i]} >> none.txt

    fi

    #luna s.lst ${i} silent=T -s 'CONTAINS sig=${airflow}'

    #HAS_EEG=$?

    #if [[ ${HAS_EEG} -eq 0 ]]; then
    #    echo "Airflow channel present"
    #else
    #    echo "No found..."
    #fi

done
