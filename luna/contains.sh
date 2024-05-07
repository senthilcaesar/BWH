#!/bin/bash

arrText=($(cat tmp/chin_already_referenced.txt))

for i in "${!arrText[@]}"; do
    idx=$(( $i + 1 ))
    val=`luna sl/chin_already_referenced.lst $idx silent=T -s 'CONTAINS sig=LChin_RChin,LChin_CChin,RChin_CChin' | grep "PRESENT" | awk '{sum+=$6;} END{print sum;}'`
    echo ${arrText[$idx]}   $val
    if [ "$val" -eq 3 ]; then
        echo ${arrText[$i]} >> threePresent.txt

    elif [ "$val" -eq 2 ]; then
        echo ${arrText[$i]} >> twoPresent.txt

    elif [ "$val" -eq 1 ]; then
        echo ${arrText[$i]} >> onePresent.txt
    fi
done
