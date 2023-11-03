#!/bin/bash

arrText=($(cat tmp/sub.txt))
rm -rf both.txt

for i in "${!arrText[@]}"; do
    idx=$(( $i + 1 ))
    val=`luna sl/xall.lst $idx silent=T -s 'CONTAINS sig=C3,C4' | grep -i PRESENT | awk '{sum+=$6;} END{print sum;}'`
    echo ${arrText[$idx]}   $val
    if [ "$val" -eq "2" ]; then
        echo ${arrText[$i]} >> bothChannelPresent.txt
    fi
done
