#!/bin/bash

i=1
while read -r case
do

    /PHShome/sq566/luna-base/luna s.lst ${i} -s STAGE eannot=/PHShome/sq566/nsrr/nap/eannots/${case}.eannot
    i=$((i+1)) 

done < cases.txt
