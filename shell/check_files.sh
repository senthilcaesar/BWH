#!/bin/bash

i=1
while read -r case
do

    FILE=/PHShome/sq566/nsrr/nap/nap/${case}/luna_suds_SOAP-tab.RData
    if test ! -f "$FILE"; then
        echo "$FILE does not exists."
    fi

done < missing.txt
