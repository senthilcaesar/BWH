#!/bin/bash

# Please run it on erisxdl1.research.partners.org
#
source /data/purcell/miniconda3/bin/activate
conda activate u-sleep
caselist=/data/purcell/projects/pops/edfs/v2-test/cases.txt

mapfile -t arr < ${caselist}

for sub in "${arr[@]}"
do
    echo ${sub}
    ut predict_one -f /data/purcell/projects/pops/edfs/v2-test/${sub}.edf -o /data/purcell/projects/pops/usleep/${sub}_usleep.hyp --model u-sleep:2.0 --overwrite
done
