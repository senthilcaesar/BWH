#!/bin/bash
caselist=$HOME/cases.txt

#BSUB -J NAP[1-599]%25
#BSUB -o %J-%I.out
#BSUB -e %J-%I.err

echo ${LSB_JOBINDEX}
case=`head -${LSB_JOBINDEX} ${caselist} | tail -1`

module load R/3.6.3
module load gcc/6.3.0

/PHShome/sq566/nsrr/nap/nap1.sh proj1 . ${case} ""


# bsub -q medium < test.sh
