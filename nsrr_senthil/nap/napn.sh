#!/bin/bash

NAP_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo $NAP_DIR

# napn.sh {run-label} {folder} {n1} {n2} {alternate-config} 

# this file is called by nap.sh, and in turn (iteratively) calls nap1.sh
# processes individuals ${m} ${m}

# this may be running locally, or on a new node

# allow one-off (per batch) calling of a generic
# configuration script ${NAP_CONFIG_SOURCE} if specified

# run label
run=$1

# folder containing all input data
input=$2

# process from row 'n', to row 'm'
n=$3
m=$4

# secondary/alternative configuration data
conf2=$5


## --------------------------------------------------------------------------------
##
## Run any arbitrary setup code (with any new environment variables
## exported to this shell)
##
## --------------------------------------------------------------------------------


if [[ ! -z "${NAP_CONFIG_SOURCE}" ]]; then
    set -o allexport
    test -f ${NAP_CONFIG_SOURCE} && source ${NAP_CONFIG_SOURCE}
    set +o allexport
fi



## --------------------------------------------------------------------------------
##
## Run through individuals, calling nap1.sh on each
##
## --------------------------------------------------------------------------------

for j in `seq $n $m`
do

 # get ID from s.lst.run
id=`awk -F"\t" ' NR==j { print $1 }' j=${j} ${input}/s.lst.run`
 
 # check ID is not empty (e.g. if s.lst.run removed)
if [ -z ${id} ]
then
   echo "problem finding matchig ID"
   exit 1
fi
 
 # run primary NAP1 script for this indiv
echo bash ${NAP_DIR}/nap1.sh ${run} ${input} ${id} {$conf2} 

done


## --------------------------------------------------------------------------------
##
## All done
##
## --------------------------------------------------------------------------------
