#!/bin/bash

start_time=$(date +%F_%H-%M-%S)

NAP_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

## --------------------------------------------------------------------------------
##
## Arguments
##
## --------------------------------------------------------------------------------

if [[ $# -lt 2 ]]; then
    echo " error: expecting 2, 3 or 4 arguments"
    echo " usage: ./nap.sh run-label folder [id|m,n|.] [alternate-config-file]"
    exit 1
fi


#
# arg 1 : run label
#

run=$1


#
# arg 2 : root data (upload) folder
#

input=$2


#
# create primary NAP output folder
#

output=${input}/nap/

mkdir -p ${output}

#
# arg 3 : (optiona) individual ID to process (or . for all in folder/s.lst), i.e. often convenient
#       : if command delimited, then treat as to-from (integer rows)
#

extract_id="."
if [[ $# -ge 3 ]]; then
    extract_id=$3
fi

#
# arg 4 : (optional) alternate configuration file
#

conf2=""
if [[ $# -eq 4 ]]; then
    conf2=$4
fi

echo
echo "NAP core arguments:"
echo "  run-label       : [ $1 ]"
echo "  folder (input)  : [ $2 ]"
echo "  IDs             : [ $3 ]"
echo "  optional config : [ $4 ]"
echo

## --------------------------------------------------------------------------------
##
## Log info
##
## --------------------------------------------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');

echo "--------------------------------------------------------------------------------" 
echo "NAP v0.06 |  process started: ${dt} "                                            
echo "--------------------------------------------------------------------------------" 

echo
echo "Running nap.sh located in folder ${NAP_DIR}"
echo

## --------------------------------------------------------------------------------
##
## Configuration file(s)
##
## --------------------------------------------------------------------------------

# default configuration file; options set here can be overwritten if
# an optional config file is specified on command line

conf=${NAP_DIR}/default.conf

# read environment variables from default (and optional) config files
# allexport ensure that environment variables are exported to the current shell

set -o allexport

# echo all commands to log from this point

set -x

if [[ ! -z "${conf2}" ]]; then
    test -f ${conf2} && source ${conf2}
fi

test -f ${conf} && source ${conf}

set +x

set +o allexport

# ensure that Luna will always have '--log' 
NAP_LUNA_ARGS="${NAP_LUNA_ARGS} --log"

## --------------------------------------------------------------------------------
##
## Create Luna sample-list, if it doesn't already exist in the upload
##  i.e. allow users to supply their own s.lst file to link EDFs and annotations,
##       or to specify alternate IDs, if the EDFs/IDs/ANNOTs can not be automatically
##       linked via a standard 'build' 
##
## --------------------------------------------------------------------------------

echo 

# create sample list s.lst if it does not already exist
# the default Luna --build option is supplemented by any specofied additional flags
# as defined in NAP_SLST_BUILD_FLAGS: e.g. -ext=-nsrr.xml

if [[ ! -f "${input}/s.lst" ]]; then
  echo "Compiling sample list :  ${input}/s.lst "
  ${NAP_LUNA} --build ${input} ${NAP_SLST_BUILD_FLAGS} | sed 's/\.\///g' > ${input}/s.lst
  # check that this worked
  if [[ ! -f "${input}/s.lst" ]]; then
    echo "could not find sample-list ${input}, bailing"
    exit 1
  elif [ ${NAP_AWS_PORTAL_MODE} -eq 1 ] && [ command -v ${NAP_AWS_CLI} &> /dev/null ]; then
    s_path="$(echo ${NAP_DIR} |  cut -d'/' -f4-5)"
    echo "Copying sample list now"
    ${NAP_AWS_CLI} s3 --profile ${NAP_AWS_PROFILE} cp s.lst s3://nap-nsrr/${s_path}"/"
  else
   echo "AWS Portal mode is off or AWS cli is not installed, skipping upload of sample list from NAP"
  fi
else
   echo "Using existing sample list :  ${input}/s.lst " 
fi

slist="${input}/s.lst"


## --------------------------------------------------------------------------------
##
## Extract only a subset from the s.lst (e.g. for testing reasons / re-runs)
##
## Copy to ${slst}.run (temporary file) which will drive napn.sh
##
## --------------------------------------------------------------------------------

if [[ "${extract_id}" != "." ]]; then

    # if comma-delimited, assume m .. n (integers)
    fields=`echo ${extract_id} | awk -F"," '{ print NF }'`
    if [ ${fields} -eq 2 ]; then
	f1=`echo ${extract_id} | awk -F"," '{ print $1 }'`
	f2=`echo ${extract_id} | awk -F"," '{ print $2 }'`
	awk ' NR >= f1 && NR <= f2 ' f1=${f1} f2=${f2} ${slist} > ${slist}.run
    else
	f1=`echo ${extract_id} | awk -F"," '{ print $1 }'`
	awk ' $1 == id ' id=${extract_id} ${slist} > ${slist}.run
    fi
    
else
    echo "Processing entire ${slist}"
    cp ${slist} ${slist}.run
fi


## --------------------------------------------------------------------------------
##
## Test that we actually have some input to process
##
## --------------------------------------------------------------------------------

lines=`wc -l ${slist}.run | awk ' { print $1 } ' `
if [[ ${lines} -eq 0 ]]; then
    echo 
    echo "Quitting... no lines to process in ${slist}.run"
    echo 
    exit 1
fi



## --------------------------------------------------------------------------------
##
## Extract each EDF and run nap.sh, aiming for N-fold parallelism
##
## --------------------------------------------------------------------------------

if [[ ${NAP_JOBN} -gt 1 ]]; then
    echo "Assuming LSF bsub queue:"
    echo "  queue  = ${NAP_LSF_QUEUE}"
    #echo "  nodes  = ${NAP_LSF_NODES}"
    echo "  rusage = ${NAP_LSF_RUSAGE}" 
fi


##
## Split up input into 'm' jobs each of size 'n', aiming for ${njobs}-fold parallelization 
##

njobs=${NAP_JOBN}

l=`awk ' NF>0 ' ${slist}.run | wc -l`
n=`awk ' NF>0 ' ${slist}.run | wc -l | awk ' function ceiling(x){return x%1 ? int(x)+1 : x}  $1 == 0 { print 0 } $1 > 0 { print ceiling( $1 / k ) } ' k=${njobs} `
m=`awk ' NF>0 ' ${slist}.run | wc -l | awk ' function ceiling(x){return x%1 ? int(x)+1 : x}  $1 == 0 { print 0 } $1 > 0 { print ceiling( $1/ ceiling( $1 / k ) ) } ' k=${njobs} `

##
## Process each of 'm' jobs, each of of 'n' individual jobs
##

let i=1
for j in `seq 1 $m`
do
let i2=$i+$n-1

if [[ $i2 -gt $l ]]; then
i2=$l
fi

b=$(echo $j | awk '{printf("%05d", $1)}') 
echo " submitting batch ${b} ... individuals m=${i} to n=${i2}"

## build and execute command line (napn.sh, which will look for s.lst.run 

cmdline="${NAP_DIR}/napn.sh ${run} ${input} $i $i2 $conf2"

# if NAP_JOBN==1, then simple local job...
if [[ ${NAP_JOBN} -eq 1 ]]; then
 echo $cmdline | bash
else  # ... else use LSF
 bsub_cmd="bsub ${NAP_LSF_QUEUE} -n 1 \
                ${NAP_LSF_RUSAGE} \
                -o ${output}/tmp/batch${b}.out \
                -e ${output}/tmp/batch${b}.err"
  # Running bsub command stored in var 'bsub_cmd' using eval 
  echo ${bsub_cmd}
  echo ${cmdline}
  echo "$cmdline" | eval ${bsub_cmd}

fi

let i=$i+$n
done


## --------------------------------------------------------------------------------
##
## Collate output
##
## --------------------------------------------------------------------------------

## NAP_COPY_OUTPUT argument is helpful in use-cases (ex: Seven Bridges)
## where there is a need to write output to the home folder

if [[ ! -z "${NAP_COPY_OUTPUT}" ]]; then
  if [[ "${NAP_COPY_OUTPUT}" == "FILE" ]]; then
    echo "Creating NAP output as tar file with run name and start time info, in the home folder"
    output_file=~/${run}'_'${start_time}'_output.tar.gz'
    tar cvzf ${output_file} -C ${output} . && rm -R ${output}
  elif [[ "${NAP_COPY_OUTPUT}" == "DIRECTORY" ]]; then
    echo "Creating an output directory with run name and start time info, in the home folder"
    output_folder=~/${run}'_'${start_time}'_output'
    mkdir -p ${output_folder}
    mv ${output}* ${output_folder}
  fi
fi

## --------------------------------------------------------------------------------
##
## Wrap up
##
## --------------------------------------------------------------------------------

dt=$(date '+%d/%m/%Y %H:%M:%S');

echo "--------------------------------------------------------------------------------" 
echo "NAP v0.06 |  process finished : ${dt} "                                            
echo "--------------------------------------------------------------------------------" 
echo 


## --------------------------------------------------------------------------------
##
## All done
##
## --------------------------------------------------------------------------------

