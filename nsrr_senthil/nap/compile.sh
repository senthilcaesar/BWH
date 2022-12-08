#!/bin/bash

# usage
# bash nsrr/nap/compile.sh output-folder < input-roots

# for 'file.txt', basically 
#  cat nap/*/file.txt > output-folder/file.txt
# but adjusting for headers (taking only once)

if [[ $# -ne 1 ]]; then
    echo "expecting one argument, the output folder:"
    echo " usage: bash nsrr/nap/compile.sh output-folder < input-roots"
    echo "quitting"
    exit 1
fi

out=$1

mkdir -p ${out}

# loop through STDIN
cat | while read F; do

    ls nap/*/${F} > ${out}/.tmp.file.lst    
    l=`wc -l ${out}/.tmp.file.lst | awk ' { print $1 } ' `

    echo "found $l files"
    
    if [[ ${l} -gt 0 ]]; then

	echo "creating ${out}/${F}"

	# get header (**assumed to be constant across files**)
	#  which will be the case for Luna output
	
	head -1 `head -1 ${out}/.tmp.file.lst` > ${out}/${F}

	for d in `cat ${out}/.tmp.file.lst`
	do
	    awk ' NR != 1 ' ${d} >> ${out}/${F}
	done
	
    fi

    # check similar lines lengths
    c=`awk -F"\t" ' { print NF } ' ${out}/${F} | sort | uniq | wc -l | awk ' { print $1 } ' `
    if [[ $c -ne 1 ]]; then
	echo "problem: found unequal numbers of tab-delimited columns in ${out}/${F}"
	exit 1
    fi

    # clean up
    rm -rf ${out}/.tmp.file.lst

done



