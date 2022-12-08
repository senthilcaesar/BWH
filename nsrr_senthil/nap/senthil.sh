#!/bin/bash


NAP_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# TODO
#  report channels dropped for out-of-range SR
#  flag and drop channels with zero variance


## --------------------------------------------------------------------------------
##
## Catch errors and clean up 
##
## --------------------------------------------------------------------------------

set -e

cleanup() {
    echo >> $LOG
    echo " *** encountered an error in NAP" >> $LOG
    echo " *** see ${ERR} for more details" >> $LOG
    echo >> $LOG

    echo >> $ERR
    echo " *** encountered an error in NAP *** " >> $ERR
    echo >> $ERR

    # leave a marker: 0 = fail
    echo -e "${id}\t0" > ${output}/${id}/nap.status
}

trap "cleanup" ERR


## --------------------------------------------------------------------------------
##
## Arguments
##
## --------------------------------------------------------------------------------

if [[ $# -lt 3 ]]; then
    echo " error: expecting 3+ arguments"
    echo " usage: ./nap1.sh run-label input-folder ID [config-file]"    
    cleanup
    exit 1
fi

#
# arg 1 : run label/group
#

run=$1

#
# arg 2 : root data (upload) folder
#

input=$2

#
# create primary NAP output folder, if doesn't already exist
#

output=${input}/nap

mkdir -p ${output}

#
# arg 3 : EDF ID in s.lst.run (i.e. to select one out of the 1+ EDFs in the upload folder)
#

id=$3

#
# arg 4 : (optional) alternate configuration file
#

conf2=""
if [[ $# -eq 4 ]]; then
    conf2=$4
fi


## --------------------------------------------------------------------------------
##
## Configuration file(s)
##
## --------------------------------------------------------------------------------

# default configuration file; options set here can be overwritten if
# an optional config file is specified on command line

conf="/PHShome/sq566/nsrr/nap/default.conf"

# echo all commands to log from this point

set -x

# read environment variables from default (and optional) config files
# allexport ensure that environment variables are exported to the current shell

set -o allexport

test -f ${conf} && source ${conf}

if [[ ! -z "${conf2}" ]]; then
    test -f ${conf2} && source ${conf2}
fi

set +o allexport



## --------------------------------------------------------------------------------
##
## Ensure that core files are wiped prior to run
##
## --------------------------------------------------------------------------------

if [ ${NAP_FRESH_START} -eq 1 ]; then
    rm -rf ${input}/nap/${id}
else
    rm -rf ${output}/${id}/tmp.*
    rm -rf ${output}/${id}/nap.log
    rm -rf ${output}/${id}/nap.err
    rm -rf ${output}/${id}/nap.status
    rm -rf ${output}/${id}/nap.issues

fi


mkdir -p ${input}/nap/${id}/


## --------------------------------------------------------------------------------
##
## Log info
##
## --------------------------------------------------------------------------------

LOG=${input}/nap/${id}/nap.log
ERR=${input}/nap/${id}/nap.err

dt=$(date '+%d/%m/%Y %H:%M:%S');


echo "--------------------------------------------------------------------------------" > $LOG
echo "NAP ${NAP_VERSION} | ${id} | process started: ${dt} "                             >> $LOG
echo "--------------------------------------------------------------------------------" >> $LOG

echo >> $LOG


## --------------------------------------------------------------------------------
##
## Create Luna sample-list, if it doesn't already exist in the upload
##  i.e. allow users to supply their own s.lst file to link EDFs and annotations,
##       or to specify alternate IDs, if the EDFs/IDs/ANNOTs can not be automatically
##       linked via a standard 'build' 
##
## --------------------------------------------------------------------------------

# create sample list s.lst if it does not already exist

if [[ ! -f "${input}/s.lst.run" ]]; then
    
    ${NAP_LUNA} --build ${input} ${NAP_SLST_BUILD_FLAGS} | sed 's/\.\///g' > ${input}/s.lst.run
    
    if [[ ! -f "${input}/s.lst.run" ]]; then
	echo
	echo "problem: could not find sample-list ${input}, bailing"
	cleanup
	exit 1
    fi
#else
#    echo "Using an existing sample list" >> $LOG
fi

echo "Logging for [ ${id} ] to [ ${input} ]" 2> $ERR

edf=`awk -v id=${id} -F"\t" ' $1 == id { print $2 } ' ${input}/s.lst.run`


# --------------------------------------------------------------------------------
#
# Log primary info
#
# --------------------------------------------------------------------------------

echo "Processing EDF ${id}" >> $LOG
echo "  - input folder:        ${input}" >> $LOG
echo "  - sample list:         ${input}/s.lst.run" >> $LOG
echo "  - EDF:                 ${edf}" >> $LOG
echo "  - primary output:      ${output}/${id}" >> $LOG
echo "  - new EDFs:            ${output}/${id}/data/" >> $LOG
echo "  - new annotations:     ${output}/${id}/annots/" >> $LOG
echo "  - issues:              ${output}/${id}/nap.issues" >> $LOG
echo "  - main log:            ${LOG}" >> $LOG
echo "  - error/verbose log:   ${ERR}" >> $LOG
echo >> $LOG
echo "Primary configuration values:" >> $LOG
echo "  - annot defs:          NAP_ANNOTS           = ${NAP_ANNOTS}" >> $LOG
echo "  - opt. annot defs:     NAP_OPT_ANNOTS       = ${NAP_OPT_ANNOTS}" >> $LOG
echo "  - signal defs:         NAP_SIGS             = ${NAP_SIGS}" >> $LOG
echo "  - opt. signal aliases: NAP_OPT_ALIASES      = ${NAP_OPT_ALIASES}" >> $LOG
echo "  - opt. signal defs:    NAP_OPT_SIGS         = ${NAP_OPT_SIGS}" >> $LOG 
echo "  - opt. signal group:   NAP_GROUP            = ${NAP_GROUP}" >> $LOG 
echo >> $LOG
echo "  - base signal defs:    NAP_BASE_SIGS        = ${NAP_BASE_SIGS}" >> $LOG
echo "  - opt. base sig. defs: NAP_BASE_OPT_SIGS    = ${NAP_BASE_OPT_SIGS}" >> $LOG
echo >> $LOG
echo "  - POPS directory:      NAP_POPS_DIR         = ${NAP_POPS_DIR}" >> $LOG
echo "  - POPS prefix:         NAP_POPS_VER         = ${NAP_POPS_VER}" >> $LOG
echo $NAP_R
echo $NAP_R_LIB
#echo "  - using CODA resources ${NAP_CODA_RESOURCE_DIR}" >> $LOG


## --------------------------------------------------------------------------------
##
## Some commands might be looking for individual-level meta-data
##
## --------------------------------------------------------------------------------

# i.e. commands expecting meta-data can add luna s.lst id=XYZ ${MD_INC} -s ' <commands> '
# i.e. below, add @${MD_INC}

MD_INC=""

if [ -n "${NAP_INDIV_METADATA}" ]; then
  MD_INC="vars=${NAP_INDIV_METADATA}" 
fi




## --------------------------------------------------------------------------------
##
## Template for adding new modules in NAP 
##
## --------------------------------------------------------------------------------

# echo "Running <command>..." >> $LOG
# <command> <input folder/id> <output> 2>> $ERR

#  ${input}            input folder (with 1 or more EDFs/annotations)
#  ${input}/s.lst      Original Luna SL
#  ${input}/s.lst.run  Row-subset (or same as) original SL
#  ${id}               ID of EDF (fileroot, id.edf)
#  ${output}           output folder
#  $ERR                standard error, redirected to:  ${output}/nap.err

# all output to ${output} should be tab-delimited files, adhering to naming
# convention described here: http://zzz.bwh.harvard.edu/luna/merge/merge/#definitions

# General NAP assumptions:
#  1) A study is a single EDF
#  2) Any annotations are in .annot, .eannot, NSRR-XML format or EDF+ Annotations format


## --------------------------------------------------------------------------------
##
## Primary outputs
##
## --------------------------------------------------------------------------------

## NSRR-harmonized EDF and annotations (.annot and .xml)
#    ${output}/${id}/data/${id}-harm.edf
#    ${output}/${id}/annots/harm.annot

## NSRR base EDF )
#    ${output}/${id}/data/${id}-base.edf

## Issues for this EDF (e.g. polarity, line noise, etc)
#     ${output}/${id}/data/nap_issues.txt 

## NAP derived tables/figures
#    ${output}/${id}

## NAP derived annotations (.annot)
#    ${output}/${id}/annots

## Output file name conventions
# <domain>_<group>_<command>.txt
# <domain>_<group>_<command>_F1_F2.txt
# <domain>_<group>_<command>_F1-C3_F2-22.txt
# <domain>_<group>_<command>_F1-C3_F2--22.txt  // allowed --> F2 is '-22'



## --------------------------------------------------------------------------------
##
## NAP procedures
##
## --------------------------------------------------------------------------------

# Harmonized EDFs
#  - harmonized channel labels / drop unrecognized channels
#  - harmonized annotationslabels / drop unrecognized annotations
#  - all annotations output as .annot and .xml
#  - staging output as .eannot
#  - EDF physical units standardized (based on EDF header)
#  - EDF fixed 1 second EDF record size
#  - EDF epoch/staging alignment enforced (if staging is present)

# Base EDFs
#  - reduce to a set of core channels
#  - re-reference EEG as needed
#  - bandpass filter EEG
#  - attempt to fix EEG polarity
#  - user can specify canonical signal list,

# SOAP/POPS
#  - automated staging based on base.edf
#  - uses csCEN alone by default
#  - user can also specify
#  - (todo: EOG and EMG staging)
#  - (todo: second round prior probability weighting)

# Issues (NI:) NAP Issues
#  - absence of staging data / poor coverage, i.e. truncated staging
#  - gross artifact (EEG only, by epoch and channel)
#  - flag potential EEG polarity issues (negative / ambiguous channels)
#  - flag likely ECG contamination
#  - flag likely line noise / EDFs with peaks in the <35 Hz range
#  - flag if likely units issue in the EEG (i.e. large values, or drift)
#  - flag if is low SOAP consistency

## --------------------------------------------------------------------------------
##
## Steps:-
##   - determine whether staging is present
##
##   - basic QC on all channels (original EDF) incl/ MTM and HEADERS
##   - run coda1 to produce summaries
##
##   - create harmonized EDF/annotations
##
##   - create base EDF/annotations
##   - EEG polarity check
##   - SOAP/POPS
##   - spindles/SO
##   - respiratory analyses (incl hypoxic burden)
##
## --------------------------------------------------------------------------------




## --------------------------------------------------------------------------------
##
## Create primary output folders: (or ensure they are blank)
##   data      EDF & original annots, e.g. manual staging
##   annots    NAP-derived annots, e.g. spindles, respiratory events, automated staging
##   issues    individual level issue reporting (single file nap_issues.txt)
##
## --------------------------------------------------------------------------------

rm -rf ${output}/${id}/data ${output}/${id}/annots
mkdir -p ${output}/${id}/data ${output}/${id}/annots

ISSUES=${output}/${id}/nap.issues
rm -rf ${ISSUES}

## --------------------------------------------------------------------------------
##
## For parsing output, use domain+group+file+{factors}.txt convention [ for dmerge tool ]
## For convenience, create some variables here, using tt-prepend
##
## --------------------------------------------------------------------------------

# Basic EDF : HEADERS, CANONICAL, FLIP
dom_core="tt-prepend=luna_core_"

# Signal stats : STATS, SIGSTATS
dom_stats="tt-prepend=luna_stats_"

# HYPNO
dom_macro="tt-prepend=luna_macro_"

# SOAP & POPS
dom_suds="tt-prepend=luna_suds_"

# spectral (PSD & MTM)
dom_spec="tt-prepend=luna_spec_"

# spindles and slow oscillations: SPINDLES SO
dom_spso="tt-prepend=luna_spso_"


    
## --------------------------------------------------------------------------------
##
## Check for presence of mapping files
##
## --------------------------------------------------------------------------------

if [ ! -f "${NAP_SIGS}" ]; then
  echo "Could not find a harmonized canonical signal definition file at : ${NAP_SIGS}" >> $LOG
  echo "Please set the NAP_SIGS variable" >> $LOG
  echo "NAP requires a canonical signal definition... bailing" >> $LOG
  cleanup
  exit 1
fi

if [ ! -f "${NAP_BASE_SIGS}" ]; then
  echo "Could not find a base canonical signal definition file at : ${NAP_BASE_SIGS}" >> $LOG
  echo "Please set the NAP_BASE_SIGS variable" >> $LOG
  echo "NAP requires a canonical signal definition... bailing" >> $LOG
  cleanup
  exit 1
fi



## --------------------------------------------------------------------------------
##
## Harmonzied and base channel reports
##
## --------------------------------------------------------------------------------

echo >> $LOG
echo "Mapping channels & annotations" >> $LOG
echo "------------------------------" >> $LOG

##
## Channels: make a toy version of the harmonized EDF (1 epoch)
##           and then use that to make the base EDF
## with drop-originals, this also reports a list of the unused signals
##

echo "Mapping harmonized channel labels, and saving a 1-epoch harmonized EDF..." >> $LOG

harmsl=${output}/${id}/harm.lst

if [ ${NAP_MAKE_DATA} -eq 1 ]; then
        
    rm -rf ${harmsl}

    ${NAP_LUNA} ${input}/s.lst.run id=${id} ${NAP_LUNA_ARGS} \
		csfile=${NAP_OPT_SIGS},${NAP_SIGS} \
		group=${NAP_GROUP} \
		outdir=${output}/${id} \
		@${NAP_OPT_ALIASES} \
		-t ${output} ${dom_core} \
		tt-append=_EDF-HARM \
        -s 'CANONICAL file=${csfile} group=${group} drop-originals &
            MASK epoch=1 & RE & 
            WRITE edf-dir=${outdir}/data/ edf-tag=TEST123XYZ drop-originals sample-list=${outdir}/harm.lst' \
        > /dev/null 2>> $ERR
fi
#
#
###
### Harmonized --> base EDF mapping report
###
#
#echo "Mapping base channel labels from the (single-epoch) harmonized EDF..." >> $LOG
#
#echo "Harm file" >> $LOG
#echo ${harmsl} >> $LOG
#
#if [ ${NAP_MAKE_DATA} -eq 1 ]; then
#
#    ${NAP_LUNA} ${harmsl} id=${id} ${NAP_LUNA_ARGS} \
#		csfile=${NAP_BASE_OPT_SIGS},${NAP_BASE_SIGS} \
#		group=${NAP_GROUP} \
#		-t ${output} ${dom_core} \
#		tt-append=_EDF-BASE \
#		-s 'CANONICAL file=${csfile} group=${group} drop-originals check ' 2>> $ERR
#    
#    # clean up the data/ folder
#    rm -rf ${harmsl}
#    rm -rf ${output}/${id}/data/*TEST123XYZ*
#fi
#
#
### --------------------------------------------------------------------------------
###
### Annotation reorts
###
### --------------------------------------------------------------------------------
#
## as certain features (e.g. split/combine class instances etc only work on
## .annot (and not XML) annots. let's first make a single .annot file and
## subsequently add skip-sl-annots=T annot-file=<this>
## do not do any conversion at this step, nb. keep any spaces: but enfore tab-delimited inputs
## write as HH:MM:SS so that things work if we remove epochs (trim because of misaligned stages)
#
#echo "Collating all annotations in a single file a.annot : ${output}/${id}/annots/a.annot" >> $LOG
#
#${NAP_LUNA} ${input}/s.lst.run id=${id} ${NAP_LUNA_ARGS} nsrr-remap=F tab-only=T silent=T \
#	    -s WRITE-ANNOTS no-specials hms file=${output}/${id}/annots/a.annot 2>> $ERR
#
#    
## now take this a.annot and create the harmonized annotation set, tracking what is mapped
## using hh:mm:ss formats
#
#echo "Creating a harmonized annotation file harm.annot : ${output}/${id}/annots/harm.annot" >> $LOG
#
#${NAP_LUNA} ${input}/s.lst.run id=${id} ${NAP_LUNA_ARGS} \
#	    skip-sl-annots=T skip-edf-annots=T \
#	    annot-file=${output}/${id}/annots/a.annot \
#	    @${NAP_ANNOTS} \
#	    @${NAP_OPT_ANNOTS} \
#	    annot-whitelist=T \
#	    align-annots=${NAP_STAGING_ANNOTS} \
#	    tab-only=T \
#	    odir=${output}/${id} \
#	    -t ${output} ${dom_core} \
#	    -s 'ALIASES & TAG A/MAPPED & ANNOTS & WRITE-ANNOTS no-specials hms file=${odir}/annots/harm.annot' 2>> $ERR
#
## now this harmonized .annot file has been created, use it below (instead of the SL-attached or EDF+ ones)
#USE_HARM_ANNOTS="skip-sl-annots=T skip-edf-annots=T annot-file=${output}/${id}/annots/harm.annot"
#
## also, track what went unmapped
#    
#${NAP_LUNA} ${input}/s.lst.run id=${id} ${NAP_LUNA_ARGS} \
#	    skip-sl-annots=T skip-edf-annots=T \
#	    annot-file=${output}/${id}/annots/a.annot \
#	    @${NAP_ANNOTS} \
#	    @${NAP_OPT_ANNOTS} \
#	    annot-unmapped=T \
#	    tab-only=T \
#	    -t ${output} ${dom_core} \
#	    -s 'TAG A2/UNMAPPED & ANNOTS' 2>> $ERR
#
## nb. using A and A2 meant we did not overwrite stuff; but we're only interested here
## in the combined table of annotation labels, so construct that here:
#
#if [  -f ${output}/${id}/luna_core_ANNOTS_ANNOT_A2.txt ]; then
#    echo "Reporting unmapped annotations..." >> $LOG
#    awk ' NR != 1 ' ${output}/${id}/luna_core_ANNOTS_ANNOT_A2.txt >> ${output}/${id}/luna_core_ANNOTS_ANNOT_A.txt
#    awk ' NR != 1 ' ${output}/${id}/luna_core_ANNOTS_ANNOT_INST_A2.txt >> ${output}/${id}/luna_core_ANNOTS_ANNOT_INST_A.txt
#    awk ' NR != 1 ' ${output}/${id}/luna_core_ANNOTS_ANNOT_INST_T_A2.txt >> ${output}/${id}/luna_core_ANNOTS_ANNOT_INST_T_A.txt
#    # & some clean up
##    rm -rf ${output}/${id}/luna_core_ANNOTS_ANNOT_A2.txt
##    rm -rf ${output}/${id}/luna_core_ANNOTS_ANNOT_INST_A2.txt
##    rm -rf ${output}/${id}/luna_core_ANNOTS_ANNOT_INST_T_A2.txt
#else
#    echo "All original annotations were mapped, good" >> $LOG
#fi
#
#
### --------------------------------------------------------------------------------
###
### Only harmonize signals/labels?  (i.e. versus a full analysis?) 
###
### --------------------------------------------------------------------------------
#
#
#if [ ${NAP_HARMONIZE_ONLY} -eq 1 ]; then
#    
#    # compile tables from the above (i.e. an early call to coda1.R
#    # this should skip anything that has not yet been created (e.g. MTM)
#    # so should be okay to run now, if not doing the full NAP
#    
#    if [ ${NAP_DO_CODA1} -eq 1 ]; then
#	echo "Compiling tables (coda1)..." >> $LOG
#	${NAP_R} ${NAP_DIR}/coda1.R ${NAP_DIR} ${NAP_CODA_RESOURCE_DIR} ${output}/${id} ${NAP_R_LIB} 0 >> $ERR 2>&1
#    fi
#    
#    echo "As NAP_HARMONIZE_ONLY is set, leaving NAP now" >> $LOG
#
#    # leave a marker: 1 = success
#    echo -e "${id}\t1" > ${output}/${id}/nap.status
#    # bye bye
#    exit 0
#    
#fi
#
#
#
### --------------------------------------------------------------------------------
###
### Original EDF  
###
### --------------------------------------------------------------------------------
#
#echo >> $LOG
#echo "Original EDF summary" >> $LOG
#echo "--------------------" >> $LOG
#
## Basic report (number of signals, annotations and duration)
#${NAP_LUNA} ${input}/s.lst.run ${NAP_LUNA_ARGS} id=${id} verbose=T -s DESC 2> ${output}/${id}/tmp.headers
#ORIG_DUR=`awk ' $1 == "duration:" { print $2 } ' ${output}/${id}/tmp.headers`
#ORIG_SEC=`awk ' $1 == "duration:" { print $4 } ' ${output}/${id}/tmp.headers`
#ORIG_NS=`awk ' $1 == "signals:" { print $2 } ' ${output}/${id}/tmp.headers`
#ORIG_NA=`awk ' $3 == "instance(s)" && $4 == "(from" ' ${output}/${id}/tmp.headers | wc -l | awk ' { print $1 } ' `
#echo "Detected $ORIG_NS channels, $ORIG_NA annotation classes & duration $ORIG_DUR" >> $LOG
#
### --------------------------------------------------------------------------------
###
### Test for the presence of existing (manual) staging in the original EDF
### (nb. temporarily turn off the trap)
###
### --------------------------------------------------------------------------------
#
#trap '' ERR; set +e
#${NAP_LUNA} ${input}/s.lst.run ${NAP_LUNA_ARGS} id=${id} ${NAP_ORIG_ANNOTS} \
#            ${USE_HARM_ANNOTS} \
#            -s 'CONTAINS stages' > ${output}/${id}/tmp.contains 2>> $ERR
#
#STAGES_ABSENT=$?
#set -e; trap "cleanup" ERR
#
#if [[ ${STAGES_ABSENT} -gt 0 ]]; then
#    echo "No existing stage information detected - will skip HYPNO, POPS, etc..." >> $LOG
#else
#    echo "Existing staging annotations detected" >> $LOG
#fi
#
## additional check on variability of staging
## i.e. if all W, or all N2, then we don't want to attempt SOAP, etc
#
#if [[ ! -f "${output}/${id}/tmp.contains" ]]; then
#    echo "internal NAP error: could not detect ${output}/${id}/tmp.contains" >> $LOG
#    cleanup
#    exit 1
#fi
#
#UNIQ_STAGES=`awk ' $5 == "UNIQ_STAGES" { print $6 } ' ${output}/${id}/tmp.contains`
#STAGE_COUNTS=`awk ' $5 == "STAGE_COUNTS" { print $6 } ' ${output}/${id}/tmp.contains`
#
#echo "Existing stage epoch counts: ${STAGE_COUNTS}" >> $LOG
#
#if [[ ${UNIQ_STAGES} -lt 2 ]]; then
#    echo "No variability in existing sleep staging: treating as if absent" >> $LOG
#    STAGES_ABSENT=1
#fi
#
#rm -rf ${output}/${id}/tmp.contains
#
## finally, ask:
## do original staging annotations align w/ EDF epochs?
#
#if [ ${STAGES_ABSENT} -eq 0 ]; then
#
#    # test for CONF in stages in the original EDF...
#    HAS_CONF=`${NAP_LUNA} ${input}/s.lst.run id=${id} ${NAP_LUNA_ARGS} skip-sl-annots=T skip-edf-annots=T annot-file=${output}/${id}/annots/a.annot  -s HYPNO | awk ' $5 == "CONF" { print $6 } '`    
#
#    if [ -z ${HAS_CONF} ]; then
#	echo "problem defining HAS_CONF: see nap1.sh" >> $LOG
#	cleanup
#	exit 1
#    fi
#    
#    if [ ${HAS_CONF} -gt 0 ]; then
#	echo "  *** original staging does not align with EDF epochs ($HAS_CONF epochs)" >> $LOG
#    else
#	echo "Original staging aligns with EDF epochs, good" >> $LOG
#    fi
#    
#fi
#
#
#
### --------------------------------------------------------------------------------
###
### Basic QC (on original EDF, w/out any channel aliasing)
###
### --------------------------------------------------------------------------------
#
#if [ ${STAGES_ABSENT} -eq 0 ]; then
#    echo "Running HYPNO..." >> $LOG
#    ${NAP_LUNA} ${input}/s.lst.run id=${id} ${NAP_LUNA_ARGS} \
#                ${USE_HARM_ANNOTS} \
#		-t ${output} ${dom_macro} \
#                -s HYPNO 2>> $ERR    
#fi
#
#echo "Running HEADERS..." >> $LOG
#${NAP_LUNA} ${input}/s.lst.run id=${id} ${NAP_LUNA_ARGS} -t ${output} ${dom_core} -s HEADERS 2>> $ERR
#
#echo "Running SIGSTATS..." >> $LOG
#${NAP_LUNA} ${input}/s.lst.run id=${id} ${NAP_LUNA_ARGS} -t ${output} ${dom_stats} -s 'SIGSTATS epoch sr-over=50' 2>> $ERR
#
#echo "Running STATS..." >> $LOG
#${NAP_LUNA} ${input}/s.lst.run id=${id} ${NAP_LUNA_ARGS} -t ${output} ${dom_stats} -s 'STATS epoch sr-under=50' 2>> $ERR
#
#
### --------------------------------------------------------------------------------
###
### MTM spectrograms (original EDF, all channels with SR > 50)
###
### --------------------------------------------------------------------------------
#
#if [ ${NAP_DO_MTM} -eq 1 ]; then
#    
#    echo "Running MTM spectrograms (EEG/EMG/EOG channels with SR of at least ${NAP_MTM_MIN_SAMPLE_RATE} Hz)..." >> $LOG
#    
#    ${NAP_LUNA} ${input}/s.lst.run id=${id} ${NAP_LUNA_ARGS} \
#		fs=${NAP_MTM_MIN_SAMPLE_RATE} \
#		-t ${output} ${dom_spec} \
#		-s 'MTM sig=${eeg},${eog},${emg} segment-sec=30 segment-inc=30 min=0.5 max=25 nw=15 epoch sr=${fs}' 2>> $ERR
#
#fi
#
#
### --------------------------------------------------------------------------------
###
### Compile SIGSTATS and MTM spectograms results into R dataframes for viewing (luna/shiny)
###
### --------------------------------------------------------------------------------
#
## i.e.     path/to/folder/text.txt
## becomes  path/to/folder/text.txt-tab.RData
##      or  path/to/folder/text.txt-fig.RData
## luna-shiny then automatically loads any *-tab.RData and *-fig.RData files
# 
## See example of MTM plot for how to save / attach images, points to .png files, 
## which can either be created in coda1.R based
##  on summary stats, or indpendently (in which case a *-fig.RData file is created
##  which just points to the existing .png
#
#if [ ${NAP_DO_CODA1} -eq 1 ]; then
#    echo "Running coda(1) - compiling SIGSTATS/MTM output into RData files..." >> $LOG
#    ${NAP_R} ${NAP_DIR}/coda1.R ${NAP_DIR} ${NAP_CODA_RESOURCE_DIR} ${output}/${id} ${NAP_R_LIB} ${NAP_DO_MTM} >> $ERR 2>&1
#fi
#
#
### --------------------------------------------------------------------------------
###
### Harmonized EDF
###
### --------------------------------------------------------------------------------
#
#echo >> $LOG
#echo "Building the harmonized EDF" >> $LOG
#echo "---------------------------" >> $LOG
#
## get EDF record size
#
#REC_DUR=`${NAP_LUNA} ${input}/s.lst.run ${NAP_LUNA_ARGS} 1 -s HEADERS | awk ' $3 == "." && $5 == "REC_DUR" { print $6 } '`
#
#if [ -z ${REC_DUR} ]; then
#    echo "problem defining REC_DUR: see nap1.sh" >> $LOG
#    cleanup
#    exit 1
#fi
#
#
##
## the harmonized EDF will have EDF record size of exactly 1 second
## this means, if making a new EDF, we need to drop channels with
## a sample rate that cannot be represented in a 1-second block
## (i.e. fractional Hz)
## if the sample rate is fractional but < 1, we can first try a ZOH
## up-sampling to get an exact 1 Hz rate (i.e. 1/30 -> expand by factor of 30)
##
## also, drop any signals with SR > 500 Hz (whether we make a new EDF or no) 
## 
#
#SL_1SEC=${input}/s.lst.run
#
#if [ $REC_DUR -eq 1 ]; then
#    echo "EDF record duration is already 1 second, good" >> $LOG
#else
#    echo "Re-writing EDF w/ a 1-sec record duration (instead of $REC_DUR)" >> $LOG
#fi
#
#echo "Re-writing EDF to drop any channels with disallowed sample rates" >> $LOG
#
#onesl=${output}/${id}/rs1.lst
#rm -rf ${onesl}
#
#${NAP_LUNA} ${input}/s.lst.run id=${id} ${NAP_LUNA_ARGS} \
#	    @${NAP_OPT_ALIASES} \
#	    outdir=${output}/${id} \
#	    f1=${NAP_HARM_SR_LWR} f2=${NAP_HARM_SR_LWR} \
#	    force-edf=T skip-edf-annots=T \
#	    -t ${output} ${dom_core} \
#	    -s 'ZOH osr=1 sr=1 &
#       	        ENFORCE-SR dur=1 range=${f1},${f2} &
#                RECORD-SIZE dur=1 edf-dir=${outdir}/data/ edf-tag=rs1 sample-list=${outdir}/rs1.lst' 2>> $ERR
#
## update $SL_1SEC to point to this (instead of the original) 
#SL_1SEC=${onesl}
## track this intermediate EDF (to be deleted at the end of the script)
#TMP_EDF=`cut -f2 ${onesl}`
#
####  TODO - note that can only change EDF record size for CONTINUOUS
####   EDFs.. need to test for this first / do we want to enforce that NAP
####   can only handle EDF and EDF+C?
#
## harmonized sample-list location
#harmsl=${output}/${id}/harm.lst
#rm -rf ${harmsl}
#
#
#
## Build the harmonized EDF (ensuring we have a 1-second EDF record
## size EDF, in $SL_1SEC)
#
## If we have manual staging present, then also ensure we align any
## staging annotations (these will be based on the mapped terms) to the
## nearest second, if fractional i.e. we assume this set {
## ${NAP_STAGING_ANNOTS} = W,N1,N2,N3,R,?,L,M,U } defines the universe
## of staging epoch annotations we can make this configurable if needed
## (but these are the mapped terms, so we should not need to)
#
## nb. we are reading the harm.annot file here, skipping any SL annotations
#
#if [ ${STAGES_ABSENT} -eq 0 ]; then
#
#    # with staging: align on epochs, and trim to include only epochs w/ staging
#    # as well as the channel remapping; nb use of single second epoch duration here
#    # (to align w/ any possible annotation realignment & EDF record size of 1.0s)
#
#    ${NAP_LUNA} ${SL_1SEC} id=${id} ${NAP_LUNA_ARGS} \
#    		csfile=${NAP_OPT_SIGS},${NAP_SIGS} \
#		group=${NAP_GROUP} \
#		outdir=${output}/${id} \
#		${USE_HARM_ANNOTS} \
#		stgannots=${NAP_STAGING_ANNOTS} \
#		cmp=${NAP_EDFZ} \
#		force-edf=T skip-edf-annots=T \
#		-t ${output} ${dom_core} \
#		-s 'CANONICAL file=${csfile} group=${group} drop-originals &
#		    EPOCH dur=1 & MASK all & MASK unmask-if=${stgannots} & RE &
#           	    WRITE edf-dir=${outdir}/data/ edfz=${cmp} edf-tag=harm with-annots sample-list=${outdir}/harm.lst' 2>> $ERR
#else
#
#    # just do the channel remapping step here
#
#    ${NAP_LUNA} ${SL_1SEC} id=${id} ${NAP_LUNA_ARGS} \
#    		csfile=${NAP_OPT_SIGS},${NAP_SIGS} \
#		group=${NAP_GROUP} \
#		outdir=${output}/${id} \
#		${USE_HARM_ANNOTS} \
#		force-edf=T skip-edf-annots=T \
#		align-annots=${NAP_STAGING_ANNOTS} \
#		cmp=${NAP_EDFZ} \
#		-t ${output} ${dom_core} \
#		-s 'CANONICAL file=${csfile} group=${group} drop-originals &
#           	    WRITE edf-dir=${outdir}/data/ edfz=${cmp} edf-tag=harm with-annots sample-list=${outdir}/harm.lst' 2>> $ERR
#    
#fi
#
#echo "Harmonized EDF         : ${output}/${id}/data/${id}-harm.edf" >> $LOG
#echo "Harmonized annotations : ${output}/${id}/annots/harm.annot" >> $LOG
#echo "Harmonized SL          : ${output}/${id}/harm.lst" >> $LOG
#
#
### --------------------------------------------------------------------------------
###
### Check whether staging & epochs align in the harmonized EDF
###
### --------------------------------------------------------------------------------
#
## test for CONF in stages...
#
#if [ ${STAGES_ABSENT} -eq 0 ]; then
#
#    HAS_CONF=`${NAP_LUNA} ${harmsl} 1 -s HYPNO | awk ' $5 == "CONF" { print $6 } '`
#
#    if [ -z ${HAS_CONF} ]; then
#	echo "problem defining HAS_CONF: see nap1.sh" >> $LOG
#	cleanup
#	exit 1
#    fi
#
#    if [ ${HAS_CONF} -gt 0 ]; then
#	echo " ** warning: staging & epochs do not align in harmonized EDF (should not happen)" >> $LOG
#    else
#	echo "Staging and EDF epochs align in the harmonized EDF, good" >> $LOG
#    fi
#fi
#
#
### --------------------------------------------------------------------------------
###
### Report some basic stats on the harmonized EDF (# signals, annots, & duration)
###
### --------------------------------------------------------------------------------
#
#${NAP_LUNA} ${harmsl} id=${id} verbose=T -s DESC 2> ${output}/${id}/tmp.headers
#HARM_DUR=`awk ' $1 == "duration:" { print $2 } ' ${output}/${id}/tmp.headers`
#HARM_SEC=`awk ' $1 == "duration:" { print $4 } ' ${output}/${id}/tmp.headers`
#HARM_NS=`awk ' $1 == "signals:" { print $2 } ' ${output}/${id}/tmp.headers`
#HARM_NA=`awk ' $3 == "instance(s)" && $4 == "(from" ' ${output}/${id}/tmp.headers | wc -l | awk ' { print $1 } ' `
#echo "Detected $HARM_NS channels, $HARM_NA annotation classes & duration $HARM_DUR" >> $LOG
#
#
### --------------------------------------------------------------------------------
###
### Store ORIG_DUR, ORIG_NS, ORIG_NA versus HARM_NA etc
###
### --------------------------------------------------------------------------------
#
## nb. 'let z=${x}-${y}' can trigger the ERR trap, as let returns non-zero exit code
## if it evaluates to 0...  so have to use awk instead to do the math...  urgh
#
##let DROP_SEC=${ORIG_SEC}-${HARM_SEC}
#DROP_SEC=`echo "${ORIG_SEC} ${HARM_SEC}" | awk ' { print $1 - $2 } '`
#DROP_NS=`echo "${ORIG_NS} ${HARM_NS}" | awk ' { print $1 - $2 } '`
#DROP_NA=`echo "${ORIG_NA} ${HARM_NA}" | awk ' { print $1 - $2 } '`
#
#echo -e "ID\tORIG_SEC\tHARM_SEC\tDROP_SEC\tORIG_NS\tHARM_NS\tDROP_NS\tORIG_NA\tHARM_NA\tDROP_NA" > ${output}/${id}/nap.nums
#echo -e "${id}\t${ORIG_SEC}\t${HARM_SEC}\t${DROP_SEC}\t${ORIG_NS}\t${HARM_NS}\t${DROP_NS}\t${ORIG_NA}\t${HARM_NA}\t${DROP_NA}" >> ${output}/${id}/nap.nums
#
#
### --------------------------------------------------------------------------------
###
### Stats on annotations; add HARM tag to distinguish from above ANNOTS outputs
###
### --------------------------------------------------------------------------------
#
#
#${NAP_LUNA} ${harmsl} id=${id} ${NAP_LUNA_ARGS} \
#	    -t ${output} ${dom_core} \
#	    stgannots=${NAP_STAGING_ANNOTS} \
#       -s ' TAG HARM/1 & 
#            ANNOTS & 
#            SPANNING annot=${stgannots} ' 2>> $ERR
#
## make ANNOT-specific files for key measures (i.e. so they can be shown easily in Sample Metrics)
## if not present for a individual, show 0 
#
## if no annotation files are attached, then luna_core_ANNOTS_ANNOT_HARM.txt will be missing;  so
## ensure it exists here, even if empty ( i.e. the lines below will correctly indicate 0 annots of each class)
#
#touch ${output}/${id}/luna_core_ANNOTS_ANNOT_HARM.txt
#
#awk -v id=${id}  ' BEGIN {n=0;t=0;m=0 } $2=="apnea" {n=$4; t=$5; m=t/n} END { printf "ID\tCOUNT\tTOT\tMEAN\n"id"\t"n"\t"t"\t"m"\n" } ' \
#    ${output}/${id}/luna_core_ANNOTS_ANNOT_HARM.txt > ${output}/${id}/luna_core_ANNOTS_ANNOT_HARM_apnea.txt 
#
#awk -v id=${id}  ' BEGIN {n=0;t=0;m=0 } $2=="hypopnea" {n=$4; t=$5; m=t/n} END { printf "ID\tCOUNT\tTOT\tMEAN\n"id"\t"n"\t"t"\t"m"\n" } ' \
#    ${output}/${id}/luna_core_ANNOTS_ANNOT_HARM.txt > ${output}/${id}/luna_core_ANNOTS_ANNOT_HARM_hypopnea.txt 
#
#awk -v id=${id}  ' BEGIN {n=0;t=0;m=0 } $2=="arousal" {n=$4; t=$5; m=t/n} END { printf "ID\tCOUNT\tTOT\tMEAN\n"id"\t"n"\t"t"\t"m"\n" } ' \
#    ${output}/${id}/luna_core_ANNOTS_ANNOT_HARM.txt > ${output}/${id}/luna_core_ANNOTS_ANNOT_HARM_arousal.txt 
#
#awk -v id=${id}  ' BEGIN {n=0;t=0;m=0 } $2=="movement" {n=$4; t=$5; m=t/n} END { printf "ID\tCOUNT\tTOT\tMEAN\n"id"\t"n"\t"t"\t"m"\n" } ' \
#    ${output}/${id}/luna_core_ANNOTS_ANNOT_HARM.txt > ${output}/${id}/luna_core_ANNOTS_ANNOT_HARM_movement.txt 
#
#awk -v id=${id}  ' BEGIN {n=0;t=0;m=0 } $2=="artifact" {n=$4; t=$5; m=t/n} END { printf "ID\tCOUNT\tTOT\tMEAN\n"id"\t"n"\t"t"\t"m"\n" } ' \
#    ${output}/${id}/luna_core_ANNOTS_ANNOT_HARM.txt > ${output}/${id}/luna_core_ANNOTS_ANNOT_HARM_artifact.txt 
#
#awk -v id=${id} ' BEGIN {n=0;t=0;m=0} $2=="desat" {n=$4; t=$5; m=t/n} END { printf "ID\tCOUNT\tTOT\tMEAN\n"id"\t"n"\t"t"\t"m"\n" } ' \
#    ${output}/${id}/luna_core_ANNOTS_ANNOT_HARM.txt > ${output}/${id}/luna_core_ANNOTS_ANNOT_HARM_desat.txt 
#
#
#
### --------------------------------------------------------------------------------
###
### Create the base EDF
###
### --------------------------------------------------------------------------------
#
#echo >> $LOG
#echo "Building the base EDF from the harmonized EDF" >> $LOG
#echo "---------------------------------------------" >> $LOG
#
## base sample-list location
#basesl=${output}/${id}/base.lst
#rm -rf ${basesl}
#
#${NAP_LUNA} ${harmsl} id=${id} ${NAP_LUNA_ARGS} \
#	    csfile=${NAP_BASE_OPT_SIGS},${NAP_BASE_SIGS} \
#	    group=${NAP_GROUP} \
#      	    outdir=${output}/${id} \
#	    cmp=${NAP_EDFZ} \
#	    fl=${NAP_BASE_FRQ_LWR} \
#	    fu=${NAP_BASE_FRQ_UPR} \
#	    ftw=${NAP_BASE_FIR_TW} \
#            fripple=${NAP_BASE_FIR_RIPPLE} \
#	    -t ${output} ${dom_core} \
#       -s 'CANONICAL file=${csfile} group=${group} drop-originals &
#           FILTER bandpass=${fl},${fu} tw=${ftw} ripple=${fripple} sig=${eeg},${eog} &
#           WRITE edf-dir=${outdir}/data edfz=${cmp} edf-tag=base with-annots sample-list=${outdir}/base.lst' 2>> $ERR
#
#echo "Written base EDF to : ${output}/${id}/data/${id}-base.edf" >> $LOG
#echo "Written base SL to : ${output}/${id}/base.lst" >> $LOG
#
#
### --------------------------------------------------------------------------------
###
### Report some basic stats on the base EDF (# signals, annots, & duration)
###
### --------------------------------------------------------------------------------
#
#${NAP_LUNA} ${basesl} id=${id} verbose=T -s DESC 2> ${output}/${id}/tmp.headers
#BASE_DUR=`awk ' $1 == "duration:" { print $2 } ' ${output}/${id}/tmp.headers`
#BASE_NS=`awk ' $1 == "signals:" { print $2 } ' ${output}/${id}/tmp.headers`
#BASE_NA=`awk ' $3 == "instance(s)" && $4 == "(from" ' ${output}/${id}/tmp.headers | wc -l | awk ' { print $1 } ' `
#echo "Detected $BASE_NS channels, $BASE_NA annotation classes & duration $BASE_DUR" >> $LOG
#
#
### --------------------------------------------------------------------------------
###
### Figure out which channels are present
###
### --------------------------------------------------------------------------------
#
#trap '' ERR; set +e
#
## does the base EDF contain any EEG signals?
#${NAP_LUNA} ${basesl} id=${id} ${NAP_LUNA_ARGS} silent=T -s 'CONTAINS sig=${eeg}'
#EEG_ABSENT=$?
#
## does the harmonizedbase EDF contain C4_M1?
#${NAP_LUNA} ${harmsl} id=${id} ${NAP_LUNA_ARGS} silent=T SOAP_SIG=${NAP_SOAP_CH} -s 'CONTAINS sig=${SOAP_SIG}'
#SOAP_ABSENT=$?
#
## do the base EDF contain csEEG signal?
#${NAP_LUNA} ${basesl} id=${id} ${NAP_LUNA_ARGS} silent=T -s 'CONTAINS sig=csEEG'
#csEEG_ABSENT=$?
#
## do the base EDF contain csCEN or csFRT signals?
#${NAP_LUNA} ${basesl} id=${id} ${NAP_LUNA_ARGS} silent=T -s 'CONTAINS sig=csCEN,csFRT'
#CEN_FRT_EEG_ABSENT=$?
#
## do the base EDF contain ECG signal?
#${NAP_LUNA} ${basesl} id=${id} ${NAP_LUNA_ARGS} silent=T -s 'CONTAINS sig=csECG'
#ECG_ABSENT=$?
#
## Naming for respiratory signals
##   csCAN : cannula
##   csTHM : thermistor
##   csNAS : nasal pressure
#
## do the base EDF contain respiratory (csCAN) signal?
#${NAP_LUNA} ${basesl} id=${id} ${NAP_LUNA_ARGS} silent=T -s 'CONTAINS sig=csCAN'
#csCAN_ABSENT=$?
#
#set -e; trap "cleanup" ERR
#
#
#if [[ ${CEN_FRT_EEG_ABSENT} -eq 2 ]]; then
#    echo "No central/frontal EEG present: skipping SOAP, PSD, spindle & polarity checks" >> $LOG
#    NAP_DO_POLARITY=0
#    NAP_DO_SPINDLES=0
#    NAP_DO_PSD=0
#fi
#
#if [[ ${SOAP_ABSENT} -ne 0 ]]; then
#    echo "Skipping SOAP (no NAP_SOAP_CH = ${NAP_SOAP_CH} signal in the harmonized EDF)" >> $LOG
#    NAP_DO_SOAP=0
#fi
#
#
#if [[ ${csEEG_ABSENT} -ne 0 ]]; then
#    echo "Skipping POPS (no signals can be mapped to csEEG)" >> $LOG
#    NAP_DO_POPS=0
#fi
#
#
#if [[ ${NAP_DO_RESPIRATORY} -eq 1 ]]; then
#    if [[ ${csCAN_ABSENT} -ne 0 ]]; then
#	echo "Skipping respiratory analysis (no signal can be mapped to csCAN)" >> $LOG
#	NAP_DO_RESPIRATORY=0
#    fi
#fi
#
#
#
### --------------------------------------------------------------------------------
###
### Summary (viewing) SIGSTATS/STATS for the harmonized EDF (i.e will have different channels
### etc from the original, potentiallu)
##y
### --------------------------------------------------------------------------------
#
#echo "Running harmonized EDF SIGSTATS..." >> $LOG
#${NAP_LUNA} ${harmsl} ${NAP_LUNA_ARGS} -t ${output} ${dom_stats} tt-append=_HARM-1 -s 'SIGSTATS epoch sr-over=50' 2>> $ERR
#
#echo "Running harmonized EDF STATS..." >> $LOG
#${NAP_LUNA} ${harmsl} ${NAP_LUNA_ARGS} -t ${output} ${dom_stats} tt-append=_HARM-1 -s 'STATS epoch sr-under=50' 2>> $ERR
#
#
#
### --------------------------------------------------------------------------------
###
### All subsequent analyses are on the base EDF (i.e. derived metrics and issues for NAP)
### The base EDF should have consistently labelled channels, annotations, referencing, etc
### It will have EDF record size of 1 and aligned staging
###
### --------------------------------------------------------------------------------
#
#
#echo >> $LOG
#echo "Sleep staging evaluation/prediction" >> $LOG
#echo "-----------------------------------" >> $LOG
#
#
### --------------------------------------------------------------------------------
###
### SOAP : requries C4_M1 by default
###
### --------------------------------------------------------------------------------
#
#
#if [[ ${NAP_DO_SOAP} -eq 1 ]]; then
#    
#    if [[ ${STAGES_ABSENT} -eq 1 ]]; then
#	echo "Skipping SOAP (no existing stage information)..." >> $LOG
#    else	
#        echo "Running SOAP on ${NAP_SOAP_CH}..." >> $LOG
#	
#	${NAP_LUNA} ${harmsl} id=${id} ${NAP_LUNA_ARGS} \
#	  	    adir=${output}/${id}/annots/ \
#		    model=${NAP_SOAP_MOD} \
#                    stgannnots=${NAP_CORE_STAGING_ANNOTS} \
#	            alias="SOAP_SIG|${NAP_SOAP_CH}" \
#		    -t ${output} ${dom_suds} \
#       	            -s 'EPOCH align=${stgannnots} &
#	                SOAP model=${model} th=8 epoch robust=0.01 annot=SOAP annot-dir=${adir}' 2>> $ERR
#  
#	# hack to unf*ck row-order formatting issue w/ -t option for some Luna commands
#	# not needed for 'E'
#		    
#	${NAP_FIXROWS} ID   < ${output}/${id}/luna_suds_SOAP.txt > ${output}/${id}/fixed_SOAP.txt
#	mv ${output}/${id}/fixed_SOAP.txt ${output}/${id}/luna_suds_SOAP.txt
#	
#	${NAP_FIXROWS} ID SS < ${output}/${id}/luna_suds_SOAP_SS.txt > ${output}/${id}/fixed_SOAP_SS.txt
#	mv ${output}/${id}/fixed_SOAP_SS.txt ${output}/${id}/luna_suds_SOAP_SS.txt
#    fi
#    
#fi
#
#
#
### --------------------------------------------------------------------------------
###
### POPS, expecting a canonical csEEG channel, param fixed to match NSRR training data 
###
### --------------------------------------------------------------------------------
#
#if [[ -d "${NAP_POPS_DIR}" && ${NAP_DO_POPS} -eq 1 ]] 
#then
#    echo "Running (single-EEG) POPS, using training data in ${NAP_POPS_DIR}" >> $LOG
#
#    # note, if PRIORS == "." then Luna skips it (i.e. so default value is set to that)
#    if [ "${NAP_POPS_PRIORS}" != "." ] ; then
#        echo "Applying elapsed-sleep priors: ${NAP_POPS_PRIORS}" >> $LOG
#    fi
#
#    # add POPS esprior, annot output
#    
#    ${NAP_LUNA} ${basesl} ${NAP_LUNA_ARGS} \
#		apath=${output}/${id}/annots \
#		ftr=${NAP_POPS_FTR} \
#		mod=${NAP_POPS_MOD} \
#                pth=${NAP_POPS_PATH} \
#		esp=${NAP_POPS_PRIORS} \
#		-t ${output} ${dom_suds} \
#                -s 'POPS features=${ftr} model=${mod} es-priors=${esp} path=${pth}' > /dev/null 2>> $ERR
#
#else
#    echo "Skipping POPS automated staging"
#fi
#
#
#
### --------------------------------------------------------------------------------
###
### Flag likely artifacts --> ${ISSUES}
###
### --------------------------------------------------------------------------------
#
#
### --------------------------------------------------------------------------------
###
### Epoch/channel level masks
###
### --------------------------------------------------------------------------------
#
#echo >> $LOG
#echo "Flagging likely artifact epochs/channels" >> $LOG
#echo "----------------------------------------" >> $LOG
#
#${NAP_LUNA} ${basesl} id=${id} ${NAP_LUNA_ARGS} \
#	    epth=${NAP_CHEP_TH} \
#       -t ${output} ${dom_core} \
#       -s ' EPOCH & 
#            CHEP-MASK sig=${eeg},${eog},${emg} flat=0.1 clipped=0.1 ep-th=${epth} & 
#            CHEP dump sig=${eeg},${eog},${emg} & 
#            CHEP epochs &
#            MASK regional=3,5 & 
#            DUMP-MASK ' 2>> $ERR
#
## for sample metrics
#awk -v id=${id} ' BEGIN { bad=0; good=0 } { ++good } $3 == 1 { ++bad } END { printf "ID\tBAD\tPROP\n"id"\t"bad"\t"bad/good"\n" } ' OFS="\t" \
#    ${output}/${id}/luna_core_DUMP-MASK_E.txt > ${output}/${id}/bad.base.epochs.summary.txt
#
#
##  
## absence of staging data / poor coverage, i.e. truncated staging
## gross artifact (EEG only, by epoch and channel)
## flag potential EEG polarity issues (negative / ambiguous channels)
## flag likely ECG contamination (e.g. EEG/ECH coherence)
## flag likely line noise / EDFs with peaks in the <35 Hz range
## flag if likely units issue in the EEG (i.e. large values, or drift)
## SOAP : flag if low SOAP consistency
#
### --------------------------------------------------------------------------------
###
### SOAP: staging/EEG inconsistencies
###
### --------------------------------------------------------------------------------
#
#if [[ -f "${output}/${id}/luna_suds_SOAP.txt" ]]; then
#    
#    k3=`cat ${output}/${id}/luna_suds_SOAP.txt | ${NAP_BEHEAD} | awk ' $1 == "K3" && $2 < 0.50 { print $2 } '`
#
#    if [[ ! -z ${k3} ]]; then
#	echo "Flagging low (<0.5) SOAP kappa (3-class) - possible staging/EEG inconsistencies" >> ${LOG} 
#	echo ${k3} | awk ' { print "SOAP" , "." , "K3="$1 } ' OFS="\t" > ${ISSUES}
#    fi
#    
#fi
#
#
#
### --------------------------------------------------------------------------------
###
### ECG-CONTAM: possible EEG/ECG contamination
###
### --------------------------------------------------------------------------------
#
#if [[ ${EEG_ABSENT} -ne 0 ]]; then
#    if [[ ${ECG_ABSENT} -ne 0 ]]; then
#	echo "Scanning for possible EEG cardiac contamination" >> ${LOG}
#	
#    fi    
#fi
#
#
### --------------------------------------------------------------------------------
###
### Flag likely line noise / peaks in EEG sleep PSD
###
### --------------------------------------------------------------------------------
#
#
#if [[ ${NAP_DO_PSD} -eq 1 ]]; then
#
#    echo "Evaluating peak/line-noise artifacts in sleep EEG PSD (<35 Hz)..." >> $LOG
#
#    ${NAP_LUNA} ${basesl} ${NAP_LUNA_ARGS} -t ${output} ${dom_spec} \
#                -s 'MASK all & MASK unmask-if=N1,N2,N3,R & RE & 
#                    TAG PKS/35 & PSD sig=${eeg} dB max=35 spectrum peaks' 2>> $ERR
#
#fi
#
#
#
#
#
### --------------------------------------------------------------------------------
###
### Flag likely EEG polarity issues? 
###
### --------------------------------------------------------------------------------
#
#if [[ ${NAP_DO_POLARITY} -eq 1 ]]; then
#
#    # initialize these (so they are not empty downstream, if the POL step is skipped)
#    echo "__dummy__" > ${output}/${id}/neg.chs
#    echo "__dummy__" > ${output}/${id}/ambig.chs
#    
#    if [[ ${STAGES_ABSENT} -eq 0 ]]; then
#	echo "Diagnosing likely EEG polarity flips in the base EDF..." >> $LOG
#	
#	# if manual staging was present, use that to find NREM segments
#	# otherwise, use the POPS-determined staging
#	
#	# write to a temporary, use pol.db rather than a text-table
#	# run this for all EEG channels
#	
#	${NAP_LUNA} ${basesl} id=${id} ${NAP_LUNA_ARGS} \
#		    -o ${output}/${id}/pol.db \
#  		    -s 'MASK all & MASK unmask-if=N2,N3 & RE &
#		    FILTER sig=csCEN,csFRT bandpass=0.3,18 tw=1 ripple=0.02 &
#                    CHEP-MASK sig=csCEN,csFRT ep-th=3 & 
#		    CHEP sig=csCEN,csFRT epochs & RE & 
#                    SPINDLES sig=csCEN,csFRT fc=15 so mag=2 all-spindles ignore-neg-peak &
#                    POL sig=csCEN,csFRT & 
#		    STATS sig=csCEN,csFRT '
#
#	# SPINDLES: expect COUPL_ANGLE 180 < X < 360 
#	# SO:       expect longer SO peak : output if NEG_DUR POS_DUR in cols 3 and 4
#	# POL       expect a **negative** (urgh...) T_DIFF statistic
#	# STATS:    expect a negative SKEW of the EEG (e.g. >1) 
#	
#	${NAP_DESTRAT} ${output}/${id}/pol.db +SPINDLES -r CH F -v COUPL_ANGLE \
#	    | awk ' NR != 1 && $4 < 180 { print $2 } ' > ${output}/${id}/neg.chs
#	
#	${NAP_DESTRAT} ${output}/${id}/pol.db +SPINDLES -r CH -v SO_NEG_DUR SO_POS_DUR \
#	    | awk ' NR != 1 && $3 > $4 { print $2 } ' >> ${output}/${id}/neg.chs
#	
#	${NAP_DESTRAT} ${output}/${id}/pol.db +POL -r CH -v T_DIFF \
#	    | awk ' NR != 1 && $3 > 0 { print $2 } ' >> ${output}/${id}/neg.chs
#	
#	${NAP_DESTRAT} ${output}/${id}/pol.db +STATS -r CH -v SKEW \
#	    | awk ' NR != 1 && $3 > 0 { print $2 } ' >> ${output}/${id}/neg.chs
#	
#	# repeat, for output
#	echo | awk ' { print "ID" , "CH" , "METHOD" , "POL" } ' OFS="\t" \
#		   > ${output}/${id}/luna_core_FLIP_CH_METHOD.txt
#    
#	${NAP_DESTRAT} ${output}/${id}/pol.db +SPINDLES -r CH F -v COUPL_ANGLE | \
#	    awk ' NR != 1 { print $1,$2,"COUPL" , ($4 < 180 ) ? "-ve" : "+ve"  } ' OFS="\t" \
#		>> ${output}/${id}/luna_core_FLIP_CH_METHOD.txt
#	
#	${NAP_DESTRAT} ${output}/${id}/pol.db +SPINDLES -r CH -v SO_NEG_DUR SO_POS_DUR  | \
#	    awk ' NR != 1 { print $1,$2,"SO", ( $3>$4 ) ? "-ve" : "+ve"  } ' OFS="\t" \
#		>> ${output}/${id}/luna_core_FLIP_CH_METHOD.txt
#	
#	${NAP_DESTRAT} ${output}/${id}/pol.db +POL -r CH -v T_DIFF | \
#	    awk ' NR != 1 { print $1,$2,"POL", ($3>0) ? "-ve" : "+ve"  } ' OFS="\t" \
#		>> ${output}/${id}/luna_core_FLIP_CH_METHOD.txt
#	
#	${NAP_DESTRAT} ${output}/${id}/pol.db +STATS -r CH -v SKEW | \
#	    awk ' NR != 1 { print $1,$2,"SKEW", ($3>0) ? "-ve" : "+ve"  } ' OFS="\t" \
#		>> ${output}/${id}/luna_core_FLIP_CH_METHOD.txt
#	
#	# get command-delimited list of channels to flip, or that are ambiguous
#	sort ${output}/${id}/neg.chs | uniq -c | awk ' $1 > 2 { print $2 } ' | \
#	    paste -s -d ',' - | awk ' NF != 0 ' > ${output}/${id}/neg.chs.tmp
#	
#	sort ${output}/${id}/neg.chs | uniq -c | awk ' $1 == 2 { print $2 } ' | \
#	    paste -s -d ',' - | awk ' NF != 0 ' > ${output}/${id}/ambig.chs
#	
#	mv ${output}/${id}/neg.chs.tmp ${output}/${id}/neg.chs 
#	
#	# if nothing to flip, set to '__dummy__'  (i.e. a string we do not expect to match a channel name)
#	flip=`wc -l ${output}/${id}/neg.chs | awk ' { print $1 }' `
#	if [[ $flip -eq 0 ]]; then
#	    echo "__dummy__" > ${output}/${id}/neg.chs
#	fi
#
#	flip=`wc -l ${output}/${id}/ambig.chs | awk ' { print $1 }' `
#	if [[ $flip -eq 0 ]]; then
#	    echo "__dummy__" > ${output}/${id}/ambig.chs
#	fi
#	
#    fi
#fi
#
#
#
#
#
### --------------------------------------------------------------------------------
###
### Spectral analyses
###
### --------------------------------------------------------------------------------
#
#
#echo >> $LOG
#echo "Sleeep micro-architecture (spectral/spindle) analyses" >> $LOG
#echo "-----------------------------------------------------" >> $LOG
#
#
#if [[ ${NAP_DO_PSD} -eq 1 ]]; then
#
#    echo "Running N2/N3 Welch PSD..." >> $LOG
#
#    ${NAP_LUNA} ${basesl} ${NAP_LUNA_ARGS} -t ${output} ${dom_spec} tt-append=_SS-N2 \
#                -s 'MASK ifnot=N2 & RE & PSD sig=csCEN,csFRT dB max=25 spectrum' 2>> $ERR
#
#    ${NAP_LUNA} ${basesl} ${NAP_LUNA_ARGS} -t ${output} ${dom_spec} tt-append=_SS-N3 \
#                -s 'MASK ifnot=N3 & RE & PSD sig=csCEN,csFRT dB max=25 spectrum' 2>> $ERR
#    
#    
#fi
#
#
#
### --------------------------------------------------------------------------------
###
### MTM spectrograms (base EDF, single EEG)  -- SKIP FOR NOW --
###
### --------------------------------------------------------------------------------
#
##if [ ${NAP_DO_MTM} -eq 1 ]; then
##    
##    echo "Running MTM spectrogram band summaries (single canonical EEG) " >> $LOG
##    
##    # create a single csEEG summary epoch-level
##    ${NAP_LUNA} ${basesl} ${NAP_LUNA_ARGS} \
##                -t ${output} ${dom_spec} tt-append=_SPSD-1 \
##                -s 'MTM tw=5 segment-sec=5 segment-inc=1 epoch max=45 sig=csEEG' 2>> $ERR
##
##fi
#
#
#
### --------------------------------------------------------------------------------
###
### Spindle/SO analyses
###
### --------------------------------------------------------------------------------
#
#
#if [[ ${NAP_DO_SPINDLES} -eq 1 ]]; then
#
#    echo "Running SPINDLES (N2)..." >> $LOG
#    
#    ${NAP_LUNA} ${basesl} id=${id} ${NAP_LUNA_ARGS} \
#		-t ${output} ${dom_spso} tt-append=_SS-N2 \
#		adir=${output}/${id}/annots/ \
#		flipchs=`cat ${output}/${id}/neg.chs` \
#		-s 'MASK ifnot=N2 & RE & 
#                    FLIP sig=${flipchs} & 
#                    CHEP-MASK sig=csCEN,csFRT ep-th=3,3 &
#                    CHEP sig=csCEN,csFRT epochs & RE &
#                    SPINDLES sig=csCEN,csFRT fc=11,15 so mag=2 nreps=1000 annot=sp-N2 annot-dir=${adir}' 2>> $ERR
#
#    echo "Running SPINDLES (N2+N3)..." >> $LOG
#
#    ${NAP_LUNA} ${basesl} id=${id} ${NAP_LUNA_ARGS} \
#		adir=${output}/${id}/annots/ \
#		flipchs=`cat ${output}/${id}/neg.chs` \
#		-t ${output} ${dom_spso} tt-append=_SS-N23 \
#		-s 'MASK all & MASK unmask-if=N2,N3 & RE & 
#                    FLIP sig=${flipchs} & 
#                    CHEP-MASK sig=csCEN,csFRT ep-th=3,3 & RE & 
#  		    CHEP sig=csCEN,csFRT epochs & RE &
#                    SPINDLES sig=csCEN,csFRT fc=11,15 so mag=2 nreps=1000 annot=sp-N23 annot-dir=${adir}' 2>> $ERR
#
#fi
#
#
#
### --------------------------------------------------------------------------------
###
### Respiratory signal analysis
###
### --------------------------------------------------------------------------------
#
#if [[ ${NAP_DO_RESPIRATORY} -eq 1 ]]; then
#    
#    echo >> $LOG
#    echo "Respiratory signal analyses" >> $LOG
#    echo "---------------------------" >> $LOG
#    
#    resp_output=${output}
#    edfname=`awk -v id=${id} -F"\t" ' $1 == id { print $2 } ' ${harmsl}`
#    xmlname=`awk -v id=${id} -F"\t" ' $1 == id { print $3 } ' ${harmsl}`
#    
#    # if using .edfz, then make a temporary .edf
#    if [[ ${NAP_EDFZ} -eq 1 ]]; then
#	echo "Creating a temporary uncompressed EDF" >> $LOG
#	cp ${edfname} ${output}/${id}/data/${id}_resp.edf.gz
#	gunzip ${output}/${id}/data/${id}_resp.edf.gz
#	edfname="${output}/${id}/data/${id}_resp.edf"
#    fi
#
#    # Converting relative path to absolute path
#    if [[ ! "${input}" = /* ]]; then
#        pwd_dir=$(pwd)
#        resp_output=${pwd_dir}/${output}
#        edfname=${pwd_dir}/${edfname}
#        xmlname=${pwd_dir}/${xmlname}
#    fi
#    
#    echo "Writing results to : ${resp_output}" >> $LOG
#    #     echo "id : ${id} " , edfname: ${edfname}, xmlname: ${xmlname}" >> $LOG
#    
#    # endotype calculations (Matlab)
#    
#    ${NAP_MATLAB} -nodisplay \
#		  -r "StartHere ${resp_output} ${id}/ ${edfname} ${xmlname}" \
#		  -sd ${NAP_DIR}"/respiratory_endotypes" \
#		  -logfile ${resp_output}/${id}/outputconvert.log 2>> $ERR
#
#    # clean up temporary *_resp.edf
#    if [[ ${NAP_EDFZ} -eq 1 ]];then
#	rm -f ${edfname}		
#    fi
#
#fi
#
#
#
### --------------------------------------------------------------------------------
###
### ECG signal analysis
###
### --------------------------------------------------------------------------------
#
## <--- todo --->
#
#
#
### --------------------------------------------------------------------------------
###
### Compile key results into R dataframes for viewing (luna/shiny)
###
### --------------------------------------------------------------------------------
#
#echo >> $LOG
#echo "Coda" >> $LOG
#echo "----" >> $LOG
# 
#echo "Compiling tables into RData files..." >> $LOG
# 
## i.e.     path/to/folder/text.txt
## becomes  path/to/folder/text.txt-tab.RData
##      or  path/to/folder/text.txt-fig.RData
## luna-shiny then automatically loads any *-tab.RData and *-fig.RData files
#
## add fextract() calls into coda2.R to create particular tables 
## also see example for PSD plots for how to save / attach images
##  these point to .png files, which can either be created in coda2.R based
##  on summary stats, or indpendently (in which case a *-fig.RData file is created
##  which just points to the existing .png
#
#${NAP_R} ${NAP_DIR}/coda2.R ${NAP_DIR} ${NAP_CODA_RESOURCE_DIR} ${output}/${id} >> $ERR 2>&1
# 
#
### --------------------------------------------------------------------------------
###
### Tidy up?
###
### --------------------------------------------------------------------------------
#
#if [ ${NAP_TIDY} -eq 1 ]; then
#
#    echo "Tidying up intermediates..." >> $LOG
#    
#    rm -f ${output}/${id}/*.txt
#    rm -f ${output}/${id}/*.txt.gz
#    rm -f ${output}/${id}/*.db
#    rm -f ${output}/${id}/tmp.headers
#    rm -f ${output}/${id}/nap.nums
#    rm -f ${output}/${id}/neg.chs
#    rm -f ${output}/${id}/ambig.chs
#    
#    rm -f ${output}/${id}/rs1.lst
#    rm -f ${TMP_EDF}
#    rm -f ${output}/${id}/annots/a.annot
#
#    rm -f *.txt
#    rm -f *.txt.gz
#    rm -f *.db
#    rm -f tmp.headers
#    rm -f nap.nums
#    rm -f rs1.lst
#    rm -f neg.chs
#    rm -f ambig.chs
#
#fi
#
#
### --------------------------------------------------------------------------------
###
### All done
###
### --------------------------------------------------------------------------------
#
#echo "All done." >> $LOG
#echo >> $LOG
#
#dt=$(date '+%d/%m/%Y %H:%M:%S');
# 
#echo "--------------------------------------------------------------------------------" >> $LOG
#echo "NAP ${NAP_VERSION} | ${id} | process completed: ${dt} "                         >> $LOG
#echo "--------------------------------------------------------------------------------" >> $LOG
#
## leave a marker: 1 = success
#echo -e "${id}\t1" > ${output}/${id}/nap.status
# 
#exit 0

