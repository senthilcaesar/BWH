#!/bin/bash

# --------------------------------------------------------------------------------
#
# ERIS-specific configuration code
#
# --------------------------------------------------------------------------------

# ERIS-required version of Matlab
ERIS_MATLAB_VERSION="matlab/2019b"

# check whether this is available via LSF 'module load'
module_exists=`module avail ${ERIS_MATLAB_VERSION} 2>&1 | wc -l`

if [[ ! ${module_exists} -eq 0 ]]; then
    module load ${ERIS_MATLAB_VERSION}
else
    echo "cannot find ${ERIS_MATLAB_VERSION}... exiting"
    exit 1
fi
