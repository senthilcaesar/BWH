# for test/development purposes, define vars to be able to re-run specific
# components of nap1.sh

id=19282_58
input=.
run=proj1

STAGES_ABSENT=0
#USE_HARM_ANNOTS

# automatic/defaults
output=${input}/nap/
harmsl=${output}/${id}/harm.lst
basesl=${output}/${id}/base.lst
NAP_LUNA_ARGS="--log"
LOG=${input}/nap/${id}/nap.log
ERR=${input}/nap/${id}/nap.err
ISSUES=${output}/${id}/nap.issues
NAP_SUDS_DIR=
NAP_SUDS_ESPRIORS=.
NAP_MTM_MIN_SAMPLE_RATE=50
dom_core="tt-prepend=luna_core_"
dom_stats="tt-prepend=luna_stats_"
dom_macro="tt-prepend=luna_macro_"
dom_suds="tt-prepend=luna_suds_"
dom_spec="tt-prepend=luna_spec_"
dom_spso="tt-prepend=luna_spso_"
NAP_STAGING_ANNOTS="N1,N2,N3,R,W,?,U,M,L"
NAP_DO_CODA1=1
NAP_LUNA=luna
NAP_R=Rscript
NAP_DIR=nsrr/nap/
NAP_CODA_RESOURCE_DIR=.
NAP_BASE_CANONICAL_SIGS=nsrr/common/resources/base.canonical.sigs
NAP_CANONICAL_GROUP=.


