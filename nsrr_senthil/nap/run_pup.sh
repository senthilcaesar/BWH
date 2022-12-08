

id="resp"
edfname="/Users/sq566/Desktop/Respiratory/resp.edf"
xmlname="/Users/sq566/Desktop/Respiratory/resp.annot"
resp_output="/Users/sq566/Desktop/nap"
NAP_DIR='/Users/sq566/Downloads/nsrr/nap'


/Applications/MATLAB_R2022b.app/bin/matlab -nodisplay \
       -r "StartHere ${resp_output} ${id}/ ${edfname} ${xmlname}" \
       -sd ${NAP_DIR}"/respiratory_endotypes" \
       -logfile ${resp_output}/${id}/outputconvert.log 2>> resp.log

