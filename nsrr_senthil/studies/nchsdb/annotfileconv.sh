#!/bin/bash
#################################################################
#
# NCHSDB annotation format - converting to NSRR
#
#################################################################

DIR=/data/nsrr/working/nchsdb/Sleep_Data/

for f in ${DIR}/*.tsv
do

# get ID
id=`basename $f | sed 's/\.tsv//g'`
echo "processing $id"

awk -F"\t" ' BEGIN { printf "# misc | Misc annots | note[txt]\n"; } \
             FNR == NR { terms[$2] = $1; next }  
             $1 == "onset" { next } 
             FNR != NR { ac = $3 ; fnd = 0 ;
                         for (t in terms) { if ( ac == t ) { ac = terms[t] ; fnd=1; break } } ;
                         if ( fnd == 1 ) print ac , "." , "." , $1 , "+"$2 , "." ;
			 else { 
                           gsub( " ", "_" , ac ) ; 
                           print "misc" , "." , "." , $1 , "+"$2 , ac 
                        } } ' OFS="\t" ${DIR}/final/annot.mapping $f > ${DIR}/annots/${id}.annot
done
