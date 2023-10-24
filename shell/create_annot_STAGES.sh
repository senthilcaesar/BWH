c="/data/nsrr/datasets/stages/original/STAGES PSGs/"
d="STNF"

while read -r id
do
 # remove spaces (comma-delim), make .annot 
 tr -d ' ' < "${c}/${d}/${id}.csv" | awk -F"," ' $3=="Wake" || $3 == "Awake" { print "W" , "." , "." , $1 , "+30" , "." } 
                                        $3=="Stage1" { print "N1" , "." , "." , $1 , "+30" , "." } 
                                        $3=="Stage2" { print "N2" , "." , "." , $1 , "+30" , "." } 
                                        $3=="Stage3" { print "N3" , "." , "." , $1 , "+30" , "." } 
                                        $3=="REM" { print "R" , "." , "." , $1 , "+30" , "." } '  \
                                        OFS="\t" > /data/purcell/projects/saps/annots/stages/${d}/${id}.annot 
                                        echo /data/purcell/projects/saps/annots/stages/${d}/${id}.annot
done < ${d}.txt
