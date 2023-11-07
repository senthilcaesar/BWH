c="/data/nsrr/datasets/stages/original/STAGES PSGs/"
d="MAYO"

while read -r id
do
    inputFile="${c}/${d}/${id}.csv"

    if test -f "$inputFile"; then

        tr -d ' ' < $inputFile | awk -F"," ' $3=="Wake"   { print "W" , "." , "." , $1 , "+30" , "." } 
                                             $3=="Stage1" { print "N1" , "." , "." , $1 , "+30" , "." } 
                                             $3=="Stage2" { print "N2" , "." , "." , $1 , "+30" , "." } 
                                             $3=="Stage3" { print "N3" , "." , "." , $1 , "+30" , "." } 
                                             $3=="REM"    { print "R" , "." , "." , $1 , "+30" , "." } '  \
                                             OFS="\t" > /data/purcell/projects/saps/annots/stages/${d}/${id}.annot 
                                             #echo /data/purcell/projects/saps/annots/stages/${d}/${id}.annot
    else
        echo "$inputFile does not exist!"
    fi

done < ${d}.txt
