for f in `ls /data/nsrr/datasets/nchsdb/sleep_data/*.tsv | xargs -n 1 basename`
do
    echo "$f"
    fannot=`echo $f | sed 's/\.tsv/\.annot/g'`
    awk -F"\t" ' NR != 1 { print $3 , $1 , "+"$2 } ' OFS="\t" /data/nsrr/datasets/nchsdb/sleep_data/${f} | tr ':' '_' > ${fannot}
done
