luna --build /data/purcell/projects/rasp/xtract/bch > sl/bch.lst
luna --build /data/purcell/projects/rasp/xtract/geis /data/purcell/projects/rasp/pops/annots/geis > sl/geis.lst
luna --build /data/purcell/projects/rasp/xtract/nimh > sl/nimh.lst
luna --build /data/purcell/projects/rasp/xtract/nyu > sl/nyu.lst
luna --build /data/purcell/projects/rasp/xtract/tch /data/purcell/projects/rasp/pops/annots/tch > sl/tch.lst

/data/nsrr/bin/dev-runner2.sh 10 sl/bch.lst . build.txt o tmp/bch "d=bch"
/data/nsrr/bin/dev-runner2.sh 10 sl/geis.lst . build.txt o tmp/geis "d=geis"
/data/nsrr/bin/dev-runner2.sh 20 sl/nimh.lst . build.txt o tmp/nimh "d=nimh"
/data/nsrr/bin/dev-runner2.sh 5 sl/nyu.lst . build.txt o tmp/nyu "d=nyu"
/data/nsrr/bin/dev-runner2.sh 30 sl/tch.lst . build.txt o tmp/tch "d=tch"

cp -avf /data/purcell/projects/rasp/xtract/bch/*.annot bch/
cp -avf /data/purcell/projects/rasp/xtract/nimh/*.annot nimh/
cp -avf /data/purcell/projects/rasp/xtract/nyu/*.annot nyu/
cp -avf /data/purcell/projects/rasp/pops/annots/geis/* geis
cp -avf /data/purcell/projects/rasp/pops/annots/tch/* tch

# when EDFZ jobs finished...
luna --build bch sl/bch.zlst
luna --build geis sl/geis.zlst
luna --build nimh sl/nimh.zlst
luna --build nyu sl/nyu.zlst
luna --build tch sl/tch.zlst

awk ' { printf "bch\t"$1" \tbch/"$1".edf.gz\n" ; printf "bch\t"$1"\tbch/"$1".annot\n" } ' OFS="\t" sl/bch.lst | awk ' { gsub( /-a/ , "" , $2) ; print "bch" , $2 , $3} ' OFS="\t" > manifests/rasp.txt
awk ' { printf "geis\t"$1"\tgeis/"$1".edf.gz\n" ; printf "geis\t"$1"\tgeis/"$1".annot\n" } ' OFS="\t" sl/geis.lst >> manifests/rasp.txt
awk ' { printf "nimh\t"$1"\tnimh/"$1".edf.gz\n" ; printf "nimh\t"$1"\tnimh/"$1".annot\n" } ' OFS="\t" sl/nimh.lst >> manifests/rasp.txt
awk ' { printf "nyu\t"$1"\tnyu/"$1".edf.gz\n" ; printf "nyu\t"$1"\tnyu/"$1".annot\n" } ' OFS="\t" sl/nyu.lst >> manifests/rasp.txt
awk ' { printf "tch\t"$1"\ttch/"$1".edf.gz\n" ; printf "tch\t"$1"\ttch/"$1".annot\n" } ' OFS="\t" sl/tch.lst >> manifests/rasp.txt

cat manifests/* > manifest.txt

# copy manifest from ERIS -> hypnos
cp ~/mount/manifest.txt /home/senthil/src/moonbeam/dat


# Update db and Phenotype
./moonbeam.cgi -add backup.db -t dat/permissions.txt -c dat/cohorts.txt -f dat/manifest.txt
./add_phenotype.sh

cp backup.db ~/moonbeam/mb.db
