# Running POPS on a single subject 
# s.lst
# learn-nsrr02	./learn-nsrr02.edf
luna s.lst 1 -o pops.db -s ' FILTER sig=EEG bandpass=0.3,35 tw=0.5 ripple=0.02
                             COPY sig=EEG tag=_NORM
                             ROBUST-NORM sig=EEG_NORM epoch winsor=0.005 second-norm=T
                             POPS alias=CEN,ZEN|EEG,EEG_NORM path=pops lib=s2 
                             & WRITE-ANNOTS annot=N1,N2,N3,R,W file=annots/^-pops.annot hms '
                
 
# Running POPS in cluster for many subjects
/data/nsrr/bin/runner2.sh 30 stages_pops.lst param/main cmd/pops.txt o tmp/pops sig=C3

cat cmd/pops.txt
FILTER sig=C3 bandpass=0.3,35 tw=0.5 ripple=0.02 
COPY sig=C3 tag=_NORM 
ROBUST-NORM sig=C3_NORM epoch winsor=0.005 second-norm=T 
POPS alias=CEN,ZEN|C3,C3_NORM path=/PHShome/sq566/nsrr/common/resources/pops lib=s2 
WRITE-ANNOTS annot=N1,N2,N3,R,W file=/data/purcell/projects/saps/annots/stages/POPS/^-pops.annot hms


# For two channel C3 and C4, Please use the below
FILTER sig=C3,C4 bandpass=0.3,35 tw=0.5 ripple=0.02
COPY sig=C3,C4 tag=_NORM
ROBUST-NORM sig=C3_NORM,C4_NORM epoch winsor=0.005 second-norm=T
POPS alias=CEN,ZEN|C3,C3_NORM equiv=CEN,ZEN|C4,C4_NORM path=/PHShome/sq566/nsrr/common/resources/pops lib=s2
WRITE-ANNOTS annot=N1,N2,N3,R,W file=/data/purcell/projects/saps/annots/stages/POPS/^-pops.annot hms
