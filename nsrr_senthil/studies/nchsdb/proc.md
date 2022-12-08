# NCH-SDB cohort

- NSRR harmonization, v1

- Created 22-Mar-2021


## Original data

All original EDFs and annotations (`.tsv`) files are located here:

```
DATA=/data/nsrr/working/nchsdb/Sleep_Data/
```

Filename/ID format: two fields, underscore separated:
```
   STUDY_PAT_ID   SLEEP_STUDY_ID 
```

Expecting n=3,984 studies/annotations

```
ls ${DATA}/*edf | wc -l
```
```
 3984
```

```
ls ${DATA}/*tsv | wc -l
```
```
3984
```


## Local working directories

Temporary scratch space 
```
mkdir tmp 
```

Harmonzied annotations
```
mkdir annots
```

Harmonzied/canonical EDFs
```
mkdir edfs
```

Final EDFs and annotations, i.e. to be postde to NSRR
```
mkdir final
```

## Checking IDs

### IDs from EDF filenames

First, get a list of IDs from the EDF filenames:

```
ls ${DATA} | grep "\.edf" | sed 's/\.edf//g' > tmp/ids
```

```
head tmp/ids 
```

```
10000_17728
10003_26038
10006_20647
...
```

Then confirm the _assumption_ that all IDs have correct form, i.e. 
two underscore-delimited fields:

```
awk -F"_" ' { print NF } ' tmp/ids | sort | uniq -c
```
```
3984 2
```

### Number of unique individuals:

```
awk -F"_" ' { print $1 } ' tmp/ids | sort | uniq | wc -l 
```
```
3673
```

### Number of PSGs per individual:

```
awk -F"_" ' { print $1 } ' tmp/ids | sort | uniq -c | awk ' { print $1 } ' | sort | uniq -c
```
```
   3400 1
    238 2
     33 3
      1 4
      1 5
```

As sanity check: 3400 + 238 + 33 + 1 + 1 = 3673 

### Check that EDF and .tsv IDs line up exactly:

```
ls ${DATA} | grep "\.tsv" | sed 's/\.tsv//g' > tmp/ids.tsv
```
```
paste tmp/ids tmp/ids.tsv | awk ' $1 != $2 ' 
```

(no output: i.e. good, this implies that all EDF and `.tsv` filenames match exactly)


## Summarize annotations

The format of `.tsv` files: one header row, three tab-delimited columns:

  - onset (seconds elapsed from EDF start)
  - duration 
  - annotation label

### One example `.tsv`:

```
head ${DATA}/10000_17728.tsv 
```
```
 onset        duration   description
 1.03515625   0.0        Montage:NCH_PSG_STANDARD, Ref
 1.03515625   0.0        Start Recording
 1.03515625   0.0        Montage:NCH_PSG_STANDARD, Ref
 1.1640625    0.0        Recording Analyzer - ECG
 1.6171875    0.0        Video Recording ON
 56.78125     0.0        Montage:Channel Test Referential
 56.78125     0.0        Montage:Channel Test Referential
 117.109375   0.0        Montage:NCH_PSG_STANDARD
 117.109375   0.0        Montage:NCH_PSG_STANDARD
 ...
```


As expected, stage annotations are 30-second epochs; also contains markers (0.0 duration) and arbitrary events (e.g. 7-sceond desat)

```
24780.9609375   30.0    Sleep stage R
24804.70703125  7.609375        EEG arousal
24810.9609375   30.0    Sleep stage R
24812.62890625  0.0     Body Position: Right
24827.36328125  7.0     Oxygen Desaturation
24839.36328125  12.0    Oxygen Desaturation
24840.9609375   30.0    Sleep stage R
24870.9609375   30.0    Sleep stage R
```

__Note:__ Many epoch annotations are not aligned with EDF records (i.e. for this individual
starts 0.9609375 seconds past the start of each integer second;  this will need to be 
addressed when working with an editted EDF)


### Dump all annotations

List most commonly-occurring first:

```
cat ${DATA}/*.tsv | awk -F"\t" ' { print $3 } ' | sort | uniq -c | sort -nr > tmp/all.annots
```

The following lists N of instances (across all individuals), then the term itself:

```
head tmp/all.annots
```
```
 1383442 Sleep stage N2
  875111 Sleep stage N3
  665637 Sleep stage W
  611181 Sleep stage R
  347244 Sleep stage ?
  215226 Oxygen Desaturation
  161644 Oximeter Event
  146038 EEG arousal
  128402 Sleep stage N1
   42179 Obstructive Hypopnea
```

### Check number of unique annotations

```
wc -l tmp/all.annots
```
```
 35816 tmp/all.annots
```

i.e. 35,815 unique labels used in the dataset as a whole

Of these, most are rare: e.g. only 940 occur 10+ times:

```
awk ' $1 > 10 ' tmp/all.annots | wc -l
```
```
    940 
```

Only 38 occur more than 3000 times (e.g. typically at least once per indiv.)

```
awk ' $1 > 3000 ' tmp/all.annots | wc -l
```
```
     38   
```

Sanity check: as calculated above, this shoould include the header row (`description`),
which should occur exactly once per file, i.e. 3,984 times: (i.e. the number of indivs)

```
grep description tmp/all.annots
```
```
  3983 description
```
 
i.e. one file apparently missing the header.   Doing a quick check of rows per file, we
see there is one annotation file that is empty, explaining the above:


```
wc -l ${DATA}/*tsv | sort -nr | tail
```
```
  96 /data/nsrr/working/nchsdb/Sleep_Data//9379_24787.tsv
  95 /data/nsrr/working/nchsdb/Sleep_Data//9880_23149.tsv
  87 /data/nsrr/working/nchsdb/Sleep_Data//4954_9088.tsv
  75 /data/nsrr/working/nchsdb/Sleep_Data//1573_24418.tsv
  67 /data/nsrr/working/nchsdb/Sleep_Data//9202_9148.tsv
  45 /data/nsrr/working/nchsdb/Sleep_Data//2569_2281.tsv
  41 /data/nsrr/working/nchsdb/Sleep_Data//8122_12256.tsv
  20 /data/nsrr/working/nchsdb/Sleep_Data//15028_10648.tsv
   6 /data/nsrr/working/nchsdb/Sleep_Data//13318_16993.tsv
   0 /data/nsrr/working/nchsdb/Sleep_Data//1483_7963.tsv
```

i.e. `1483_7963` has an empty annotation file

nb. some other individuals have very restricted annotations: e.g.

```
cat ${DATA}/13318_16993.tsv
```
```
onset       duration  description
0.47265625  0.0       Montage:NCH Standard with TcCO2_NDx, Ref
0.47265625  0.0       Montage:NCH Standard with TcCO2_NDx, Ref
0.47265625  0.0       Start Recording
0.4765625   0.0       Started Analyzer - ECG
12.5        0.0       Video Recording ON
```

nb. suggests some partial/truncated studies, and/or studies that do not have 
staging available, etc, so we should look out for those.


### Manual review of common annotations and mapping to NSRR terms

Create file `final/annot.mappings` that specifies the annotation remappings, with format:

 - rows = 2 tab-delimited columns
 - NSRR term   -->  NCHSDB term

The NSRR terms should all be listed in a CENTRAL place, under the NSRR repository
that is COMMON ACROSS STUDIES.  

Rremapping .tsv to .annot, with key terms swapped in: all non-recognized terms --> `misc` catch-all annotation class

This reads a dictionary from `final/annot.mapping` and then applies it 
to the *.tsv files to translate terms;  it outputs the .annot formatted
file to the local `annots/` folder 

Option: can either:

 - keep as tab delimited and keep spaces in free-text misc annots
 - OR: can make whitespace delimited, and change spaces to (e.g.) underscores 

by default, Luna (and many other tools) will default to whitespace delimiters
and so this translation makes reading (under default settings) easier

(remove the `gsub()` line below if not wanting this conversion)

```
for f in `ls ${DATA}/*tsv`
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
                        } } ' OFS="\t" final/annot.mapping $f > annots/${id}.annot

done
```

Check that all staging annotations appear to conform to standard 30-seconds:

```
awk -F"\t" ' $1 == "N1" || $1 == "N2" || $1 == "N3" || $1 == "R" || $1 == "W" || $1 == "?" { print $1 , $5 } ' annots/* | sort | uniq -c
```

(output not shown, as too long).  Basically, most are 30-seconds long, i.e. no 10-second/20-second or other variable (multi-epoch) schemes used.
We do see a tail of extra <30 epochs, but _presumably_ these are at the end of each file (i.e. last epoch truncated).  

### Compile a new sample-list

The `--build` option will match the EDFs and the .annot files in `annots/`, matching
on ID (i.e. base file name). 

```
luna --build ${DATA} annots > s.lst
```

```
 wrote 3984 EDFs to the sample list
  3984 of which had 1 linked annotation files
```

```
head s.lst 
```
```
 10000_17728  /data/nsrr/working/nchsdb/Sleep_Data/10000_17728.edf  annots/10000_17728.annot
 10003_26038  /data/nsrr/working/nchsdb/Sleep_Data/10003_26038.edf  annots/10003_26038.annot
 10006_20647  /data/nsrr/working/nchsdb/Sleep_Data/10006_20647.edf  annots/10006_20647.annot
 10009_25600  /data/nsrr/working/nchsdb/Sleep_Data/10009_25600.edf  annots/10009_25600.annot
 10012_22912  /data/nsrr/working/nchsdb/Sleep_Data/10012_22912.edf  annots/10012_22912.annot
 10015_10063  /data/nsrr/working/nchsdb/Sleep_Data/10015_10063.edf  annots/10015_10063.annot
```

### Check coverage of staging annotations

Here we

 - check whether the EDFs are actually discontinuous (not really necessary, all seem to be EDF+C) but doesn't hurt (`SEGMENTS`)
 - check coverage of annotations (staging) over the recordings (`SPANNING`)
 
Test on one sample:

```
luna s.lst 1 -o out.db -s ' SEGMENTS & SPANNING annot=N1,N2,N3,R,W,? '
```

Send all to the LSF:
```
/data/nsrr/bin/runner.sh 20 s.lst . cmd/segments.txt o out/segments out/segments
```

nb. we come across the error below (empty/bad .tsv); so, this doesn't look at all files, but will give a sense of things

```
destrat out/segments.batch00*.db +SEGMENTS  -v NSEGS  | awk ' $2 != 1 ' 
```
(no output, i.e confirming that all EDFs are of only a single, contiguous period of time, not EDF+D)

In terms of annotations, the `NSEGS` variable basically gives similar info to the `SEGMENTS` command, but w.r.t to the interval spanned by a set of annotations. 

```
destrat out/segments.batch00*.db +SPANNING -v NSEGS  | awk ' $2 != 1 ' 
```

```
ID	        NSEGS
13318_16993	0
4480_5926	2
```

i.e. two people have a value other than 1

First person has empty .tsv.  (we saw this above)

For the for second, with 2 segments, implies a gap in staging:

```
destrat out/segments.batch00*.db -i 4480_5926  +SPANNING | awk ' NR == 1 || $1 == "4480_5926" ' | behead
```

```
luna s.lst id=4480_5926 -o out.db -s ANNOTS
```
```
destrat out.db +ANNOTS -r ANNOT/N1,N2,N3,R,W,? INST T
```

```
10500.53515625  30.0    Sleep stage W
10530.53515625  30.0    Sleep stage W
10617.6875      0.0     Headbox Reconnected
10617.6875      0.0     Started Analyzer - ECG
10618.01953125  0.0     Video Recording ON
10620.53515625  30.0    Sleep stage W
10622.69140625  0.0     Started Analyzer - Sleep Events
10622.921875    0.0     Started Analyzer - Data Trends
10622.921875    0.0     Started Analyzer - Sleep Events
10622.921875    0.0     Started Analyzer - ECG
10650.53515625  30.0    Sleep stage W
```

i.e. an actual pause in recording found (one minute w/ no annotations).   This will mean that this file is written as an EDF+D (i.e. as truly discontinuous).

TODO: We could simply append two W annotations onto the .annot 

## EDFs: review channels

We found this message seemingly for all studies read in:

```
EDF+ [/data/nsrr/working/nchsdb/Sleep_Data/10006_20647.edf] did not contain any time-track: adding...
```

Inspect EDF header manually:

```
xxd -l 256 -c 32 /data/nsrr/working/nchsdb/Sleep_Data/10000_17728.edf 
```
```
 0000000: 3020 2020 2020 2020 3030 3030 3030 3020 5820 3031 2d4a 414e 2d32 3030 3120 5858  0       0000000 X 01-JAN-2001 XX
 0000020: 5858 582c 5858 5858 2020 2020 2020 2020 2020 2020 2020 2020 2020 2020 2020 2020  XXX,XXXX                        
 0000040: 2020 2020 2020 2020 2020 2020 2020 2020 2020 2020 2020 2020 5374 6172 7464 6174                          Startdat
 0000060: 6520 3031 2d4a 414e 2d32 3030 3120 5820 5820 5820 2020 2020 2020 2020 2020 2020  e 01-JAN-2001 X X X             
 0000080: 2020 2020 2020 2020 2020 2020 2020 2020 2020 2020 2020 2020 2020 2020 2020 2020                                  
 00000a0: 2020 2020 2020 2020 3031 2e30 312e 3031 3230 2e30 332e 3037 3636 3536 2020 2020          01.01.0120.03.076656    
 00000c0: 4544 462b 4320 2020 2020 2020 2020 2020 2020 2020 2020 2020 2020 2020 2020 2020  EDF+C                           
 00000e0: 2020 2020 2020 2020 2020 2020 3838 3936 2020 2020 3420 2020 2020 2020 3235 2020              8896    4       25  
```

i.e. no EDF Annotations tracks:
 
```
0000100: 454f 4720 4c4f 432d 4d32 2020 2020 2020 454f 4720 524f 432d 4d31 2020 2020 2020  EOG LOC-M2      EOG ROC-M1      
0000120: 454d 4720 4368 696e 322d 4368 696e 3120 4545 4720 4633 2d4d 3220 2020 2020 2020  EMG Chin2-Chin1 EEG F3-M2       
0000140: 4545 4720 4634 2d4d 3120 2020 2020 2020 4545 4720 4333 2d4d 3220 2020 2020 2020  EEG F4-M1       EEG C3-M2       
0000160: 4545 4720 4334 2d4d 3120 2020 2020 2020 4545 4720 4f31 2d4d 3220 2020 2020 2020  EEG C4-M1       EEG O1-M2       
0000180: 4545 4720 4f32 2d4d 3120 2020 2020 2020 4545 4720 435a 2d4f 3120 2020 2020 2020  EEG O2-M1       EEG CZ-O1       
00001a0: 454d 4720 4c4c 6567 2d52 4c65 6720 2020 4543 4720 454b 4732 2d45 4b47 2020 2020  EMG LLeg-RLeg   ECG EKG2-EKG    
00001c0: 536e 6f72 6520 2020 2020 2020 2020 2020 5265 7370 2050 5441 4620 2020 2020 2020  Snore           Resp PTAF       
00001e0: 5265 7370 2041 6972 666c 6f77 2020 2020 5265 7370 2054 686f 7261 6369 6320 2020  Resp Airflow    Resp Thoracic   
0000200: 5265 7370 2041 6264 6f6d 696e 616c 2020 5370 4f32 2020 2020 2020 2020 2020 2020  Resp Abdominal  SpO2            
0000220: 5261 7465 2020 2020 2020 2020 2020 2020 4574 434f 3220 2020 2020 2020 2020 2020  Rate            EtCO2           
0000240: 4361 706e 6f20 2020 2020 2020 2020 2020 5265 7370 2052 6174 6520 2020 2020 2020  Capno           Resp Rate       
0000260: 432d 666c 6f77 2020 2020 2020 2020 2020 5469 6461 6c20 566f 6c20 2020 2020 2020  C-flow          Tidal Vol       
0000280: 5072 6573 7375 7265 2020 2020 2020 2020 2020 2020 2020 2020 2020 2020 2020 2020  Pressure                        
```

i.e. can (and should) be read as standard EDF (add `force-edf=T`) to Luna command line; nb. this speeds up traversal of all EDFs too...

Generate list of channel names; note, this command will also implicitly
check that all .annot files are correctly formated & readable.

```
luna s.lst force-edf=T -o tmp/headers.db -s HEADERS
```

Encountered an error with an annotation for one study:

```
 Processing: 1483_7963 [ #1036 ]
   forcing read as EDF [else remove force-edf=1]
  duration: 07.56.32 | 28592 secs ( clocktime 21.54.48 - 05.51.19 )

  signals: 26 (of 26) selected in a standard EDF file:
   Patient_Event | EOG_LOC-M2 | EOG_ROC-M1 | EMG_Chin1-Chin2 | EEG_F3-M2 | EEG_F4-M1 | EEG_C3-M2 | EEG_C4-M1
   EEG_O1-M2 | EEG_O2-M1 | EEG_CZ-O1 | EMG_LLeg-RLeg | ECG_EKG2-EKG | Snore | Resp_PTAF | Resp_Airflow
   Resp_Thoracic | Resp_Abdominal | SpO2 | Rate | EtCO2 | Capno | Resp_Rate | C-flow
   Tidal_Vol | Pressure
 error : expecting 6/4/3 columns, but found 5
   (use tab-only option to ignore space delimiters)
```

As expected: .annot files are text:
```
file annots/* | head
```
```
annots/10000_17728.annot: ASCII text
annots/10003_26038.annot: ASCII text
annots/10006_20647.annot: ASCII text
annots/10009_25600.annot: ASCII text
annots/10012_22912.annot: ASCII text
annots/10015_10063.annot: ASCII text
```

Something odd about this one: it is a binary file: i.e. 

```
file annots/1483_7963.annot 
```
```
annots/1483_7963.annot: data
```

Checking the original file also: a file full of null (00) characters

```
xxd ${DATA}/1483_7963.tsv | less 
```
```
0000000: 0000 0000 0000 0000 0000 0000 0000 0000  ................
0000010: 0000 0000 0000 0000 0000 0000 0000 0000  ................
0000020: 0000 0000 0000 0000 0000 0000 0000 0000  ................
0000030: 0000 0000 0000 0000 0000 0000 0000 0000  ................
0000040: 0000 0000 0000 0000 0000 0000 0000 0000  ................
0000050: 0000 0000 0000 0000 0000 0000 0000 0000  ................
```

Generate an exclude list to skip this file

```
echo "1483_7963" > tmp/excludes
```

Re-run with this EDF/annot pair excluded

```
luna s.lst force-edf=T exclude=tmp/excludes -o tmp/headers.db -s HEADERS
```

### Check/compile EDF channel names

Extract table of channel label and sample rate per EDF

```
destrat tmp/headers.db +HEADERS -r CH -v SR > tmp/channels
```

nb; by default, Luna swaps spaces with underscores in EDF channel names when reading

```
ID           CH               SR
10000_17728  EOG_LOC-M2       256
10000_17728  EOG_ROC-M1       256
10000_17728  EMG_Chin2-Chin1  256
10000_17728  EEG_F3-M2        256
10000_17728  EEG_F4-M1        256
```

Tabulate how many unique channel labels

```
awk ' NR!=1 { print $2 } ' tmp/channels | sort | uniq -c | sort -nr > tmp/channels2
```

```
wc -l tmp/channels2
```
```
151 
```

Full tabulation here:

```
  3970 EEG_O2-M1
  3970 EEG_O1-M2
  3970 EEG_C3-M2
  3969 Resp_Rate
  3969 Rate
  3969 EtCO2
  3969 Capno
  3968 EEG_F3-M2
  3967 Snore
  3961 EEG_C4-M1
  3959 EEG_F4-M1
  3932 EOG_LOC-M2
  3930 EOG_ROC-M1
  2823 Pressure
  2820 Resp_Thoracic
  2820 Resp_Abdominal
  2819 Resp_Airflow
  2819 ECG_EKG2-EKG
  2818 SpO2
  2818 EEG_CZ-O1
  2817 Tidal_Vol
  2816 Resp_PTAF
  2816 EMG_LLeg-RLeg
  2790 C-flow
  2721 Patient_Event
  2635 EMG_Chin1-Chin2
  1417 TcCO2
  1152 EEG_Cz-O1
  1152 C-Flow
  1151 OSAT
  1148 XFlow
  1148 Resp_Chest
  1148 Resp_Abdomen
  1146 Snore_DR
  1146 Resp_Flow
  1146 Flow_DR
  1146 EMG_RLEG+-RLEG-
  1146 EMG_LLEG+-LLEG-
  1146 EMG_CHIN1-CHIN2
  1146 ECG_LA-RA
  1146 C-Pressure
   136 EEG_Chin1-Chin2
    37 EEG_LOC-M2
    31 EEG_ROC-M1
    21 EMG_Chin3-Chin2
    21 EMG_Chin2-Chin1
    12 EEG_Therm
    12 EEG_Spare
    12 EEG_Snore
    12 EEG_RLeg2
    12 EEG_RLeg1
    12 EEG_Press
    12 EEG_O2
    12 EEG_O1
    12 EEG_M2
    12 EEG_M1
    12 EEG_LLeg2
    12 EEG_LLeg1
    12 EEG_F4
    12 EEG_F3
    12 EEG_EKG2
    12 EEG_EKG1
    12 EEG_E2
    12 EEG_E1
    12 EEG_Chin3
    12 EEG_Chin2
    12 EEG_Chin1
    12 EEG_Chest
    12 EEG_C4
    12 EEG_C3
    12 EEG_Abd
    12 EEG_40
    12 EEG_33
    12 EEG_32
    12 EEG_31
    12 EEG_30
    12 EEG_29
    12 EEG_28
    12 EEG_27
    12 EEG_26
    12 EEG_25
    12 EEG_24
    12 EEG_23
    12 EEG_22
    12 EEG_21
    12 EEG_20
     9 EEG_ROC-M2
     9 EEG_F4-M2
     9 EEG_C4-M2
     7 PTAF
     6 EMG_CHIN1-CHIN3
     5 EEG_Chin3-Chin2
     3 EMG_LLEG-RLEG
     3 EEG_EKG2-EKG
     2 SNORE_DR
     2 Resp_FLOW-Ref
     2 EMG_RAT1-RAT2
     2 EMG_LAT1-LAT2
     2 EEG_EKG-RLeg
     2 ECG_ECGL-ECGR
     1 Thoracic
     1 T6
     1 T5
     1 T4
     1 T3
     1 ROC
     1 RLeg
     1 Resp_Airflow+-Re
     1 PZ
     1 PR
     1 PPG
     1 Position
     1 Pleth
     1 P4
     1 P3
     1 OZ
     1 OSat
     1 O2
     1 O1
     1 M2
     1 M1
     1 LOC
     1 LLeg
     1 FZ
     1 FPZ
     1 Fp2
     1 Fp1
     1 F8
     1 F7
     1 F4
     1 F3
     1 EMG_RLEG-RLEG2
     1 EMG_LLEG-LLEG2
     1 EMG_Chin1-Chin3
     1 EKG2
     1 EKG
     1 EEG_Chin1-Chin3
     1 DC8
     1 DC4
     1 DC3
     1 CZ
     1 Chin3
     1 Chin2
     1 Chin1
     1 C4
     1 C3
     1 Airflow
     1 Abdominal
     1 40
     1 39
     1 38
```

From the above, create a alias set,  e.g. 

```
echo -e  'alias\tO1|EEG_O1-M2' >  final/sigs.alias
echo -e  'alias\tO2|EEG_O2-M1' >> final/sigs.alias

echo -e  'alias\tC3|EEG_C3-M2' >> final/sigs.alias
echo -e  'alias\tC4|EEG_C4-M1' >> final/sigs.alias

echo -e  'alias\tF3|EEG_F3-M2' >> final/sigs.alias
echo -e  'alias\tF4|EEG_F4-M1' >> final/sigs.alias

echo -e  'alias\tLOC|EOG_LOC-M2' >> final/sigs.alias
echo -e  'alias\tROC|EOG_ROC-M1' >> final/sigs.alias
```

```
cat final/sigs.alias
```
```
alias  O1|EEG_O1-M2
alias  O2|EEG_O2-M1
alias  C3|EEG_C3-M2
alias  C4|EEG_C4-M1
alias  F3|EEG_F3-M2
alias  F4|EEG_F4-M1
alias  LOC|EOG_LOC-M2
alias  ROC|EOG_ROC-M1
```


TODO: expand to respiratory channels

TODO: create fuller signals alias file for NCH that includes more than the 8 above


Quick check it converts labels as expected 

```
luna s.lst 1 force-edf=T exclude=tmp/excludes @final/sigs.alias -s DESC
```
```
 EDF filename      : /data/nsrr/working/nchsdb/Sleep_Data/10000_17728.edf
 ID                : 10000_17728
 Clock time        : 20.03.07 - 05.56.10
 Duration          : 09:53:04
 # signals         : 25
 Signals           : LOC[256] ROC[256] EMG_Chin2-Chin1[256] F3[256] F4[256] C3[256]
                     C4[256] O1[256] O2[256] EEG_CZ-O1[256] EMG_LLeg-RLeg[256] ECG_EKG2-EKG[256]
                     Snore[256] Resp_PTAF[256] Resp_Airflow[256] Resp_Thoracic[256] Resp_Abdominal[256] SpO2[256]
                     Rate[256] EtCO2[256] Capno[256] Resp_Rate[256] C-flow[256] Tidal_Vol[256]
                     Pressure[256]
```


### Check sample rates

```
awk ' NR!=1 { print $3 } ' tmp/channels | sort | uniq -c
```
```
     29 255.999733159887
  85475 256
  14902 400
   5753 512
```

Hmm.. see what's going on w/ the non-integer SR

```
grep 255.9997 tmp/channels 
```

 - just a single indiv. (i.e. w/ 29 channels, all data this SR)
 - ID = 5053_2167

Check EDF record durations:

```
destrat tmp/headers.db +HEADERS -v REC_DUR | awk ' NR!=1 { print $2 } ' | sort | uniq -c 
```
```
      2 12
    793 2
      1 3.59766
   3187 4
```

Note one individual w/ strange EDF record duration... might this per chance be 5053_2167...?

```
luna s.lst id=5053_2167 -s SUMMARY
```

```
  Rec. dur. (s)  : 3.59766

 Signal 2 : [EOG_LOC-M2]
         # samples per record : 921
         transducer type      : 
         physical dimension   : uV
         min/max (phys)       : 8711/-8711
         EDF min/max (phys)   : 8711/-8711
         min/max (digital)    : -32768/32767
         EDF min/max (digital): -32768/32767
         pre-filtering        :                                                                                 
```

i.e. Yes.  5053_2167 has 921 samples per record; This _almost_ implies
256 Hz but 921/256 = 3.59765625, not 3.59766 ... Thus the imperfect
sample rate of 255.999 not 256.

 - nb. silly record duration to use, given imprecise floating-point represent (8 chars) 
     is not guaranteed to be able to track this precisely in EDF...
 - could try re-ALIGNing with the others, but a concern that it will not be read properly 
 -  i.e. given one can only read an integer number of samples, it 
       might read only 255 out of every 256 sample-points
 -  i.e. this is effectively a broken EDF, so we will exclude

Add this person to the exclude list:

```
echo "5053_2167" >> tmp/excludes
```


## Check staging annotations

For convenience, let's aggregate all project-specific options in a single 
parameter file:

```
echo -e  'force-edf\tT' > param
echo -e  'exclude\ttmp/excludes' >> param 
cat final/sigs.alias >> param
```

As always, try 1 first:

```
luna s.lst @param 1 -o tmp/hypno.db -s HYPNO 
```

```
 CMD #1: HYPNO
   options: sig=*
  set epochs to default 30 seconds, 1186 epochs
  *** found 37 epoch(s) of 1186 with conflicting spanning annotations
  *** check that epochs and annotations align as intended
  *** see EPOCH 'start-annot' or 'offset' options
```

Nb. This info is also stored in the CONF variable of the HYPNO output:

```
destrat tmp/hypno.db +HYPNO -v CONF
```
```
ID           CONF
10000_17728  37
```

What's going on here?  Examine the annotation files: 


```
 # misc | Misc annots | note[txt]
 misc    .       .       1.03515625      +0.0    Montage:NCH_PSG_STANDARD,_Ref
 misc    .       .       1.03515625      +0.0    Start_Recording
 misc    .       .       1.03515625      +0.0    Montage:NCH_PSG_STANDARD,_Ref
 misc    .       .       1.1640625       +0.0    Recording_Analyzer_-_ECG
 misc    .       .       1.6171875       +0.0    Video_Recording_ON
 misc    .       .       56.78125        +0.0    Montage:Channel_Test_Referential
 misc    .       .       56.78125        +0.0    Montage:Channel_Test_Referential
 misc    .       .       117.109375      +0.0    Montage:NCH_PSG_STANDARD
 misc    .       .       117.109375      +0.0    Montage:NCH_PSG_STANDARD
 misc    .       .       2801.015625     +0.0    Impedance_at_5_kOhm
 misc    .       .       3150.890625     +0.0    Impedance_at_5_kOhm
 misc    .       .       3950.5859375    +0.0    left_foot_stim
 misc    .       .       3958.57421875   +0.0    right_foot_stim
 misc    .       .       3964.08203125   +0.0    eyes_closed
 misc    .       .       3972.16796875   +0.0    eyes_open
 W       .       .       4080.9609375    +30.0   .
 lights_off      .       .       4087.16796875   +0.0    .
 N1      .       .       4110.9609375    +30.0   .
 N2      .       .       4140.9609375    +30.0   .
 misc    .       .       4141.2890625    +0.0    sleep_onset
 N2      .       .       4170.9609375    +30.0   .
```

Extracting the staging annots only from the above (below), we see that each is 30
seconds, but starts at an arbitrary sample-point. That is, it is not
aligned with EDF records

```
 W       .       .       4080.9609375    +30.0   .
 N1      .       .       4110.9609375    +30.0   .
 N2      .       .       4140.9609375    +30.0   .
 N2      .       .       4170.9609375    +30.0   .
```

This is a problem potentially, for two reasons:

First, in Luna, by default epochs are defined starting from 0 seconds 
when reading annotations to determine what stage an epoch should be designated as, 
this is a problem if > 1 (conflicting) stage annotations are present, as will happen
if epochs and stage annotations are not aligned

This can be remedied easily by adding the 'align' option to EPOCH: 
instead of starting EPOCHs at 0.0, it will start counting at the start of
the first of the listed annotations

```
luna s.lst @param 1 -o tmp/hypno.db -s 'EPOCH align=?,N1,N2,N3,W,R & HYPNO'
```
```
destrat tmp/hypno.db +HYPNO -v CONF
```
```
ID           CONF
10000_17728  0
```


This will obviously impact estimates of stage duration, etc, as conflicting epochs
are otherwise set to missing (`?`) so important to do.


_However, more fundamentally_ and arising as a property of EDF (i.e. rather than 
any Luna-specific workflow), this also means that EDFs are not aligned with EDF records.

This will make it difficult to restructure the EDF, e.g. if wanting to write out a new EDF
(with harmonized labels, canonical signals, re-sampled sample rates, or artifactual/bad epochs
removed, etc. and keep as EDF format) 

i.e. can only write whole EDF records to EDF.  If the record size is 4.0 (or 1.0, or 12.0 seconds, etc)
then this will not be possible, if wanting to take an epoch that spans, say, 22.565523 to 52.565523 seconds

Therefore, we will use the `ALIGN` command to make a new EDF with fixed EDF record size, and that 
perfectly aligns annotations and EDF records;  this will also mean that default Luna epochs (i.e. starting 
at 0.0 seconds) will align too without the need for further work. 

We'll do this in two steps: 

1) make canonical EDFs and change EDF record size to 1.0 across all EDFs/channels

```
cat cmd/canonical.txt
```
```
 % Create canonical signals 

 CANONICAL file=sigs.canonical group=NCH 

 % Drop all channels except these new canonical signals

 SIGNALS req=csC3,csC4,csO1,csO2,csF3,csF4
 
 % force as a standard EDF (instead of EDF+) if possible
 
 EDF

 % change EDF record size to 1 second; this also writes the EDF immediately

 RECORD-SIZE dur=1 edf-dir=edfs/ edf-tag=can6 sample-list=cs.lst with-annots
```

Canonical signals are defined here:

```
cat sigs.canonical
```
```
NCH C3 EEG_C3-M2 . 256 uV
NCH C4 EEG_C4-M1 . 256 uV
NCH F3 EEG_F3-M2 . 256 uV
NCH F4 EEG_F4-M1 . 256 uV
NCH O1 EEG_O1-M2 . 256 uV
NCH O2 EEG_O2-M1 . 256 uV
```

i.e. just use original EDF labels, and create six new (EEG) canonical signals,  which will be csC3, csC4, etc

Test on one sample:

```
luna s.lst 1 @param < cmd/canonical.txt
```

Run on all samples:

```
rm -rf cs.lst

/data/nsrr/bin/runner.sh 20 s.lst param cmd/canonical.txt o out/mk-canon6 out/mk-canon6
```

This creates a new sample list : `cs.lst` 

```
 EDF filename      : edfs/10000_17728-can6.edf
 ID                : 10000_17728
 Clock time        : 20.03.07 - 05.56.10
 Duration          : 09:53:04
 # signals         : 6
 Signals           : csC3[100] csC4[100] csF3[100] csF4[100] csO1[100] csO2[100]
```

```
head cs.lst
```
```
 10000_17728  edfs/10000_17728-can6.edf   annots/10000_17728.annot
 1096_18565   edfs/1096_18565-can6.edf    annots/1096_18565.annot
 10003_26038  edfs/10003_26038-can6.edf   annots/10003_26038.annot
 10975_5821   edfs/10975_5821-can6.edf    annots/10975_5821.annot
 10978_6367   edfs/10978_6367-can6.edf    annots/10978_6367.annot
```

i.e. this sample-list points to the new EDFs and also tracks the old annotations

Step 2: in the new set of EDFs (with uniform EDF record sizes), re-align with staging annotations


The only requirement is that any selected annotations to realign to are:

 - not overlapping
 - have duration that is an integer multiple of the EDF record size

```
cat cmd/align.txt
```
```
 ALIGN align=N1,N2,N3,R,W,? 
       edf-tag=aligned 
       edf-dir=final/ 
       annot-out=final/^-can6-aligned.annot
```

nb. as jobs running in parallel and is a fast job, good practice to write to 
separate sample-lists and then concatenate. i.e. otherwise, the sample-list 
file may corrupt if two processes are writing simultaneously.  (Unlikely to 
happen, but can if thousands of files and dozens of parallel jobs...) 

Test on 1 / run on all:

```
luna cs.lst 1 @param < cmd/align.txt
```

```
/data/nsrr/bin/runner.sh 20 cs.lst param cmd/align.txt o out/mk-align out/mk-align
```

Build a new sample list for these realigned EDFs

```
luna --build final > cs-aligned.lst
```

```
wrote 3958 EDFs to the sample list
  3958 of which had 1 linked annotation files
```
```
head cs-aligned.lst
```
```
10000_17728-can6-aligned  final/10000_17728-can6-aligned.edf  final/10000_17728-can6-aligned.annot
10003_26038-can6-aligned  final/10003_26038-can6-aligned.edf  final/10003_26038-can6-aligned.annot
10006_20647-can6-aligned  final/10006_20647-can6-aligned.edf  final/10006_20647-can6-aligned.annot
10009_25600-can6-aligned  final/10009_25600-can6-aligned.edf  final/10009_25600-can6-aligned.annot
...
```

nb. as IDs taken from EDF root, want to change those back

```
awk ' { gsub( "-can6-aligned" , "" , $1 ) ; print $1 , $2 , $3 } ' OFS="\t" cs-aligned.lst > final.lst
```

```
head final.lst
```
```
10000_17728  final/10000_17728-can6-aligned.edf  final/10000_17728-can6-aligned.annot
10003_26038  final/10003_26038-can6-aligned.edf  final/10003_26038-can6-aligned.annot
10006_20647  final/10006_20647-can6-aligned.edf  final/10006_20647-can6-aligned.annot
10009_25600  final/10009_25600-can6-aligned.edf  final/10009_25600-can6-aligned.annot
10012_22912  final/10012_22912-can6-aligned.edf  final/10012_22912-can6-aligned.annot
...
```


## Sanity check assembled EDFs/ANNOTs

At this point, we've assembled everything; we have N=3,920 EDFs/ANNOTs
in the final set. i.e. some dropped if no EEG, etc

### Check these are all standard (i.e. continuous/contiguous) recordings now

First, check whether these are standard EDFs or not. Based on the
`SPANNING` analysis above, this should be the case (except for one
person), but let's check in the final dataset anyway.

When writing an EDF, by default Luna will 'downcast' from EDF+D to 
standard EDF if possible, i.e. no discontinuities and no `EDF Annotations`
channels (other than the EDF record time-track)

```
luna final.lst -o tmp/qc1.db -s 'SEGMENTS & SPANNING annot=N1,N2,N3,R,W,?'
```

Looking at the number of _segments_ per individual: third col is `NSEGS` so find any values > 1 :
```
destrat tmp/qc1.db +SEGMENTS | awk ' $3 != 1 '
```
```
ID	        NGAPS	NSEGS
4480_5926	1	2
```

i.e. this finds a single person who has two segments (and will be an
EDF+D).  All other recordings should be as standard EDF.  This was the
person we identified above.  It doesn't really matter, but just to
clean things up, we might want to go back and insert two "Sleep stage
W" (or unknown stage) just so we get a continuous standard EDF out of
the process, and don't have to worry about a single EDF+D.


### Check coverage of annotations

Because of the realignment, _by definition_ we should expect that all epochs have a corresponding sleep stage
(i.e. it is how they were selected...).   Still, does not hurt to check:

```
destrat tmp/qc1.db +SPANNING 
```
```
destrat tmp/qc1.db +SPANNING -v SPANNED_PCT | awk ' { print $2 } '  | sort | uniq -c
```
```
      1 SPANNED_PCT
   3920 100
```
i.e. all EDFs have 100% of their duration spanned by one of N1,N2,N3,R,W,?.

Also, we can confirm that none of these sleep annotations overlap with each other: (the `ANNOT_OVERLAP` variable
is set to `YES` if so):
```
destrat tmp/qc1.db +SPANNING -v ANNOT_OVERLAP  | awk ' { print $2 } ' | sort | uniq -c
```
```
     1 ANNOT_OVERLAP
  3920 NO
```

### Check sleep macro-architecture

```
luna final.lst -o out/hypno.db -s HYPNO
```

e.g. 2734 people w/ at least 6 hours of sleep recording:
```
destrat out/hypno.db +HYPNO -v TST  | awk ' NR=1 && $2 > 6*60 ' | wc -l
```
```
2734
```

People w/ short or no recording sleep

```
destrat out/hypno.db +HYPNO -v TST  | awk ' NR=1 && $2 < 30 ' 
```
```
16009_6094	0
16261_17596	21.5
16933_26176	0
17053_2659	24.5
17839_2824	3.5
18853_18163	26.5
1885_21658	0
2065_5257	24.5
235_25954	0
3559_22582	0
9559_12535	0
9766_20215	17.5
```

nb. these are just some of the basic sanity checks we might do
(i.e. to just confirm that most recordings are mostly staged as sleep,
etc.)  However, fuller descrptions of sleep macro-architecture should
be performed later on the final dataset, i.e. a separte activity from
the this more basic data reformatting/sanity-checking exercise.


### Look at EDF headers

```
luna final.lst -o out/headers.db -s HEADERS
```

As expected, we have one EDF+D file because of the gapped staging (can fix).
```
destrat out/headers.db +HEADERS -v EDF_TYPE  | awk ' NR!=1 { print $2 } ' | sort  | uniq -c
```
```
  3957 EDF
     1 EDF+D
```

All have the correct 6 signals
```
destrat out/headers.db +HEADERS -v NS  | awk ' NR!=1 { print $2 } ' | sort  | uniq -c
```
```
   3957 6
      1 7
````
(note: the one extra w/ 7 is the extra `EDF Annotations` channel that an EDF+D has)



### A final sanity check of realigned signal data 

i.e. did the realignment/shifting process faithfully track annotations and signals, etc.

Do for a single ID (=10000_17728).

Old EDF start time (`luna -s DESC`): 
```
 duration: 09.53.04 | 35584 secs ( clocktime 20.03.07 - 05.56.10 )
```

New (re-aligned) EDF times:
```
 duration: 08.45.00 | 31500 secs ( clocktime 21.11.07 - 05.56.06 )
```

From the output of the `align.txt` procedure for this individual, we see the shift in time (from the console log):
```
  applying a offset of -4080.96 to all annotations when writing out
```
i.e. because that was when the first wake epoch was observed (from the .annot):
```
W       .       .       4080.9609375    +30.0   .
```

Look at the annotations to check clock-time is preserved despite shifting everything.  Pick the first arousal annotation.  From the original .tsv:

```
grep arousal ${DATA}/10000_17728.tsv | head -1 
```
```
6917.6484375	12.90625	EEG arousal
```

From the initial .annot:
```
grep arousal annots/10000_17728.annot | head -1 
```
```
arousal	.	.	6917.6484375	+12.90625	.
```


Find this first annotation from both the original and realigned datasets: first in elapsed seconds
```
luna s.lst     id=10000_17728 -s WRITE-ANNOTS annot=arousal file=f1.annot
luna final.lst id=10000_17728 -s WRITE-ANNOTS annot=arousal file=f2.annot
```
which gives (from the orig, then final respectively):
```
arousal	       .	      .	 6917.648     6930.554	    .	
arousal	       .	      .	 2836.688     2849.594	    .
```
i.e. as expected, we have a shift of times (6917.648 - 2836.688 = 4080.96 seconds)


Now, look at the clocktime encoding (i.e. add `hms`) and check this looks okay: (i.e. the elapsed time is different
but we expect the clocktime to be (more or less... see note below) the same:

```
luna s.lst     id=10000_17728 -s WRITE-ANNOTS annot=arousal hms file=f1.annot
luna final.lst id=10000_17728 -s WRITE-ANNOTS annot=arousal hms file=f2.annot
```
```
arousal	       .	      .	 21:58:24.648 21:58:37.554  .
arousal	       .	      .	 21:58:23.688 21:58:36.594  .
```

i.e. basicallty the same, although there is shift 0.96 seconds backwards in time.
This is an unavoidable consequence of EDF starttime format only allowing integer
seconds in the header: i.e. it is impossible to accurately encode the case when
the study "really" starts at some fraction of a second, if we are to rely on EDF headers (which,
practically, we do want/need to do).

Overall, this expected small shift in clock-times is not a big deal,
but just noting here that it is a consequence of realignment.  Just
noting it here so we don't get confused later...


## Spectral analysis of EEG 

As a quick check of the EEG - just for now looking at C3:

```
luna final.lst -o peaks100.db -s 'MASK ifnot=N2 & RE & EPOCH require=20 & PSD sig=csC3 spectrum peaks max=100'
```

```
luna final.lst -o peaks25.db -s 'MASK ifnot=N2 & RE & EPOCH require=20 & PSD sig=csC3 spectrum peaks '
```

[ to be completed ... ] 

## Summary

TODO: brief summary of the above process, summarizing main observations and any lessons learned, etc.

 - basic channel/annotation label harmonization (--> updating NSRR central lists)
 - gaps in staging
 - empty annotation files
 - checks of sample rates and EDF duration
 - alignment of annotations and EDFs 
 - EEG line noise


Actions (to be) completed:

 - [x] Check basic file count, integrity, etc
 - [x] Realign EDF records to annotation
 - [x] Create canonical	EEG EDFs w/ harmonized staging annotations
 - [x] Initial spectral analysis of sleep EEG
 - [x] Extract a broader list of harmonized annotations (beyond staging) and map to NSRR terms
 - [x] Update NSRR master list of annotation terms
 - [ ] Extract canonical respitory signals
 - [ ] Consider options for dealing w/ line noise artifact in these recordings
 - [ ] Create a	brief README summarizing these steps (and pointing to this file in the repo)
 - [ ] Post EDFs, annotations and the README on NSRR


