# WSC Cohort
 - NSRR Harmonization v1
 - Created April 6 2021

## Original Data
All original EDF's and annotations are located here
```
TWIN_DIR=/data/nsrr/working/wsc-scoring-annotations/twin
GAMMA_DIR=/data/nsrr/working/wsc-scoring-annotations/gamma
```

EDF Filename/ID format: three fields, hyphen separated
```
COHORT_NAME VISIT_ID SLEEP_STUDY_ID
```
Expecting n=743 studies/annotations TWIN

Expecting n=1827 studies/annotations GAMMA

```
ls ${TWIN_DIR}/*.edf | wc -l
743
ls ${GAMMA_DIR}/*.edf | wc -l
1827
```

Scoring Filename format: three fields, hyphen separated
```
TWIN: wsc-visitID-studyIDallScore.txt (ex: wsc-visit3-46950allScore.txt)

ls ${TWIN_DIR}/*allScore.txt | wc -l
743
```

```
GAMMA: wsc-visitID-studyIDsco.txt (ex: wsc-visit1-10119sco.txt)
       wsc-visitID-studyIDstg.txt (ex: wsc-visit1-10119stg.txt)
       
ls ${GAMMA_DIR}/*sco.txt | wc -l
1827
ls ${GAMMA_DIR}/*stg.txt | wc -l
1827
```


## Local working directories
Temporary scratch space
```
mkdir t_tmp 
mkdir g_tmp
```
Harmonzied annotations
```
mkdir twin_annots
mkdir gamma_annots
```
Final EDFs and annotations, i.e. to be posted to NSRR
```
mkdir twin_final
mkdir gamma_final
```

## Checking ID's
First, get a list of IDs from the EDF filenames:
```
ls ${TWIN_DIR} | grep "\.edf" | sed 's/\.edf//g' > t_tmp/ids
ls ${GAMMA_DIR} | grep "\.edf" | sed 's/\.edf//g' > g_tmp/ids
```

Then confirm the assumption that all IDs have correct form, i.e. three hyphen-delimited fields:
```
awk -F"-" ' { print NF } ' t_tmp/ids | sort | uniq -c
743 3

awk -F"-" ' { print NF } ' g_tmp/ids | sort | uniq -c
1827 3
```

Check that EDF and scoring file ID's line up exactly:
```
ls ${TWIN_DIR} | grep "allScore\.txt" | sed 's/allScore\.txt//g' > t_tmp/ids.allScore
paste t_tmp/ids t_tmp/ids.allScore | awk ' $1 != $2 ' 

ls ${GAMMA_DIR} | grep "sco\.txt" | sed 's/sco\.txt//g' > g_tmp/ids.sco
ls ${GAMMA_DIR} | grep "stg\.txt" | sed 's/stg\.txt//g' > g_tmp/ids.stg
paste g_tmp/ids g_tmp/ids.sco | awk ' $1 != $2 '
paste g_tmp/ids g_tmp/ids.stg | awk ' $1 != $2 '
```
(no output: i.e. good, this implies that all EDF and scoring file names match exactly)



## Annotation Summary

All the below mapping are listed in [nsrr-mapping](https://gitlab-scm.partners.org/zzz-public/nsrr/-/blob/master/common/annotationdefinitions.md) page

### GAMMA
Gamma annots summary: stg files have all stages information. And sco files have all annots other than stages.
```
cat /data/nsrr/working/wsc-scoring-annotations/gamma/*sco.txt | awk -F"\t" ' {print $5} ' | sort | uniq -c | sort -nr > g_tmp/gamma.annots
```
Here we can see the output of Gamma Annots (except stages)
```
head g_tmp/gamma.annots 
1836061 
 226701 SaO2
 153936 Hypopnea 
 152270 LMA 
 114107 LM 
  30707 Hypopnea
  27154 Obs Apnea
  20434 LM
   7109 SPON Arousal
   6982 Spon Arousal
   6881 Resp Arousal
   6145 RESP Arousal
   5131 Central Apnea 
   5105 OBS Apnea
   1829 Marker Text
   1618 Obst. Apnea
   ...
```
Looks like there are 2 extra header rows (1829 Marker Text), let's find these files:
```
for i in ${GAMMA_DIR}/*sco.txt
do
  header_count=`grep -e "Marker" $i | wc -l`
  if [[ ! ${header_count} -eq 1 ]];then
    echo "File with header row other than 1: ${i}"
  fi
done
```
Files with header row other than 1 are:

 - /data/nsrr/working/wsc-scoring-annotations/gamma/wsc-visit1-95797sco.txt
 - /data/nsrr/working/wsc-scoring-annotations/gamma/wsc-visit1-97405sco.txt

Manually verified that these two files have double header rows

Let's dump all Gamma Stages now,
```
cat /data/nsrr/working/wsc-scoring-annotations/gamma/*stg.txt | awk -F"\t" ' {print $5} ' | sort | uniq -c | sort -nr > g_tmp/gamma.stages
```
here we can see the output of Gamma Stages 
 ```
 884641 2
 303191 0
 225076 5
 135180 1
  99151 7
  90451 3
   4989 4
   1826 User-Defined Stage
     33 6
````
Looks like one file is missing header, Let's search that file:

```
for i in /data/nsrr/working/wsc-scoring-annotations/gamma/*stg.txt
do
  is_header=`head -1 ${i} | grep -e "User-Defined Stage" | wc -l`
  if [[ ${is_header} -eq 0 ]];then
    echo "File missing header is ${i}"
  fi
done
```
File missing header is /data/nsrr/working/wsc-scoring-annotations/gamma/wsc-visit1-62741stg.txt

### TWIN
Twin annots summary: Since data is tab and hyphen delimited, we have to perform separate arousal, resp event and stage calculations
```
cat /data/nsrr/working/wsc-scoring-annotations/twin/*allScore.txt | tr -d '\r' | sed 's/ - /\t/g' | awk -F"\t" ' {print $2} ' | grep -v -e "AROUSAL" -e "RESPIRATORY EVENT" -e "STAGE"  > tmp.txt
cat /data/nsrr/working/wsc-scoring-annotations/twin/*allScore.txt | tr -d '\r' | sed 's/ - /\t/g' | awk -F"\t" ' {print $2,$4} ' OFS="\t" | grep -e "AROUSAL" -e "RESPIRATORY EVENT"  >> tmp.txt
cat /data/nsrr/working/wsc-scoring-annotations/twin/*allScore.txt | tr -d '\r' | sed 's/ - /\t/g' | awk -F"\t" ' {print $2,$3} ' OFS="\t" | grep -e "STAGE"  >> tmp.txt
cat tmp.txt | sort | uniq -c | sort -nr > t_tmp/all.annots
```

Let's look at the Twin annots
```
head t_tmp/all.annots
 157976 LM
  84620 DESATURATION
  76455 RESPIRATORY EVENT	HYPOPNEA
  41739 STAGE	N2
  37579 AROUSAL	LM
  30258 STAGE	N1
  24687 STAGE	W
  17082 POSITION
  14975 STAGE	N3
  14212 AROUSAL	RESPIRATORY EVENT
  13583 AROUSAL	SPONTANEOUS
   8903 STAGE	R
   7089 RESPIRATORY EVENT	OBSTRUCTIVE APNEA
   2326 NEW MONTAGE
    984 RESPIRATORY EVENT	CENTRAL APNEA
    967 STAGE	NO STAGE
    ...
```
It is a huge list with many annotations representated as scoring notes. So total count of twin annots are:
```
cat t_tmp/all.annots | awk '{s+=$1}END{print s}'
555915
```

## Manual review of common annotations and mapping to NSRR terms

### Twin 
Mappings for annot file creation
| annot | NSRR |
| ------ | ------ |
|DESATURATION|desat|
|AROUSAL SPONTANEOUS| arousal_spontaneous|
|AROUSAL RESPIRATORY EVENT|   arousal_respiratory|
|AROUSAL LM|  arousal_lm|
|RESPIRATORY EVENT HYPOPNEA|  hypopnea|
|RESPIRATORY EVENT OBSTRUCTIVE APNEA| apnea_obstructive|
|RESPIRATORY EVENT CENTRAL APNEA| apnea_central|
|RESPIRATORY EVENT MIXED APNEA|   apnea_mixed|
|LM|  lm|
|LIGHTS OUT|  lights_off|
|LIGHTS ON|   lights_on|
|PAUSED|  paused|
|START RECORDING| startrecording|

All other annotation events will be mapped as "misc" (with spaces converted to underscore)

Mapping for eannot file creation:
|eannot|NSRR|
|------|------|
|STAGE N1| N1|
|STAGE N2| N2|
|STAGE N3| N3|
|STAGE W|W|
|STAGE R|R|
|STAGE NO STAGE|?|


### GAMMA

Mappings for annot file creation:
| annot | NSRR |
| ------ | ------ |
|SaO2|desat|
|SponArousal or sponarousal or SPONArousal|arousal_spontaneous|
|Arousal|arousal_standard|
|RespA or resparousal or RespArousal or RESPArousal|arousal_respiratory|
|LMA|arousal_lm|
|Hypopnea or CentralHypopnea or Obst.Hypopnea|hypopnea|
|OA or ObsApnea or Obst.Apnea or ObstApnea or OBSApnea or Apnea|apnea_obstructive|
|CA or CentralApnea|apnea_central|
|MA or MixedApnea|apnea_mixed|
|LM|lm|
All other annotation events will be mapped as "misc" (with spaces converted to underscore)


Mapping for eannot file creation:
|eannot|NSRR|
|------|------|
|0|W|
|1|N1|
|2|N2|
|3|N3|
|4|N3|
|5|R|
|6|movement|
|7|?|

Now run .annot and .eannot file creation script with output in g_final/ and t_final/ folders respectively



## Compile new sample list
The --build option will match the EDFs and the .annot files in annots/, matching
on ID (i.e. base file name).

```
luna --build ${GAMMA_DIR} g_annots > g_s.lst
```
wrote 1827 EDFs to the sample list
  1827 of which had 2 linked annotation files

```
luna --build ${TWIN_DIR} t_annots > t_s.lst
```
wrote 743 EDFs to the sample list
  743 of which had 2 linked annotation files


## Check coverage of staging annotations
Here we
 - check whether the EDFs are actually discontinuous (not really necessary, all seem to be EDF+C) but doesn't hurt (SEGMENTS)
 - check coverage of annotations (staging) over the recordings (SPANNING)

Write a command file as,
```
cat cmd/segments.txt
SEGMENTS
SPANNING annot=N1,N2,N3,R,W,?,movement
```

Run the job:
```
/data/nsrr/bin/runner.sh 20 t_s.lst . cmd/segments.txt o t_out/segments t_out/segments
/data/nsrr/bin/runner.sh 20 g_s.lst . cmd/segments.txt o g_out/segments g_out/segments
```

```
destrat g_out/segments.batch00*.db +SEGMENTS  -v NSEGS  | awk ' $2 != 1 '
```
no output, looks fine

```
 destrat g_out/segments.batch00*.db +SPANNING -v NSEGS  | awk ' $2 != 1 '
```
no output, looks fine

```
destrat t_out/segments.batch00*.db +SEGMENTS  -v NSEGS  | awk ' $2 != 1 '
```
no output, looks fine

```
destrat t_out/segments.batch00*.db +SPANNING  -v NSEGS  | awk ' $2 != 1 '
ID	NSEGS
wsc-visit4-45113	0
```
i.e. One people have a value other than 1. 

Add above conflicting file into excludes list
```
echo "wsc-visit4-45113" > t_tmp/excludes
```



## Check and compile EDF channels

For gamma stduy, We already know that there are following issues with:

1. Negative dur values:
 - wsc-visit1-10191
 - wsc-visit1-11162
 - wsc-visit1-41115
 - wsc-visit1-75614
 - wsc-visit1-82488
 - wsc-visit1-89175
 - wsc-visit2-12325
 - wsc-visit2-64948

2. Missing start time:
 - wsc-visit2-77724

So, add them to exclude list when calculating headers: g_tmp/excludes
```
echo "wsc-visit1-10191" > g_tmp/excludes
echo "wsc-visit1-11162" >> g_tmp/excludes
echo "wsc-visit1-41115" >> g_tmp/excludes
echo "wsc-visit1-75614" >> g_tmp/excludes
echo "wsc-visit1-82488" >> g_tmp/excludes
echo "wsc-visit1-89175" >> g_tmp/excludes
echo "wsc-visit2-12325" >> g_tmp/excludes
echo "wsc-visit2-64948" >> g_tmp/excludes
echo "wsc-visit2-77724" >> g_tmp/excludes
```

Let's compile EDF's with annotations using sample list,
```
luna g_s.lst force-edf=T exclude=g_tmp/excludes -o g_tmp/headers.db -s HEADERS
```

Found One error
```
Processing: wsc-visit2-89696 [ #1624 ]
error : expecting 689 epoch annotations, but found 841
```
wsc-visit2-89696  EDF clock time and stg file epoch based clock doesn't align
EDF: duration: 05.44.30 | 20670 secs ( clocktime 00.16.25 - 06.00.54 ) -> 689 epochs
stg: 841 Epochs

Add this file to exclude list and re-run,
```
echo "wsc-visit2-89696" >> g_tmp/excludes
```

Upon re-run, Command ran without any errors. Now we have all the Gamma Header information extracted
```
luna t_s.lst force-edf=T  -o t_tmp/headers.db -s HEADERS
```
Command executed without any errors, so now we have all the Twin Header information

Now, extract table of channel label and sample rate per EDF

```
destrat g_tmp/headers.db +HEADERS -r CH -v SR > g_tmp/channels
destrat t_tmp/headers.db +HEADERS -r CH -v SR > t_tmp/channels
```

Tabulate how many unique channel labels

```
awk ' NR!=1 { print $2 } ' g_tmp/channels | sort | uniq -c | sort -nr > g_tmp/channels_uniq
wc -l g_tmp/channels_uniq
34
```

Full tabulation here
```
   1817 SNORE
   1817 RIB_CAGE
   1817 ECG
   1817 ABDOMEN
   1813 CHIN_EMG
   1807 POSITION
   1803 R_EOG
   1803 L_EOG
   1701 LEG_EMG
   1679 SAO2
   1510 L_OCC
   1510 L_CENT
   1450 NASAL_PRES
   1014 ORAL_FLOW
   1014 NASAL_FLOW
    457 FLOW.1
    457 FLOW
    367 NASAL_PRESS
    344 AIRFLOW
    307 OCC
    307 FRONTAL
    307 CENTRAL
    112 SaO2
    112 LEGS
     35 AIRFLOW.1
     26 SA02
     14 REOG
     14 LEOG
     10 Position
      4 LEG
      4 FLOW.2
      4 CHIN
      2 AIR_FLOW
```

```
awk ' NR!=1 { print $2 } ' t_tmp/channels | sort | uniq -c | sort -nr > t_tmp/channels_uniq
wc -l t_tmp/channels_uniq
57
```

Full tabulation here
```
    743 SUM
    743 Chest
    743 Airflow
    743 Abd
    742 NasalP
    741 EKG1-EKG2
    736 SpO2
    735 SNORE
    735 POS
    734 LLEG1-RLEG1
    734 LEOG-M2
    730 REOG-M1
    730 Chin1-Chin2
    728 F3-M2
    724 O1-M2
    722 C3-M2
     75 PZ-M2
     75 CZ-M2
     74 FZ-M2
     11 O1-AVG
     10 Chin1-Chin3
      9 REOG-AVG
      8 F3-AVG
      7 SaO2
      6 C4-M1
      4 REOG-M2
      4 O2-M1
      4 F4-M1
      4 F4-AVG
      3 O2-M2
      3 Chin3-Chin2
      3 C4-M2
      2 RLEG1-RLEG2
      2 PZ-AVG
      2 LLEG2-RLEG2
      2 LLEG2-RLEG1
      2 LLEG1-RLEG2
      2 FZ-AVG
      2 CZ-AVG
      2 CPRES
      2 CFLOW
      1 PZ-CZ
      1 O1-M1
      1 LLEG1-LLEG2
      1 LLEG1-EKG2
      1 FZ-M1
      1 F4-M2
      1 F3-M1
      1 F3-AVG.1
      1 EKG1-AVG
      1 CZ-M1
      1 C4-AVG
      1 C3-M1
      1 C3-AVG.1
      1 Airflow.1
```


From the above, create a alias set into g_final/sigs.alias and t_final/sigs.alias



## Check Sample rate
Gamma: 
```
awk ' NR!=1 { print $3 } ' g_tmp/channels | sort | uniq -c
  27308 100
   1764 128
```

Twin:
```
 awk ' NR!=1 { print $3 } ' t_tmp/channels | sort | uniq -c
  12113 200
```

## Check EDF duration

Gamma:
```
destrat g_tmp/headers.db +HEADERS -v REC_DUR | awk ' NR!=1 { print $2 } ' | sort | uniq -c
```

```
      8 10
      2 6
     79 7
    654 8
```
No floating values, output looks fine

Twin:
```
destrat t_tmp/headers.db +HEADERS -v REC_DUR | awk ' NR!=1 { print $2 } ' | sort | uniq -c
```

```
   1614 1
```
No floating values, output looks fine


## Check sleep macro-architecture

For convenience, let's aggregate all project-specific options in a single
parameter file:

GAMMA:
```
echo -e  'force-edf\tT' > g_param
echo -e  'exclude\tg_tmp/excludes' >> g_param
cat final/sigs.alias >> g_param
```
As always, try 1 first:
```
luna g_s.lst @g_param 1 -o g_tmp/hypno.db -s HYPNO
destrat g_tmp/hypno.db +HYPNO -v CONF
```
No conflicts, proceed to submit batch job`

TWIN:
```
echo -e  'force-edf\tT' > t_param
echo -e  'exclude\tt_tmp/excludes' >> t_param
cat sigs.alias >> t_param
```


As always, try 1 first:
```
luna t_s.lst @t_param 1 -o t_tmp/hypno.db -s HYPNO
destrat t_tmp/hypno.db +HYPNO -v CONF
```
No conflicts, proceed to submit batch job`


Now submit LSF job for all files,
```
/data/nsrr/bin/runner.sh 20 g_s.lst g_param cmd/hypno.txt o g_out/hypno g_out/hypno
destrat g_out/hypno.batch00*.db +HYPNO -v CONF
```
Output looks fine, no conflicts

```
/data/nsrr/bin/runner.sh 20 t_s.lst t_param cmd/hypno.txt o t_out/hypno t_out/hypno
destrat t_out/hypno.batch00*.db +HYPNO -v CONF
```
Output looks fine, no conflicts


## Canonical EDF's

With the conversion of Original EDF's to Canonical EDF's, we will be performing the following:
  - Regenerate Canonical Signals as listed in sigs.canonical file with 'cs' prefix
  - Resample Canonical Signals as listed in sigs.canonical file
  - Drop Signals that are not listed in sigs.canonical file

Submit LSF jobs with the following input,

TWIN
```
cat t_canonical.txt
```

```
ANON
CANONICAL file=nsrr/studies/wsc/sigs.canonical group=WSC
SIGNALS keep=csEEG,csEMG,csLOC,csROC,csECG,csCAN,csTRM,csTHX,csABD,csOXY,csPOS
MINMAX sig=csEEG,csLOC,csROC
WRITE edf-dir=t_canonical/ edf-tag=canonical with-annots sample-list=t_canonical/canonical.lst
```
Run the job
```
/data/nsrr/bin/runner.sh 40 t_s.lst sigs.alias cmd/t_canonical.txt o t_out/canonical t_out/canonical
```

GAMMA

```
cat g_canonical.txt
```

```
ANON
CANONICAL file=nsrr/studies/wsc/sigs.canonical group=WSC
SIGNALS keep=csEEG,csEMG,csLOC,csROC,csECG,csCAN,csTRM,csTHX,csABD,csOXY,csPOS
MINMAX sig=csEEG,csLOC,csROC
WRITE edf-dir=g_canonical/ edf-tag=canonical with-annots sample-list=g_canonical/canonical.lst
```
Run the job
```
/data/nsrr/bin/runner.sh 40 g_s.lst sigs.alias cmd/g_canonical.txt o g_out/canonical g_out/canonical
```

All the canonical EDF's and sample lists are available in t_canonical and g_canonical folders

## Notes:
- In GAMMA, there are occurences of "Central Hypopnea" and "Onstructive Hypopnea" in file "wsc-visit2-33842sco.txt" and "wsc-visit1-59707sco.txt" respectively. These occurences are actually meant to be "Hypopnea". This is taken care of in the annotation file conversion script.
- PLM and PLME events are present in GAMMA study. But WSC team has not scored any PLM events so we have tagged these events under "misc" category in annotation file conversion script.


## Summary
In overall, we have completed the following:
- [x] Input file check (consistency across all original scoring files)
 -- File count with specific extension(s)
 -- Filename: Pattern 
 -- Filename: Number of Columns
 -- Delimiter
- [x] Basic channel/annotation label harmonization 
- [x] checks of sample rates and EDF duration
- [x] Sleep Macro Architecture
- [x] Canonicals generation
- [ ] Post Original EDFs, Canonical EDF's, annotations, excludes, annots/stages list, Channels and the README on NSRR

