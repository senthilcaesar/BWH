### ERIS Folder structure

| Folder | Description |
|----|----|
| `/data/nsrr/working/msp-edfs/` | Original deposited data |
| `/data/nsrr/working/msp/` | Working directory for processing |

We also create the following key folders in the working directory:

| Folder name | Contents |
|----|----|
| `files/` | Other relevant, original non-data files (e.g. descriptive PDFs, meta-data as XLS/PPTs, etc) |
| `sl/` | All derived sample lists and other ID lists (e.g. ID exclusions) | 
| `tmp/` | Scratch folder, e.g. including all raw outputs from Luna runs (`.db` files prior to extraction to `res/` etc) |
| `derived/` | Any intermediate EDF/annotation files (i.e. steps between original import and running NAP) |  
| `cmd/` | Any Luna scripts |
| `res/` | Final compiled/merged/derived _higher-level_ outputs, e.g. typically compiled from `tmp/` |  
| `nap/` | NAP automatically generates this folder with NAP output and harmonized EDFs/annotations
| `dist/` | Final distribution folder, i.e. copied EDFs/annotations from `nap/` that will be imported into NSRR |

```
cd /data/nsrr/working/msp
mkdir files sl tmp derived res nap dist cmd
```

## 1) Key design factors


There are 106 individuals, each with one EDF:

```
ls ../msp-edfs/*.edf | wc -l
106
```
Key demographic data (including site) are in
`files/msp-dataset-0.1.0.pre-racecat.csv`. Here we'll
make a version with temporary IDs used to process files here, with
`id_` prefixes, and save a reformatted tab-delimited file.
```R
d <- read.csv("files/msp-dataset-0.1.0.pre-racecat.csv")
d$ID <- gsub( "^", "id_",  d$id ) 
d <- d[ , c( "ID" , "mat_age" , "inf_sex" , "mat_race"  ) ]

# check for missing data
table( complete.cases( d ) )

# reformat
names(d) <- c("ID","maternalage","infantsex","maternalrace") 
d$infantsex <- ifelse( d$infantsex == 1 , "M" , "F" )
write.table( d , file="files/demo.txt" , sep="\t" , row.names=F, quote=F, col.names=T )
```
Note that `ID` numbers are not sequential due to excluded participants
```
head files/demo.txt

ID	maternalage	infantsex	maternalrace
id_50	  32	        F	          White
id_53	  20	        F	          White
id_55	  23	        F	          White
id_75	  33	        F	          White
id_94	  30	        F	          White
id_125	  28	        F	          White
id_127	  29	        M	          White
id_129	  39	        M	          White
id_1	  25	        M	          African American
```
Infant sex information is missing for `id_60`
```
id_60	  27	        NA	          African American
```

Each individual has one csv and one xml file, e.g.:
```
ls ../msp-edfs/ | head
```
```
S001.edf
S001_event.csv
S001_Import.xml
S002.edf
S002_event.csv
S002_Import.xml
S003.edf
S003_event.csv
S003_Import.xml
```
```
ls ../msp-edfs/*_event.csv | wc -l
106

ls ../msp-edfs/*_Import.xml | wc -l
106
```

Other pertinent points include:

  - The raw data is scored and annotated by a Registered PSG Technologist (RPSGT); Exported annotations by RemLogic into XML files
  - There are 3 Maternal race/ethnicity (self-reported) within MSP (labels `White` `African American` `Asian/Indian`
    that correspond to IDs with the digit 1, 2, 3 respectively)
  - Equipment and device used to capture the raw polysomnography data `Natus RemLogic N7000`
  - | Sleep Stage |  |
    |-----|-----|
    | `Wake` | 0 |
    | `N1` |	1	|
    | `N2` |	2	|
    | `N3` |	3	|
    | `N4` |	4	|
    | `REM` |	5	|

## 2) EDF checks

We first build an EDF-only sample list (`sl/s0.lst`)

```
luna --build ../msp-edfs/ > sl/s0.lst
```
```
wrote 106 EDFs to the sample list
  106 of which had 0 linked annotation files
```
Checking the validity of all EDFs:
```
luna sl/s0.lst -s DESC
```
(the above completes with no issues detected - e.g. no corrupt/truncated EDF files)
Running `HEADERS` to summarize EDF headers:

```
luna sl/s0.lst -o tmp/headers.db -s HEADERS signals
```
All files are standard EDFs (i.e. continuous and without EDF Annotations)

```
destrat tmp/headers.db +HEADERS -v EDF_TYPE  | cut -f2 | sort | uniq -c
```
```
106 EDF
```
The longest studies are around 11 hours:
```
destrat tmp/headers.db +HEADERS -v TOT_DUR_SEC | sort --key=2 -nr | head -5
```
```
id_S022	42000
id_S004	39600
id_S043	38400
id_S024	38400
id_S005	38400
```
There is one short recordings (5 hours)
```
destrat tmp/headers.db +HEADERS -v TOT_DUR_SEC | sort --key=2 -n | awk ' { print $1,$2  } ' | head -6
```

<pre>
ID      TOT_DUR_SEC
<b>id_S099    19200</b>
id_S032   25200
id_S127   26400
id_S044   28800
id_S051   28800
</pre>


