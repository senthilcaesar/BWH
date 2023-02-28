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


There are 1,105 individuals, each with one EDF:

```
ls ../msp-edfs/*.edf | wc -l
106
```
Key demographic data (including site) are in
`files/msp-dataset-0.1.0.pre-racecat.csv`. Here we'll
make a version with temporary IDs used to process files here, with
`id_` prefixes, and save a reformatted tab-delimited file:
```R
d <- read.csv("files/msp-dataset-0.1.0.pre-racecat.csv")
d <- d[order(d$id),]
d$ID <- gsub( "^", "id_",  d$id ) 
d <- d[ , c( "ID" , "mat_age" , "inf_sex" , "mat_race"  ) ]

# check for missing data
table( complete.cases( d ) )

# reformat
names(d) <- c("ID","age","sex","race") 
d$sex <- ifelse( d$sex == 1 , "M" , "F" )
write.table( d , file="files/demo.txt" , sep="\t" , row.names=F, quote=F, col.names=T )
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
