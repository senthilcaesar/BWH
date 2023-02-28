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
