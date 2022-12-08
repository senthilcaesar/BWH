# Suggested conventions for upload to NSRR/NAP

The NSRR can accept data in a variety of formats and conventions.
However, adhering to as many of the following guidelines as possible
will faciliate processing via [NAP](http://remnrem.net/nap/)
and reduce the chance of error.

This evolving document covers:

- [Physiologic signal data](#physiologic-signal-data) formats and conventions

- [Annotation](#annotations) formats and conventions

- (Under-development) channel/annotation mapping [webtool](#channel-and-annotation-mapping-tool)

- Participant/recording [identifiers](#identifiers)

- File and folder [naming conventions](#file-and-folder-naming-conventions)


## Physiologic signal data

### Formats

- `EDF` is the standard format for NSRR and NAP, and so is strongly preferred

- `EDF+` is also accepted; however, although EDF+ annotations are extracted
  and used, currently NAP does not automatically handle _discontinuous_
  recordings (i.e. `EDF+D`)

- Importantly: __only open (i.e. non-proprietary) formats are supported)__

- ASCII text files (or R/Matlab matrices/tables) can be accepted,
  although EDF is strongly preferred

- Other open formats (e.g. `BDF(+)`) may be added in future, as needed

- If using Luna to pre-process EDFs, we accept the gzip compressed
  Luna `.edfz` format.  In general, EDFs can and should be compressed
  before being transferred (i.e. simple ZIP or gzip-compressed tar
  archive).

### Specifications/pre-processing

 - Sample rates should be kept as high as possible, e.g. for EEG, EMG, ECG, etc ideally at least 200 Hz
 
 - Do not apply any additional filtering when exporting files: i.e. output raw signals

### Labels

 - Referential signals can be either with respect to the recording reference
 (e.g. `Fpz`, etc) or re-referenced (e.g. `C3-M2`).  The former should
 be labelled `C3-REF` or e.g. `C3-FpZ`, for example. If re-referenced, the
 label should reflect this with the `-` character (i.e. `C3-M2`).

 - for standard PSG assuming a linked or contralateral mastoid
 reference, include `M1-REF` and `M2-REF` (or `A1`, `A2`) channels
 explicitly, i.e. if the channels are not already re-referenced.

 - channels do _not_ need to be in the same order across different
   EDFs, or have the same sample rates; likewise, not all channels
   need be present in all EDFs transferred

 - avoid EDFs with identical channel names; NAP will make these names
 unique (by adding `.1`, `.2`) etc but this introduces an level of
 ambiguity better avoided.   If the EDF does contain channels with the
 same label, please indicate clearly how they differ and give suggested names;
 if they differ by transducer type in the EDF header, please state that, etc.
  
 - (minor) although this goes against the EDF specification, avoid
   spaces in channel names (i.e. `EEG C3-A2`).  (NAP will swap spaces
   with underscores automatically; however, the additional verbosity
   makes working with script/command-line tools unnecessarily
   awkward).  Specifications of what channels are should be part of
   the meta-data that accompanies a study (see below) rather than
   encoded in the 16 characters of an EDF channel.



## Annotations

Unlike EDF, there is no equivalent standard for annotation data.  (EDF+ Annotations are
difficult to work with - e.g. constrained by EDF reecord size, and hard to modify or append
to without re-writing the entire EDF, and so should be avoided.)
 

### Format

Bottom line is that we strongly prefer any ASCII plain-text regular
format, i.e. that is easy to parse automatically by machine.

 - internally, we are using the simple `.annot` format used by Luna,
   as described
   [here](http://zzz.bwh.harvard.edu/luna/ref/annotations/#annot-files)

 - the `.annot` format is just a tab-separated value (`.tsv`) file (although fields can be whitespace delimited too).  As such, it is typically easy for NSRR to convert and text-based format to `.annot`.  As such, it is not necessary to export annotations (for both staging and other events) as `.annot`, so long as a standard (and ideally plain-text) format is used.  That is, `.csv` is preferred over `.xlsx`, etc. 

 - the `.annot` format is quite generic in whether times are specified as epochs, elapsed seconds, or in hh::mm::ss format, also whether intervals only have a start, a start and explicit stop time, or a start time and a duration; also, it is flexible with respect to including any associated meta-data (e.g. the extent of desaturation)

### Labels

Any set of consistent labels can be used to identify events; please
send meta-data that gives full descriptions of what any abbreviations
mean. Part of NAP is to map any annotations (and variants thereof) to
a standard, harmonized set of annotation labels.   These are currently being
reviewed, by a page will be made public to specify these soon.

Ensure that no PHI is present in the annotations (or EDF headers).

### Staging

 - Ideally staging should align with the start of the EDF recording, although this is not strictly necessary. 

 - When possible, include explicit `LightsOn` and `LightsOff` annotations (or use the `L` epoch code) 

 - Annotate any breaks that occurred in the recording 

 - Provide a legend if using non-standard terms (e.g. numeric values 0,1,2,3,etc...) 


## Channel and annotation mapping tool

_Note: this section is a place-holder - please re-visit in a few weeks..._

We have an _under-development_ and (currently undocumented) web tool that can be used to indicate
how a set of channel labels and annotations map onto the NSRR/NAP
harmonized sets:
[http://18.188.74.28/cgi-bin/cgi-mapper](http://18.188.74.28/cgi-bin/cgi-mapper)

Here is a [pre-populated](http://18.188.74.28/cgi-bin/cgi-mapper?f1=C3^C4^M1^M2^LOC^ROC^ECG2^ECG1^EMG1^EMG2^EMG3^L_Leg^R_Leg^AIRFLOW^THOR_EFFORT^ABDO_EFFORT^SNORE^SUM^POSITION^OX_STATUS^PULSE^SpO2^NASAL_PRES^PlethWV^Light^HRate&f2=ASDA_arousal|Arousal_(ASDA)^Hypopnea|Hypopnea^Limb_movement_-_left|Limb_Movement_(Left)^Limb_movement_-_right|Limb_Movement_(Right)^REM_sleep|5^SpO2_artifact|SpO2_artifact^SpO2_desaturation|SpO2_desaturation^Stage_1_sleep|1^Stage_2_sleep|2^Stage_3_sleep|3^Wake|0&f3=&f4=&inp=1) version version (i.e. click submit) to test.


## PSG meta-data

In addition to the general NSRR meta-data forms associated with uploading data: specfic to the EDF/annotation data: 

 - understanding that variations between individual recordings are to be expected, it is still useful to have
  an overview of the canonical/typical montage/sample rate etc used

 - as noted above, descriptions of the terms used in annotation files
   are important; also, noting which annotations are key, and which
   (if any) are ignorable

 - contact details for a person who can respond to any technical
   questions about the recordings


## Identifiers

Identifiers for individuals/studies should follow these conventions where possible:

- no [personal identifiable information (PII)](https://www.dol.gov/general/ppii)

- no special characters (e.g. `$`, `(`, `^`, `/` or `\`, etc)

- no whitespace characters in IDs

- use only alphanumeric, period, underscore and hyphen/minus characters

- start IDs with a letter rather than a digit code; certainly avoid
  starting with leading zeros (i.e. `0001`)

- be consistent with case: `ID1` is _not_ the same as `id1` 

- a good strategy is for IDs to take the form
 `cohort-individual-recording`, e.g. `shhs-00021-2` would denote the
 study (SHHS), the individual (`00021`) and the night/recording (`2`)
 in a way that is easy to parse by machine and by eye.

- if using the `-` character to denote these three fields, ensure that
  no individual component has a `-` also (i.e. there should always be
  exactly three `-`-delimited fields.  If IDs are existing (e.g. for
  individuals) and contain `-`, then chooise `_` or `.` as the
  delimiting character: e.g. `shhs-study.00021.2`

- in general, try to use a consistent scheme across all individuals in
  the study

- if individual IDs are numeric, then add leading zeros (as long as
  the entire ID does not start with leading zeros): i.e. in the above
  example `shhs-00021-1` is preferred over `shhs-21-1`.  The reason
  for this is that sorting IDs will preserve the numeric ordering if
  leading zeros are present (i.e. otherwise `shhs-11-1` comes before
  `shhs-2-1` in a lexicographic (dictionary order) sort.

- if the same subject has a repeated study, new codes can be added:
  e.g. `shhs-00021-1rep` to denote a repeated version of that first
  study; maintaining the simple `study-individual-recording` structure
  makes it easy to extract individual level IDs (i.e. `shhs-00021`) to
  link to time-invariant, individual level variables.


## File and folder naming conventions

Ideally, name files to match the IDs of the corresponding
study/individual/recording, with a single extension speciyfing the
type of file. This makes it easier to match up EDFs and
scoring/annotation files automatically.

For example, for individual/recording `shhs-00021-1` the EDF and annotation files will be:

 - EDF `shhs-00021-1.edf`

 - Annotations: `shhs-00021-1.txt`, `shhs-00021-1.xml` or `shhs-00021-1.tsv`

If there are multiple versions of the one type of file, it is better
to use different folders to denote that: e.g. if scoring from two
manual stagers were available:

 - advised: `scorer1/shhs-00021-1.tsv` and `scorer2/shhs-00021-1.tsv`
   where `scorer1` and `scorer2` are two folders/directories (note: in
   Windows, these would be `scorer1\` and `scorer2\`)

 - not advised: e.g. `shhs-00021-1_S1.tsv` and `shhs-00021-1_S2.tsv` -
   as this breaks the convention of being able to match EDFs and
   annotations unambiguously based on the `{ID}.{extension}` rule

Similarly, rather than two scorers, the file structure might track
separately staging versus respiratory annotations
(e.g. `staging/shhs-00021-1.tsv` and `resp/shhs-00021-1.tsv`).  Folder
names can be anything, but ideally will not contain spaces or special
characters.


Aside from the convention to use folder structure to denote different "versions" of the same type of file, there are no specifications on folder structure.  For example, all files could be in one folder, or individual folders, or some combination.    For example, all the below are valid (i.e. as the all preserve the {ID}.{extension} rule): for two individuals:

Single folder:
```
data/shhs-00021-1.edf
data/shhs-00021-1.tsv

data/shhs-00022-1.edf
data/shhs-00022-1.tsv
```

Split by file type:
```
data/edfs/shhs-00021-1.edf
data/edfs/shhs-00022-1.edf

data/annots/shhs-00021-1.tsv
data/annots/shhs-00022-1.tsv
```

Split by file type (with two types of annotations):
```
data/edfs/shhs-00021-1.edf
data/edfs/shhs-00022-1.edf

data/annots/staging/shhs-00021-1.tsv
data/annots/staging/shhs-00022-1.tsv

data/annots/resp/shhs-00021-1.tsv
data/annots/resp/shhs-00022-1.tsv
```

One folder = one individual (nb. here with repeated measures for both individuals)
```
data/shhs-00021/shhs-00021-1.edf
data/shhs-00021/shhs-00021-1.tsv
data/shhs-00021/shhs-00021-2.edf
data/shhs-00021/shhs-00021-2.tsv

data/shhs-00022/shhs-00022-1.edf
data/shhs-00022/shhs-00022-1.tsv
data/shhs-00022/shhs-00022-2.edf
data/shhs-00022/shhs-00022-2.tsv
```

One folder = one recording: 
```
data/shhs-00021-1/shhs-00021-1.edf
data/shhs-00021-1/shhs-00021-1.tsv
data/shhs-00021-2/shhs-00021-2.edf
data/shhs-00021-2/shhs-00021-2.tsv

data/shhs-00022-1/shhs-00022-1.edf
data/shhs-00022-1/shhs-00022-1.tsv
data/shhs-00022-2/shhs-00022-2.edf
data/shhs-00022-2/shhs-00022-2.tsv
```

The specific folder names (`data`, `edfs`, `annots`, etc) are
arbitrary.  However, if organizing the data by folders corresponding
to individuals or recordings, then the folder name should match the
implied individual / recording ID, i.e.  if using the above identifier
convention, this will be `shhs-00021/` (for the individual), or
`shhs-00021-1` (for the particular recording), etc.


