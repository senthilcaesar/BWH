# Harmonizing NSRR polysomnography data

The NSRR provides data from a range of cohorts, collected by a variety of research teams in different populations at different times and settings, often using different recording devices and data labelling/export protocols.  To facilitate cross-cohort analysis, we have established a series of checks and data cleaning protocols with the goal of harmonizing each cohort's physiologic signal and annotation data.  As many NSRR cohorts were themselves collected as multi-site studies, harmonization will also make individual cohorts easier to work with, by adopting a consistent nomenclature and flagging potential issues.  

Moving forwards, new cohorts will be released in two forms: 
 - an original dataset presented _as is_, 
 - and a corresponding NSRR _harmonized_ version.  

For most users, we advise working with the _harmonized_ dataset if available.  As well as harmonizing newly released cohorts, we also aim to retrospectively harmonize all existing major NSRR polysomnography datasets over the coming months.

## What are the goals of harmonization?  Why is it needed?

Combining different datasets can be logistically and technically challenging, in part due to the overhead of dealing with variable standards between (and within)  cohorts.  At best, this can be tedious; at worst, it can increase chances of error or misinterpretation.  Many of the differences between studies are only superficial, however, and can in principle be avoided.  For example, NSRR cohorts contain many redundant, interchangeable terms for the same basic signals: across six common cohorts (SHHS, MrOS, SOF, CFS, CHAT and CCSHS) are there 328 distinct channel (EDF) labels and 79 annotation labels.  However, many terms are only trivially different alternatives, for example:

 - at least six different labels for abdominal effort channels (`ABD`,  `Abdo`,  `ABDO EFFORT`,  `Abdominal` & `ABDO RES`)
 - inconsistent use of spaces and capitalization (`RIGHT LEG1`, `Rleg`, `R Leg`, `RLEG`, `Rleg1`, `RLeg1`)
 - diverse labels for the same concept: e.g. `Arousal ()`,  `Arousal|Arousal ()`  `Arousal|Arousal (Arousal)`, `Arousal|Arousal (Standard)`,  `Arousal|Arousal (STANDARD)`, `Arousal(Asda)`, `Arousal (ASDA)`, etc

Although such inconsistencies do not impede manual inspection or scoring, they can be more problematic when developing automated pipelines. Likewise, different recordings may contain the same fundamental type of information, for instance, for a left-central EEG referenced to the contralateral mastoid, but with the information represented differently:

| | Record 1 | Record 2 | Record 3 |
|---|---|---|---|
| Channel label | `EEG C3-M2` | `C3` | `EEG2` |
| Reference | none (already referenced) | `A2` | `REF` |
| Sample rate | 256 Hz | 128 Hz | 150 Hz |
| Units | uV | mV | V

Although these three signals can of course be relabelled, re-referenced, rescaled or resampled to yield similar representations easily enough, performing such steps across hundreds or thousands of individuals, when multiple variations of the above exist within and between cohorts and potentially impact many channels, can quickly become burdersome.

One key goal of NSRR harmonization, therefore, is to remap channel and annotation nomenclature to a single, common set. Beyond label harmonization, there are several other components of the harmonization process that aim to deliver more robust, analysis-ready recordings:  

 - __Core sanity checks on all files:__ signal files (European Data Format, EDF) are checked for corrupted header information or file truncation

 - __Fix units, referencing, sample rates, EDF record sizes:__ when transformations do not meaningfully change the core information content of the signal, we adopt a single convention for representing common signal types: e.g. following above EEG example, central EEGs will be re-referenced to the contralateral mastoid, labelled `C3_M2`, scaled to `uV` units and resampled to 128 Hz.  Labels such as `C3_M2` are chosen with ease of downstream data manipulation in mind, e.g. if using command-line packages such as [Luna](https://zzz.bwh.harvard.edu/luna/) or programming tools such as Matlab or R, it is easier to specify variables based on channel labels that do not contain special characters (e.g. arithmetic symbols) or spaces.
  
 - __Dropping undocumented, flat, corrupted or duplicated signals:__ if signals are obviously lacking valid physiologic information, we drop them from the EDF to reduce file sizes and streamline subsequent data processing

 - __Check for protected health information (PHI):__  following manual review, PHI (in filenames, EDF headers or annotations) is removed 

 - __Review signal distribution properties:__ for key raw signals and selected derived metrics (e.g. stage durations, spectral power), we generate descriptive statistics (raw signal means, SD, range), visualizations (e.g. histograms for key metrics, EEG power spectra plotted both per individual and at a group level) and where appropriate _true positive_ control analyses (e.g. ensuring an expected age-related trend is seen if one would be expected).  For key derived we check for outliers and any evidence of extreme distributions (e.g. bimodalities that may indicate some subset of the data has been mislabelled, e.g. labelling _V_ instead of _uV_, etc.  For multi-site studies, we additionally look for evidence of large inter-site differences, conditional on key demographics, that may be indicative of site-specific technical issues. 

 - __Subsetting or "flattening" gapped (discontinous) recordings:__ although most NSRR recordings are continuous, we also receive "discontinous" recordings containing _gaps_, either as EDF+D files or multiple standard EDF files per night.  In one instance, each night's recording was split across hundreds of 5-minute EDFs; in a second instance, single nights were exported to several partially overlapping EDFs, which sometimes contained gaps.  For many purposes, it is easier to work with simple, continuous signals, however.  Within [Luna](https://zzz.bwh.harvard.edu/luna), we have implemented routines for merging multiple EDFs as well as reducing gapped EDFs to continuous recordings, via the [EDF-MINUS](https://zzz.bwh.harvard.edu/luna/ref/manipulations/#edf-minus) command. 


 - __Standard formats for annotation data:__   Despite some inherent limitations of the EDF format itself, and a lack of adherence by many EDF-generating software tools to the precise EDF specification, the existence of a widely used community standard (i.e. EDF) is incredibly beneficial.  There is no equivalent for _annotation data_, however.  Here _annotation data_ means any type of interval- or event-based descriptors attached to signals, including sleep staging and respiratory events.  Although EDF+ nomimally supports annotations, in practice it is far from optimal.  Annotation data are typically received in _ad hoc_ in-house formats (e.g. CSV or XML but with no universal schema/column structure) and may use different conventions for representing dates and time (clock-time versus offsets in elapsed seconds, event durations versus stop times, 12 versus 24 hour clocks, etc).   We have developed a simple, machine-readble text-based format and [tools](https://zzz.bwh.harvard.edu/luna/ref/annotations/) to read and write it. Importantly, the format is designed 1) to allow flexibility in how events are specified, 2) to be easy to generate by hand, by spreadsheet or programmatically, and 3) is easy to edit, visualize and analyze, using a regular format that can be read by common tools such as R, Excel, etc as well as command-line text-processing tools such as `awk`, `grep`, etc.  This single `.annot` format will be used across all harmonized datasets.

 - __Alignment of staging and signals:__ Given annotation data in a common file format, one step of the harmonization process is to check the integrity of any (manual) sleep staging: e.g. 1) mapping labels to a standard set (`N1`, `N2`, `N3`, `R` and `W` as well as `?` and `L`) and checking for atypical configurations (e.g. reduced staging with only `N2` and `W` epoch observed), 2) checking that most of the EDF recording interval has overlapping associated staging information, and vice versa, 3) using a statistical approach ([SOAP](https://zzz.bwh.harvard.edu/ref/soap)) to check the consistency of staging and signals (e.g. EEG), and potentially detect and correct temporal misalignment (e.g. if staging and signals were accidentically shifted _N_ epochs relative to each other) and finally, 4) ensuring that EDF records, epochs and staging annotation boundaries are temporally aligned, which makes downstream analysis and manipulation of EDFs easier.

__Collectively, these steps seek to make harmonized datasets more robust and easier to work with, without fundamentally altering the core information content of recordings.__  

When thinking about data cleaning and harmonization broadly, there are of course numerous additional procedures that one could potentially employ.  For example, considering only polysomnography EEGs:

  - signals may vary in terms of prior (hardware) filtering, which may or may not have been documented (either at the study level, or in the pre-filtering field of each EDF header)
  - signals may be mislabelled (e.g. swapping O1 and Fp1 channels)
  - signals may have flipped polarity (e.g. as may be evident by [considering](https://zzz.bwh.harvard.edu/luna/vignettes/nsrr-polarity) phase-based analysis or the waveforms of NREM transients)
  - signals may have marked electrical line noise or cardiac contamination
  - etc.  

Although it may in some cases be possible to detect and even (largely) correct such issues (all of which have been observed in incoming NSRR recordings), we generally __do not__ attempt this secondary level of pre-processing.  Oftentimes, issues may be ambiguous to detect (e.g. whether or not polarity is flipped) meaning that any fix may only be applied to a subset of records.  Furthermore, in some cases a "fix" may actually introduce new and potentially trickier problems.  Consequently, in most cases we intend only to _flag_ such issues as and when they arise, as the optimal solution will likely depend on the exact nature of the planned subsequent analysis.   These decisions may vary by cohort and context however: if we do perform some of these secondary modifications of the data, this will be explicitly noted in the harmonization notes.  


__In other words: NSRR harmonization aims to make files easier to work with by removing unncessary differences between studies where possible.  This does _not_ mean that harmonized datasets can always be assumed to be directly comparable: inherent differences in populations, study designs and recording protocols still need to be factored into any analytic plan that considers multiple cohorts.__


## Why still distribute original _as is_ datasets?

Harmonization adds value to the original data, through the identification of gross problems as well as improved consistency in labelling.  So why would we still want to distribute the original data if the harmonized set is better?   There are several good reasons for this:

 - part of the harmonization process is to drop unusual, uncommon (e.g. only seen in a few individuals), undocumented, or potentially redundant channels, annotations and/or recordings.   There may be cases where a user would rather include these and work through any issues themselves.
 - naturally, we can never completely guarantee that the harmonization process might not introduce some secondary issues in the data
 - to support reproducible and transparent research, keeping the original data available can a) help in interpreting any pre-existing published analyses and b) allow users to adopt their own harmonization procedures, as well as c) make the NSRR harmonization process more transparent.


## An example: the Pediatric Adenotonsillectomy Trial of Snoring (PATS) study

The PATS cohort is a dataset [recently available dataset on NSRR](https://sleepdata.org/datasets/pats).  PATS was a multi-centered, randomized, single-blinded, 12-month interventional study that compares the impact of adenotonsillectomy (AT) on measures of behavior, quality of life, and healthcare utilization in children aged 3-12 with mild sleep-disordered breathing (SDB), conducted from 2016 to 2022.  Here we step through some of the checks and observations made throughout the PATS harmonization process:

### EDF Checks

Out of 1,010 total EDFs, 749 were standard continuous EDFs, while 261 were EDF+C files (uninterrupted recordings but with EDF+ annotations).  All files adhered to a standard record size of 1 second and most maintained typical sample rates.  A few channels with unusual sample rates of 258 Hz, 276 Hz, and 287 Hz were noted.  One EDF file (`pats-813652-baseline`) was corrupt and therefore removed.

### Protected Health Information (PHI)

Reviewing EDF headers, 261 files contained technician-written annotations, which were scanned for potential PHI. After determining that the EDF+ annotations did not contain critical information (beyond the staging and event-annotation available in the corresponding text annotation files), we used the <a href="https://zzz.bwh.harvard.edu/luna/" target="_blank">Luna </a>commands to drop the EDF+ annotations, and export the files as standard EDFs (thereby also ensuring that no PHI data were leaked).  Additionally, all EDF header date values were anonymized, set to the default date of January 1, 1985 as specified by the EDF standards. 

### Signal Review and harmonization

The Harmonized EDFs contains the following list of channel labels.

| Harmonized Channel labels         |
|--------------|
| C3_M2        |
| C4_M1        |
| F3_M2        |
| F4_M1        |
| O1_M2        |
| O2_M1        |
| HR           |
| LOC          |
| ROC          |
| LAT          |
| RAT          |
| ECG          |
| EMG          |
| cap          |
| SpO2         |
| thorax       |
| EtCO2        |
| pleth        |
| abdomen      |
| airflow      |
| nasal_pres   |
| thermistor   |


- **EEG selection:** 

   All individuals had EEG channels present featuring central (`C3, C4`) frontal (`F3, F4`) and Occipital (`O1, O2`) 
   each referenced to the contralateral mastoid `M1` and `M2`.

- **EMG Referencing:** 

   Individuals had either one, two or three EMG channels, with some 
   referenced to the recording reference, but others with bipolar EMG derivations. 
   For 782 individuals, we re-referenced EMGs, whereas 227 already meet the standard
   referencing protocol. For consistency, re-referenced EMGs
   were uniformly renamed: `LChin_RChin` to `EMG1`, `LChin_CChin` to
   `EMG2`, and `RChin_CChin` to `EMG3`.

- **ECG Signal Correction:**

   Issues with the polarity of ECG signals in some subjects were
   addressed by flipping the signals back to their correct
   orientation. Further, signals labeled `ECG3`, which are the actual
   ECG readings, are renamed to `ECG` for uniformity.

- **Position Channel Review:**

   Among 1010 individuals, approximately 200 had position
   channels.  However, after reviewing the signal distribution, these data were
   considered unreliable due to lack of calibration, and so dropped from the harmonized set.
   
### Channel naming QC

We observed a moderate degree of variability in channel naming conventions across PATS sites, 
particularly for oxygen saturation, airflow and CO<sub>2</sub> levels. 
This inconsistency includes variations in case sensitivity and
additional characters such as underscores or suffixes.

 - oxygen saturation channels included: `SpO2`, `SpO2x`, `SPO2__2`, `SpO2xx`, `SPO2`
 - airflow channels included: `Airflow`, `AirFlow`, `AirflowXX`, `Airflowx`, `Airflow_x`, `Airflow2`
 - carbon dioxide levels included: `EtCO2`, `ETCO2`, `EtCO2_2`, `NK_EtCO2`, `EtCO2___NK`, `EtCO2_XT`, `EtCO2_Neo`, `NK_EtCO2_Wave`, `EtCO2__NK`

Visualizing channel distributions, it was evident that a subset of these channels were simply empty or corrupt channels. Valid signals selected for inclusion in the harmonized sets were `SpO2`, `Airflow`, `AirFlow`, `EtCO2`, and `EtCO2_2`.


### Annotations

Annotation files <a
href="https://zzz.bwh.harvard.edu/luna/ref/annotations/#luna-annotations"
target="_blank">(.annot) </a> are generated from the original XML
files. The event `Periodic Breathing` has been removed. The event
labeled as `Unsure` is, in PATS, mapped to `Mixed_Apnea`.

Here are some of the most common events and their mappings:

| Original label       | Harmonized _class_ label | Harmonized _instance_ label |
|----------------------|------------------|---------------------|
| N1                   | N1               | .                   |
| N2                   | N2               | .                   |
| N3                   | N3               | .                   |
| W                    | W                | .                   |
| R                    | R                | .                   |
| EtCO2_artifact       | artifact         | EtCO2               |
| SpO2_artifact        | artifact         | SpO2                |
| SpO2_desaturation    | desat            | .                   |
| Arousal              | arousal          | .                   |
| Arousal_ASDA         | arousal          | .                   |
| Limb_Movement_Right  | LM               | right               |
| Limb_Movement_Left   | LM               | left                |
| PLM_Right            | PLM              | right               |
| PLM_Left             | PLM              | left                |
| Hypopnea             | hypopnea         | .                   |
| Obstructive_Apnea    | apnea            | obstructive         |
| Central_Apnea        | apnea            | central             |
| Unsure               | apnea            | mixed               |


### Staging QC

Luna's <a
href="https://zzz.bwh.harvard.edu/luna/ref/hypnograms/#stage"
target="_blank">STAGE </a> command was used to relabel sleep stages as 
`W` (Wake), `N1`, `N2`, `N3`, `R` (REM sleep), `?`, and `L`. This step
also confirmed that there are no _conflicts_ in sleep
staging (i.e. inconsistent overlapping stage annotations). Additionally, the <a
href="https://zzz.bwh.harvard.edu/luna/ref/annotations/#spanning"
target="_blank">SPANNING </a> command checked that the staging
annotations cover almost the entire duration of the recordings with no
EDFs missing more than 1% of coverage. The <a
href="https://zzz.bwh.harvard.edu/luna/ref/hypnograms/#hypno"
target="_blank">HYPNO </a> command was used to ensure that all studies
contain at least some recorded sleep stages; one individual
was found to have less than one hour of recorded sleep. Generally,
most individuals were assigned at least four different sleep stage
labels.  Overall, no major issues were detected with the existing 
manual staging in PATS. 




### Secondary data modifications

Finally, based on visual review of the data, we performed two _second-level_ harmonization steps in PATS:

 - to detect and flip likely EEG polarity flips, which appeared to impact a large proportion of the sample
 
 - to use empirical averaged ECG waveforms to select a single, consistently oriented ECG per individual, given the range of ECG montages in the original data


#### ECG selection

```
TOOD: give sense of a) the original issue, b) the procedure, c) show some graphics to illustrate, d) generate and point to a clear entry in an ISSUES file that shows what was selected for whom
```

#### EEG polarity Check

Quality review revealed that there were 289 individuals whose central EEG channels (`C3_M2` and `C4_M1`) required correction
due to polarity issues.   

The corrective action involves identifying these specific EEG channels for the affected individuals, flipping the polarity, and then saving the corrected files in a format suitable for distribution by the National Sleep Research Resource (NSRR). To ensure
the accuracy of this process, a luna <a href="https://zzz.bwh.harvard.edu/luna/ref/manipulations/#flip" target="_blank">FLIP </a> command was executed on all 289 individual's data. Upon visual inspection, we confirmed that the polarity of the channel has indeed been flipped as intended.

## Listing known issues in harmonized datasets

```
TBD
```

## Summary

```
TODO: rather than abruptly finish, close w/ a brief summary and notes and future directions, etc
```
 